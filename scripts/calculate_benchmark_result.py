import argparse
import os
from collections import defaultdict

import pandas as pd
from posebusters import PoseBusters
from tqdm import tqdm

POSEBUSTER_TEST_COLUMNS = [
    # accuracy #
    "rmsd_≤_2å",
    # chemical validity and consistency #
    "mol_pred_loaded",
    "mol_true_loaded",
    "mol_cond_loaded",
    "sanitization",
    "molecular_formula",
    "molecular_bonds",
    "tetrahedral_chirality",
    "double_bond_stereochemistry",
    # intramolecular validity #
    "bond_lengths",
    "bond_angles",
    "internal_steric_clash",
    "aromatic_ring_flatness",
    "double_bond_flatness",
    "internal_energy",
    # intermolecular validity #
    "minimum_distance_to_protein",
    "minimum_distance_to_organic_cofactors",
    "minimum_distance_to_inorganic_cofactors",
    "volume_overlap_with_protein",
    "volume_overlap_with_organic_cofactors",
    "volume_overlap_with_inorganic_cofactors",
]


def get_group_info(dataset: str, dataset_folder: str) -> pd.DataFrame:
    group_dict = defaultdict(list)
    for item_name in os.listdir(dataset_folder):
        item_dir = os.path.join(dataset_folder, item_name)
        if not os.path.isdir(item_dir):
            continue
        group_dict["PDB_CCD_ID"].append(item_name)
        if dataset == "posex_cross_dock":
            group_path = os.path.join(dataset_folder, item_name, "group_id.txt")
            with open(group_path, "r") as f:
                lines = f.readlines()
                lines = [line.strip() for line in lines]
            group_dict["PDB_GROUP"].append(lines[0])
            group_dict["GROUP"].append(lines[1])
        elif dataset in ["posex_self_dock", "posex_supp"]:
            group_dict["PDB_GROUP"].append(item_name)
            group_dict["GROUP"].append(item_name)
        else:
            raise RuntimeError()
    df_group = pd.DataFrame(group_dict)
    return df_group

def main(args: argparse.Namespace):
    docking_data = pd.read_csv(args.input_file)
    bust_dict = defaultdict(list)

    total_samples = len(docking_data["PDB_CCD_ID"])
    for pdb_ccd_id in tqdm(docking_data["PDB_CCD_ID"]):
        mol_true = os.path.join(args.dataset_folder, f"{pdb_ccd_id}/{pdb_ccd_id}_ligands.sdf")
        if args.model_type == "alphafold3":
            pdb_ccd_id = pdb_ccd_id.lower()
        mol_cond = os.path.join(args.model_output_folder, f"{pdb_ccd_id}/{pdb_ccd_id}_model_protein_aligned.pdb")
        mol_pred = os.path.join(args.model_output_folder, f"{pdb_ccd_id}/{pdb_ccd_id}_model_ligand_aligned.sdf")
        if not os.path.exists(mol_pred):
            print(f"File {mol_pred} does not exist")
            continue
        bust_dict["PDB_CCD_ID"].append(pdb_ccd_id.upper())
        bust_dict["mol_pred"].append(mol_pred)
        bust_dict["mol_true"].append(mol_true)
        bust_dict["mol_cond"].append(mol_cond)

    bust_data = pd.DataFrame(bust_dict)
    print("Number of Benchmark Data: ", total_samples)
    print("Number of Posebusters Data: ", len(bust_data))
    save_folder = os.path.dirname(args.input_file)

    # Calculate posebusters result
    buster = PoseBusters(config="redock", top_n=None, max_workers=None)
    buster.config["loading"]["mol_true"]["load_all"] = True
    bust_results = buster.bust_table(bust_data, full_report=True)
    df_index = bust_data[["mol_pred", "PDB_CCD_ID"]].rename(columns={"mol_pred": "file"})
    bust_results = pd.merge(bust_results.reset_index(level="file"), df_index, on='file').drop(columns="file")
    if args.dataset in ["posex_self_dock", "posex_cross_dock", "posex_supp"]:
        df_group = get_group_info(args.dataset, args.dataset_folder)
        df_group_sim = pd.read_csv(os.path.join(args.dataset_folder, "qtm.csv"))
        bust_results = pd.merge(bust_results, df_group, on="PDB_CCD_ID")
        bust_results = pd.merge(bust_results, df_group_sim, on="GROUP")
        if args.dataset == "posex_cross_dock":
            total_samples = df_group.GROUP.unique().shape[0]
    if args.relax == "true":
        res_path = os.path.join(save_folder, f"{args.dataset}_benchmark_result_{args.model_type}_relax.csv")
    else:
        res_path = os.path.join(save_folder, f"{args.dataset}_benchmark_result_{args.model_type}.csv")
    bust_results.to_csv(res_path, index=False)
    if args.dataset in ["posex_self_dock", "posex_cross_dock", "posex_supp"]:
        test_data = bust_results[POSEBUSTER_TEST_COLUMNS].copy()
        bust_results.loc[:, "pb_valid"] = test_data.iloc[:, 1:].all(axis=1)
        bust_results = bust_results.groupby("PDB_GROUP").agg({"rmsd": "mean", "pb_valid": "mean", "GROUP": "first"})
        bust_results = bust_results.groupby("GROUP").agg({"rmsd": "mean", "pb_valid": "mean"})
        accuracy = len(bust_results[bust_results["rmsd"] <= 2.0]) / total_samples
        print(f"RMSD ≤ 2 Å: {accuracy * 100:.2f}%")
        valid_data = bust_results[(bust_results["rmsd"] <= 2) & (bust_results["pb_valid"] >= 0.5)]
        print(f"RMSD ≤ 2 Å and PB Valid: {len(valid_data) / total_samples * 100:.2f}%")
    else:
        # Calculate accuracy    
        accuracy = len(bust_results[bust_results["rmsd_≤_2å"] == True]) / total_samples
        print(f"RMSD ≤ 2 Å: {accuracy * 100:.2f}%")

        # Calculate posebusters test result
        test_data = bust_results[POSEBUSTER_TEST_COLUMNS].copy()
        test_data.loc[:, "pb_valid"] = test_data.iloc[:, 1:].all(axis=1)
        valid_data = test_data[test_data["rmsd_≤_2å"] & test_data["pb_valid"]]
        print(f"RMSD ≤ 2 Å and PB Valid: {len(valid_data) / total_samples * 100:.2f}%")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", type=str, required=True, help="Path to the benchmark input file")
    parser.add_argument("--dataset_folder", type=str, required=True, help="Path to the dataset folder")
    parser.add_argument("--model_output_folder", type=str, required=True, help="Path to the model output folder")
    parser.add_argument("--dataset", type=str, required=True, help="Dataset name")
    parser.add_argument("--model_type", type=str, required=True, help="Model type")
    parser.add_argument("--relax", type=str, required=True, help="relax mode (true or false)")
    args = parser.parse_args()

    main(args)
