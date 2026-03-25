from collections import defaultdict
from io import StringIO
from math import e
from pathlib import Path

import numpy as np
import openmm.app as app
from openmm.app.element import hydrogen
from openmm.app.internal.pdbx.reader.PdbxReader import PdbxReader
from pdbfixer.pdbfixer import PDBFixer, Sequence, substitutions
from rdkit import Chem

from .pdb_helper import (
    METALS,
    PROTEIN_RESIDUES,
    WATERS,
    add_conformer_to_mol,
    add_h_to_receptor_mol,
    format_4letter,
)

substitutions["OTY"] = "TYR"
substitutions["4AF"] = "ALA"
substitutions["SNN"] = "ASN"
substitutions["HP9"] = "PHE"


def get_sequences_info(cif_fn: Path):
    with open(cif_fn, "r") as file:
        reader = PdbxReader(file)
        data = []
        reader.read(data)
        block = data[0]

    # Load the sequence data.

    sequenceData = block.getObj("pdbx_poly_seq_scheme")
    sequences = {}
    nan_sequences = {}

    if sequenceData is not None:
        asymIdCol= sequenceData.getAttributeIndex("asym_id")
        entityIdCol = sequenceData.getAttributeIndex("entity_id")
        residueCol = sequenceData.getAttributeIndex("mon_id")
        authresCol = sequenceData.getAttributeIndex("auth_mon_id")
        for row in sequenceData.getRowList():
            asymId = row[asymIdCol]
            entityId = row[entityIdCol]
            residue = row[residueCol]
            auth_residue = row[authresCol]
            uniq_key = f"{entityId}.{asymId}"
            if uniq_key not in sequences:
                sequences[uniq_key] = []
                nan_sequences[uniq_key] = []
            sequences[uniq_key].append(residue)
            nan_sequences[uniq_key].append(
                None if auth_residue == "?" else auth_residue
            )

    atomsiteData = block.getObj("atom_site")
    visited_chains = set()
    complete_sequences = []
    uncomplete_sequences = []
    if atomsiteData is not None:
        label_asymIdCol = atomsiteData.getAttributeIndex('label_asym_id')
        auth_asymIdCol = atomsiteData.getAttributeIndex("auth_asym_id")
        entityIdCol = atomsiteData.getAttributeIndex("label_entity_id")
        for row in atomsiteData.getRowList():
            asymId = row[label_asymIdCol]
            entityId = row[entityIdCol]
            uniq_key = f"{entityId}.{asymId}"
            auth_asymId = row[auth_asymIdCol]

            if (uniq_key in sequences) and (asymId not in visited_chains):
                complete_sequences.append(Sequence(auth_asymId, sequences[uniq_key]))
                uncomplete_sequences.append(Sequence(auth_asymId, nan_sequences[uniq_key]))
                visited_chains.add(asymId)

    return complete_sequences, uncomplete_sequences


def trim_none(lst):
    if not lst:  # 空列表直接返回
        return []
    # 正向找第一个非None
    start = 0
    while start < len(lst) and lst[start] is None:
        start += 1
    # 反向找最后一个非None
    end = len(lst) - 1
    while end >= 0 and lst[end] is None:
        end -= 1
    # 切片截取有效区间（注意end+1）
    return lst[start : end + 1] if start <= end else []


class PDBFixerWrapper(PDBFixer):
    def _downloadNonstandardDefinitions(self):
        """If the file contains any nonstandard residues, download their definitions and build
        the information needed to add hydrogens to them.
        """
        app.Modeller._loadStandardHydrogenDefinitions()
        resnames = set(residue.name for residue in self.topology.residues()) - set(
            METALS
        )
        definitions = {}
        for name in resnames:
            if name not in app.Modeller._residueHydrogens:
                # Try to download the definition.
                ccdDefinition = self._downloadCCDDefinition(name)
                if ccdDefinition is None:
                    continue

                # Record the atoms and bonds.
                atoms = [
                    (atom.atomName, atom.symbol.upper(), atom.leaving)
                    for atom in ccdDefinition.atoms
                ]
                bonds = [(bond.atom1, bond.atom2) for bond in ccdDefinition.bonds]
                definitions[name] = (atoms, bonds)
        return definitions

    def addMissingHydrogens(self, pH=7.0, forcefield=None) -> None:
        extraDefinitions = self._downloadNonstandardDefinitions()
        variants = [
            (
                self._describeVariant(res, extraDefinitions)
                if res.name not in ["ACE", "NME"]
                else None
            )
            for res in self.topology.residues()
        ]
        modeller = app.Modeller(self.topology, self.positions)
        modeller.addHydrogens(pH=pH, forcefield=forcefield, variants=variants)
        self.topology = modeller.topology
        self.positions = modeller.positions

    def removeHeterogens(self, keepWater=True, keepMetal=True):
        """"""
        """Remove all heterogens from the structure.

        Parameters
        ----------
        keepWater : bool, optional, default=True
            If True, water molecules will not be removed.

        Returns
        -------
        a list of Residue objects that were removed

        Examples
        --------

        Remove heterogens in Abl structure complexed with imatinib.

        >>> fixer = PDBFixer(pdbid='2F4J')
        >>> fixer.removeHeterogens(keepWater=False)

        """

        keep = set(PROTEIN_RESIDUES)
        if keepWater:
            keep = keep.union(WATERS)
        if keepMetal:
            keep = keep.union(METALS)
        toDelete = []
        for residue in self.topology.residues():
            if residue.name not in keep:
                toDelete.append(residue)
        modeller = app.Modeller(self.topology, self.positions)
        modeller.delete(toDelete)
        self.topology = modeller.topology
        self.positions = modeller.positions
        return toDelete

    def applyMutations(self, mutations, chain_id):
        """Apply a list of amino acid substitutions to make a mutant protein.

        Parameters
        ----------
        mutations : list of strings
            Each string must include the resName (original), index,
            and resName (target).  For example, ALA-133-GLY will mutate
            alanine 133 to glycine.
        chain_id : str
            String based chain ID of the single chain you wish to mutate.

        Notes
        -----

        If a target residue is not a standard amino acid, and if no template
        has been registered for it with registerTemplate(), this function
        attempts to look it up from the PDB and create a new template for it.

        We require three letter codes to avoid possible ambiguitities.
        We can't guarantee that the resulting model is a good one; for
        significant changes in sequence, you should probably be using
        a standalone homology modelling tool.

        Examples
        --------

        Find nonstandard residues.

        >>> fixer = PDBFixer(pdbid='1VII')
        >>> fixer.applyMutations(["ALA:57 :GLY"], "A")
        >>> fixer.findMissingResidues()
        >>> fixer.findMissingAtoms()
        >>> fixer.addMissingAtoms()
        >>> fixer.addMissingHydrogens(7.0)

        """
        # Retrieve all residues that match the specified chain_id.
        # NOTE: Multiple chains may have the same chainid, but must have unique resSeq entries.
        resSeq_to_residue = (
            dict()
        )  # resSeq_to_residue[resid] is the residue in the requested chain corresponding to residue identifier 'resid'
        for chain in self.topology.chains():
            if chain.id == chain_id:
                for residue in chain.residues():
                    resSeq_to_residue[f"{int(residue.id)}{residue.insertionCode}"] = (
                        residue
                    )

        # Make a map of residues to mutate based on requested mutation list.
        residue_map = (
            dict()
        )  # residue_map[residue] is the name of the new residue to mutate to, if a mutation is desired
        for mut_str in mutations:
            old_name, resSeq, new_name = mut_str.split(":")

            if resSeq not in resSeq_to_residue:
                raise (
                    KeyError(
                        "Cannot find chain %s residue %d in system!"
                        % (chain_id, resSeq)
                    )
                )

            residue = resSeq_to_residue[resSeq]  # retrieve the requested residue

            if residue.name != old_name:
                raise (
                    ValueError(
                        "You asked to mutate chain %s residue %d name %s, but that residue is actually %s!"
                        % (chain_id, resSeq, old_name, residue.name)
                    )
                )

            if self._getTemplate(new_name) is None:
                # Try to download a template from the PDB.
                self.downloadTemplate(new_name)
                if self._getTemplate(new_name) is None:
                    raise (
                        KeyError(
                            "Cannot find residue %s in template library!" % new_name
                        )
                    )

            # Store mutation
            residue_map[residue] = new_name

        # If there are mutations to be made, make them.
        if len(residue_map) > 0:
            deleteAtoms = []  # list of atoms to delete

            # Find atoms that should be deleted.
            for residue in residue_map.keys():
                replaceWith = residue_map[residue]
                residue.name = replaceWith
                template = self._getTemplate(replaceWith)
                standardAtoms = set(atom.name for atom in template.topology.atoms())
                for atom in residue.atoms():
                    if (
                        atom.element in (None, hydrogen)
                        or atom.name not in standardAtoms
                    ):
                        deleteAtoms.append(atom)

            # Delete atoms queued to be deleted.
            modeller = app.Modeller(self.topology, self.positions)
            modeller.delete(deleteAtoms)
            self.topology = modeller.topology
            self.positions = modeller.positions

    def findMissingResidues(self, uncomplete_seq):
        chains = [c for c in self.topology.chains() if len(list(c.residues())) > 0]
        chainWithGaps = {}

        # Find the sequence of each chain, with gaps for missing residues.
        visited_chainids = set()
        for chain in chains:
            for sequence in uncomplete_seq:
                if chain.id == sequence.chainId:
                    if chain.id not in visited_chainids:
                        chainWithGaps[chain] = trim_none(sequence.residues)
                        visited_chainids.add(chain.id)

        # Try to find the chain that matches each sequence.
        chainSequence = {}
        chainOffset = {}
        for sequence in self.sequences:
            for chain in chains:
                if chain.id != sequence.chainId:
                    continue
                if chain in chainSequence:
                    continue
                for offset in range(
                    len(sequence.residues) - len(chainWithGaps[chain]) + 1
                ):
                    if all(
                        a == b or b == None
                        for a, b in zip(
                            sequence.residues[offset:], chainWithGaps[chain]
                        )
                    ):
                        chainSequence[chain] = sequence
                        chainOffset[chain] = offset
                        break
                if chain in chainSequence:
                    break

        # Now build the list of residues to add.

        self.missingResidues = {}
        for chain in self.topology.chains():
            if chain in chainSequence:
                offset = chainOffset[chain]
                sequence = chainSequence[chain].residues
                gappedSequence = chainWithGaps[chain]
                index = 0
                for i in range(len(sequence)):
                    if (
                        i < offset
                        or i >= len(gappedSequence) + offset
                        or gappedSequence[i - offset] is None
                    ):
                        key = (chain.index, index)
                        if key not in self.missingResidues:
                            self.missingResidues[key] = []
                        residueName = sequence[i]
                        if residueName in substitutions:
                            residueName = substitutions[sequence[i]]
                        self.missingResidues[key].append(residueName)
                    else:
                        index += 1


def assign_chains(fixer: PDBFixerWrapper):
    protein_chains_idxs = {}
    hetatom_chains_idxs = {}
    other_chains_idxs = []
    for chain_id, chain in enumerate(fixer.topology.chains()):
        tmp_res_names = [residue.name for residue in chain._residues]
        tmp_res_names_set = list(set(tmp_res_names))
        if len(tmp_res_names) > 20 and all(
            [x in PROTEIN_RESIDUES + list(METALS) for x in tmp_res_names_set]
        ):
            protein_chains_idxs[chain.index] = chain.id

        elif all([x in WATERS + METALS for x in tmp_res_names_set]):
            hetatom_chains_idxs[chain.index] = chain.id

        else:
            other_chains_idxs.append(chain.index)

    return (
        protein_chains_idxs,
        hetatom_chains_idxs,
        list(other_chains_idxs),
    )


def pdb_to_fixer(pdb_fn, cif_fn=None):
    with open(pdb_fn, "r") as f:
        fixer = PDBFixerWrapper(pdbfile=f)

    # remove hetero atoms by residue name
    fixer.findNonstandardResidues()
    # 记录一些非标氨基酸
    std_residues_mapping = {
        residue.name: std_resname for residue, std_resname in fixer.nonstandardResidues
    }
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=True, keepMetal=True)
    (
        protein_chain_idxs,
        hetatom_chains_idxs,
        other_chains_idxs,
    ) = assign_chains(fixer)
    fixer.removeChains(chainIndices=other_chains_idxs)

    # 将首尾变成ACE/NME
    chain_mutations = defaultdict(list)
    not_same_chain_idxs = []
    for chain in fixer.topology.chains():
        if chain.id in protein_chain_idxs.values():
            start_residue, end_residue = None, None
            for residue in chain._residues[:10]:
                if residue.name in PROTEIN_RESIDUES:
                    start_residue = residue
                    chain_mutations[chain.id].append(
                        f"{start_residue.name}:{start_residue.id}{start_residue.insertionCode}:ACE"
                    )
                    break

            for residue in chain._residues[::-1][:10]:
                if residue.name in PROTEIN_RESIDUES:
                    end_residue = residue
                    chain_mutations[chain.id].append(
                        f"{end_residue.name}:{end_residue.id}{end_residue.insertionCode}:NME"
                    )
                    break

            if (start_residue is not None) or (end_residue is not None):
                chain_len1 = len(
                    [res for res in chain._residues if res.insertionCode == " "]
                )
                chain_len2 = int(end_residue.id) - int(start_residue.id) + 1

                if chain_len1 != chain_len2:
                    not_same_chain_idxs.append(chain.id)

    uncomplete_seqs = []
    if len(not_same_chain_idxs) > 0:
        complete_sequences, uncomplete_sequences = get_sequences_info(cif_fn)
        complete_seqs = []
        for sequence, nan_sequence in zip(complete_sequences, uncomplete_sequences):
            if sequence.chainId in not_same_chain_idxs:
                if len(std_residues_mapping) > 0:
                    new_residues = [
                        std_residues_mapping[x] if x in std_residues_mapping else x
                        for x in sequence.residues
                    ]
                    sequence.residues = new_residues

                    new_residues1 = [
                        std_residues_mapping[x] if x in std_residues_mapping else x
                        for x in nan_sequence.residues
                    ]
                    nan_sequence.residues = new_residues1
                complete_seqs.append(sequence)
                uncomplete_seqs.append(nan_sequence)
        fixer.sequences = complete_seqs

    fixer.findMissingResidues(uncomplete_seqs)
    print(f"{fixer.missingResidues=}")

    print(f"Remove end missing residues, and it seems to be not necessary...")
    to_delete_keys = set()
    for key, values in fixer.missingResidues.items():
        if key[1] == 0:
            to_delete_keys.add(key)

        chain = fixer.topology._chains[key[0]]
        if key[1] == len(chain._residues):
            to_delete_keys.add(key)

        new_values = []
        for name in values:
            if name in substitutions:
                new_values.append(substitutions[name])
            elif name in PROTEIN_RESIDUES:
                new_values.append(name)
            else:
                new_values.append("ALA")
        fixer.missingResidues[key] = new_values

    for key in to_delete_keys:
        del fixer.missingResidues[key]
    print(f"{fixer.missingResidues=}")

    # add capped residues
    for chain_name, mutates in chain_mutations.items():
        fixer.applyMutations(mutates, chain_name)

    fixer.findMissingAtoms()
    print(f"{fixer.missingAtoms=}")
    fixer.addMissingAtoms()
    # fixer.addMissingHydrogens(7.0)

    # with open("test.pdb", "w") as f:
    #     app.PDBFile.writeFile(fixer.topology, fixer.positions, file=f, keepIds=True)
    return fixer


def fixer_into_protein_mol(fixer: PDBFixer, residue_tables):
    connect_table, charge_table, hydrogen_table = residue_tables
    atom_mapper = dict()

    new_positions = []
    histide_residues = {}
    hip_set = set()

    new_mol = Chem.RWMol()
    for chain in fixer.topology.chains():
        residue_records = {}
        for atom in chain.atoms():
            key = [
                atom.residue.name,
                atom.residue.id,
                atom.residue.chain.id,
                atom.name,
                atom.index,
                atom.residue.insertionCode,
            ]
            # if int(key[1])>=60 and int(key[1])<=65:
            #     print(key)

            if key[3] in ["OXT"]:
                continue
            if tuple(key[:3]) not in residue_records and key[0] not in ["ACE"]:
                residue_records[tuple(key[:3])] = len(residue_records) + 1

            # 去除第一个残基的多余的H2，H3
            if residue_records.get(tuple(key[:3]), False) == 1 and key[3] in [
                "H2",
                "H3",
            ]:
                continue

            rdatom = Chem.Atom(atom.element.atomic_number)
            mi = Chem.AtomPDBResidueInfo()
            # if key[0] == "FE":
            #     key[0] = "FE2"
            mi.SetResidueName(f"{key[0]:>3s}")
            mi.SetResidueNumber(int(key[1]))
            mi.SetChainId(key[2])
            if atom.residue.insertionCode == "":
                mi.SetInsertionCode(" ")
            else:
                mi.SetInsertionCode(atom.residue.insertionCode)

            if key[0] in ["NME"] and key[3] == "C":
                mi.SetName(format_4letter("CH3"))
            # elif key[3] == "FE":
            #     mi.SetName(format_4letter("FE2"))
            else:
                mi.SetName(format_4letter(key[3]))
            rdatom.SetMonomerInfo(mi)

            try:
                chg = charge_table[key[0]][key[3]]
            except:
                chg = 0

            rdatom.SetFormalCharge(chg)
            atom_idx = new_mol.AddAtom(rdatom)
            atom_mapper[key[4]] = atom_idx

            tmp_coord = fixer.positions[key[4]]
            new_positions.append([tmp_coord.x, tmp_coord.y, tmp_coord.z])
            if key[0] == "HIS":
                if tuple(key[:3]) not in histide_residues:
                    histide_residues[tuple(key[:3])] = []
                histide_residues[tuple(key[:3])].append(atom_idx)
                if key[3] == "HE2":
                    hip_set.add(tuple(key[:3]))

    # 更新HIP残基的状态
    if len(hip_set) > 0:
        for hip_key in list(hip_set):
            for atom_idx in histide_residues[hip_key]:
                atom = new_mol.GetAtomWithIdx(atom_idx)
                mi = atom.GetPDBResidueInfo()
                mi.SetResidueName("HIP")
                chg = charge_table["HIP"][mi.GetName().strip()]
                atom.SetFormalCharge(chg)
                atom.SetMonomerInfo(mi)

    for bond in fixer.topology.bonds():
        atom1 = bond.atom1
        atom2 = bond.atom2
        if atom1.index in atom_mapper and atom2.index in atom_mapper:
            at1 = atom_mapper[atom1.index]
            at2 = atom_mapper[atom2.index]
            key1 = (atom1.residue.name, atom1.residue.id, atom1.residue.chain.id)
            key2 = (atom2.residue.name, atom2.residue.id, atom2.residue.chain.id)

            if key1 == key2:
                try:
                    res_name = "HIP" if key1 in hip_set else key1[0]
                    bond_value = connect_table[res_name][(atom1.name, atom2.name)]
                except:
                    bond_value = 1

            else:
                bond_value = 1

            bond = new_mol.GetBondBetweenAtoms(at1, at2)
            if bond is None:
                if bond_value == 1:
                    new_mol.AddBond(at1, at2, Chem.BondType.SINGLE)
                elif bond_value == 2:
                    new_mol.AddBond(at1, at2, Chem.BondType.DOUBLE)
                else:
                    print(f"warning {bond_value=}")
                    new_mol.AddBond(at1, at2, Chem.BondType.SINGLE)

    new_mol = new_mol.GetMol()
    problems = Chem.DetectChemistryProblems(new_mol)
    for problem in problems:
        cur_idx = problem.GetAtomIdx()
        atom = new_mol.GetAtomWithIdx(cur_idx)
        pdb_info = atom.GetPDBResidueInfo()
        print(
            pdb_info.GetResidueName(),
            pdb_info.GetResidueNumber(),
            pdb_info.GetChainId(),
            pdb_info.GetName(),
        )
        print("connected atoms:")
        for nbr_atom in atom.GetNeighbors():
            nbr_atom: Chem.Atom
            pdb_info1 = nbr_atom.GetPDBResidueInfo()
            print(
                pdb_info1.GetResidueName(),
                pdb_info1.GetResidueNumber(),
                pdb_info1.GetChainId(),
                pdb_info1.GetName(),
            )

    new_positions = np.array(new_positions) * 10.0
    new_mol = add_conformer_to_mol(new_mol, new_positions)
    new_mol = Chem.RemoveAllHs(new_mol, sanitize=False)
    Chem.SanitizeMol(new_mol)
    new_mol = add_h_to_receptor_mol(new_mol, hydrogen_table)

    return new_mol


def protein_mol_to_file(protein_mol: Chem.Mol, out_fn: str | Path = None):
    # Chem.MolToPDBFile(protein_mol, "debug.pdb")
    receptor_omm = app.PDBFile(StringIO(Chem.MolToPDBBlock(protein_mol)))
    if out_fn is not None:
        with open(out_fn, "w") as f:
            app.PDBFile.writeFile(
                receptor_omm.topology, receptor_omm.positions, f, keepIds=True
            )


def ligand_rdmol_to_file(ligand_mol: Chem.Mol, out_fn: str) -> None:
    with Chem.SDWriter(out_fn) as writer:
        ligand_mol = Chem.RemoveHs(ligand_mol)
        writer.write(ligand_mol)
    return
