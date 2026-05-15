import os
import csv
import argparse
import numpy as np
from Bio.PDB import PDBParser, is_aa
from Bio.PDB.SASA import ShrakeRupley


def normalize_resname(resname):
    return str(resname).strip().upper()


def load_npz_resnames(npz_path):
    data = np.load(npz_path, allow_pickle=True)
    if "res_name3" not in data:
        raise KeyError("npz has no 'res_name3'")
    res_name3 = [normalize_resname(x) for x in data["res_name3"]]
    return res_name3


def extract_chain_id_from_filename(pdb_path):

    base = os.path.basename(pdb_path)
    name, _ = os.path.splitext(base)

    if "_" in name:
        parts = name.split("_")
        if len(parts[-1]) == 1:
            return parts[-1]

    return None


def compute_residue_sasa_from_pdb(pdb_path, chain_id=None):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_path)

    sr = ShrakeRupley()
    sr.compute(structure, level="R")

    model = structure[0]

    if chain_id is None:
        chains = list(model.get_chains())
        if len(chains) != 1:
            raise ValueError(
                f"pdb has multiple chains ({[c.id for c in chains]}), "
                f"have to cite chain_id"
            )
        chain = chains[0]
    else:
        if chain_id not in model:
            raise ValueError(f"Цепь '{chain_id}' not foind in PDB")
        chain = model[chain_id]

    residues_info = []

    for residue in chain:
        # residue.id = (hetflag, resseq, icode)
        hetflag, resseq, icode = residue.id

        # Только стандартные аминокислоты
        if not is_aa(residue, standard=True):
            continue

        resname = normalize_resname(residue.resname)
        sasa = float(getattr(residue, "sasa", 0.0))

        residues_info.append({
            "resname": resname,
            "resseq": int(resseq),
            "icode": str(icode).strip(),
            "sasa": sasa,
        })

    return residues_info


def compare_npz_and_pdb(npz_resnames, pdb_residues):
    pdb_resnames = [r["resname"] for r in pdb_residues]

    print("\nLength comparison")
    print(f"len(res_name3 from npz) = {len(npz_resnames)}")
    print(f"len(aa residues from pdb) = {len(pdb_resnames)}")

    min_len = min(len(npz_resnames), len(pdb_resnames))
    mismatches = []

    for i in range(min_len):
        if npz_resnames[i] != pdb_resnames[i]:
            mismatches.append((i, npz_resnames[i], pdb_resnames[i]))

    print("\n=== ПЕРВЫЕ 15 ОСТАТКОВ ===")
    for i in range(min(15, min_len)):
        print(f"{i:3d}: npz={npz_resnames[i]:>3s} | pdb={pdb_resnames[i]:>3s}")

    print("\nOutput")
    if len(npz_resnames) != len(pdb_resnames):
        print("Lengths do not match.")
    else:
        print("Lengths match.")

    if mismatches:
        print(f"Inconsistencies in residue names: {len(mismatches)}")
        print("First ten inconsistencies:")
        for i, a, b in mismatches[:10]:
            print(f"  idx={i}: npz={a}, pdb={b}")
    else:
        print("Names match by order.")

    return mismatches


def save_sasa_csv(pdb_residues, out_csv):
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["index0", "resname", "resseq", "icode", "sasa"])

        for i, r in enumerate(pdb_residues):
            writer.writerow([i, r["resname"], r["resseq"], r["icode"], r["sasa"]])


def save_sasa_npy(pdb_residues, out_npy):
    sasa_array = np.array([r["sasa"] for r in pdb_residues], dtype=np.float64)
    np.save(out_npy, sasa_array)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--npz", required=True, help="Path to .npz")
    parser.add_argument("--pdb", required=True, help="Path .pdb")
    parser.add_argument("--chain", default=None, help="Chain ID")
    args = parser.parse_args()

    npz_resnames = load_npz_resnames(args.npz)

    chain_id = args.chain
    if chain_id is None:
        chain_id = extract_chain_id_from_filename(args.pdb)

    print("NPZ:", args.npz)
    print("PDB:", args.pdb)
    print("Chain:", chain_id)

    pdb_residues = compute_residue_sasa_from_pdb(args.pdb, chain_id=chain_id)
    mismatches = compare_npz_and_pdb(npz_resnames, pdb_residues)

    base = os.path.splitext(os.path.basename(args.pdb))[0]
    out_csv = f"{base}_residue_sasa.csv"
    out_npy = f"{base}_residue_sasa.npy"

    save_sasa_csv(pdb_residues, out_csv)
    save_sasa_npy(pdb_residues, out_npy)

    print("\n=== ФАЙЛЫ СОХРАНЕНЫ ===")
    print("CSV:", out_csv)
    print("NPY:", out_npy)

    if len(npz_resnames) == len(pdb_residues) and not mismatches:
        print("\nSuccess.")
    else:
        print("\nFailure, have to check inconsistencies in residue order/length.")


if __name__ == "__main__":
    main()