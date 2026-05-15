import argparse
import os
import numpy as np

import pyrosetta
from pyrosetta import rosetta

def init_pyrosetta():
    pyrosetta.init("-mute all")

def weighted_sum(emap, weights, terms):
    s = 0.0
    for t in terms:
        w = float(weights[t])
        if w != 0.0:
            s += w * float(emap[t])
    return s

def build_tensor_for_pose(pose, protein_only=True):
    scorefxn = rosetta.core.scoring.get_score_function()
    scorefxn(pose)

    weights = scorefxn.weights()
    ST = rosetta.core.scoring.ScoreType

    vdw_terms = [ST.fa_atr, ST.fa_rep]
    elec_terms = [ST.fa_elec]
    hb_terms = [ST.hbond_sr_bb, ST.hbond_lr_bb, ST.hbond_bb_sc, ST.hbond_sc]

    active_terms = [t for t in ST.__members__.values() if float(weights[t]) != 0.0]
    known = set(vdw_terms + [ST.fa_sol] + elec_terms + hb_terms)
    mech_terms = [t for t in active_terms if t not in known]

    if protein_only:
        idx_pose_to_tensor = {}
        tensor_to_pose = []
        k = 0
        for i in range(1, pose.total_residue() + 1):
            if pose.residue(i).is_protein():
                idx_pose_to_tensor[i] = k
                tensor_to_pose.append(i)
                k += 1
        N = k
    else:
        N = pose.total_residue()
        idx_pose_to_tensor = {i: i - 1 for i in range(1, N + 1)}
        tensor_to_pose = list(range(1, N + 1))

    T = np.zeros((N, N, 5), dtype=np.float32)

    egraph = pose.energies().energy_graph()
    if hasattr(egraph, "find_edge"):
        get_edge = egraph.find_edge
    elif hasattr(egraph, "get_edge"):
        get_edge = egraph.get_edge
    else:
        raise RuntimeError("energy_graph has neither find_edge nor get_edge")

    for i in range(1, pose.total_residue() + 1):
        for j in range(i + 1, pose.total_residue() + 1):
            if protein_only and ((i not in idx_pose_to_tensor) or (j not in idx_pose_to_tensor)):
                continue

            edge = get_edge(i, j)
            if edge is None:
                continue

            a = idx_pose_to_tensor[i] if protein_only else i - 1
            b = idx_pose_to_tensor[j] if protein_only else j - 1

            emap = edge.fill_energy_map()

            vdw = weighted_sum(emap, weights, vdw_terms)
            ele = weighted_sum(emap, weights, elec_terms)
            hb = weighted_sum(emap, weights, hb_terms)
            mech = weighted_sum(emap, weights, mech_terms)

            T[a, b, 0] = vdw
            T[a, b, 2] = ele
            T[a, b, 3] = mech
            T[a, b, 4] = hb

            T[b, a, 0] = vdw
            T[b, a, 2] = ele
            T[b, a, 3] = mech
            T[b, a, 4] = hb

    energies = pose.energies()
    for pose_i in tensor_to_pose:
        a = idx_pose_to_tensor[int(pose_i)] if protein_only else int(pose_i) - 1
        emap_i = energies.residue_total_energies(int(pose_i))
        T[a, a, 1] = float(emap_i[ST.fa_sol])

    return T, np.array(tensor_to_pose, dtype=np.int32)


def extract_parsed_structure_for_tensor(pose, tensor_to_pose):
    pdbinfo = pose.pdb_info()

    res_chain = []
    res_num = []
    res_icode = []
    res_name3 = []
    ca_xyz = []

    for pose_i in tensor_to_pose:
        r = pose.residue(int(pose_i))
        res_name3.append(r.name3())

        if pdbinfo is not None:
            res_chain.append(pdbinfo.chain(int(pose_i)))
            res_num.append(pdbinfo.number(int(pose_i)))
            res_icode.append(pdbinfo.icode(int(pose_i)))
        else:
            res_chain.append("?")
            res_num.append(int(pose_i))
            res_icode.append("")

        if r.has("CA"):
            v = r.xyz("CA")
            ca_xyz.append([float(v.x), float(v.y), float(v.z)])
        else:
            ca_xyz.append([np.nan, np.nan, np.nan])

    return {
        "tensor_to_pose": np.array(tensor_to_pose, dtype=np.int32),
        "res_chain": np.array(res_chain, dtype=object),
        "res_num": np.array(res_num, dtype=np.int32),
        "res_icode": np.array(res_icode, dtype=object),
        "res_name3": np.array(res_name3, dtype=object),
        "ca_xyz": np.array(ca_xyz, dtype=np.float32),
    }

def detect_pdb_gaps(pose, protein_only=True):
    pdbinfo = pose.pdb_info()
    if pdbinfo is None:
        return np.empty((0, 7), dtype=object)

    records = []
    prev = None

    for i in range(1, pose.total_residue() + 1):
        r = pose.residue(i)
        if protein_only and (not r.is_protein()):
            continue

        ch = pdbinfo.chain(i)
        num = pdbinfo.number(i)
        ic = pdbinfo.icode(i)

        if prev is not None:
            pch, pnum, pic, pi = prev
            if ch == pch and (num - pnum) > 1:
                records.append((ch, pnum, pic, num, ic, pi, i))

        prev = (ch, num, ic, i)

    return np.array(records, dtype=object)


def contact_matrix_from_ca(ca_xyz, cutoff):
    ca_xyz = np.asarray(ca_xyz, dtype=np.float32)
    diff = ca_xyz[:, None, :] - ca_xyz[None, :, :]
    dist = np.sqrt(np.sum(diff ** 2, axis=-1))

    C = (dist <= cutoff).astype(np.int8)
    np.fill_diagonal(C, 0)
    return C


def build_contact_matrices(ca_xyz, cutoffs=(4.0, 6.0, 8.0)):
    contacts = {}
    for cutoff in cutoffs:
        contacts[f"C{int(cutoff)}"] = contact_matrix_from_ca(ca_xyz, cutoff)
    return contacts


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb", required=True, help="Path to input PDB file")
    parser.add_argument("--out", required=True, help="Output .npz path")
    parser.add_argument("--protein_only", action="store_true", help="Use only protein residues")
    args = parser.parse_args()

    pose = pyrosetta.pose_from_pdb(args.pdb)
    T, tensor_to_pose = build_tensor_for_pose(pose, protein_only=args.protein_only)
    parsed = extract_parsed_structure_for_tensor(pose, tensor_to_pose)
    gaps = detect_pdb_gaps(pose, protein_only=args.protein_only)
    contacts = build_contact_matrices(parsed["ca_xyz"], cutoffs=(4.0, 6.0, 8.0))

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)

    np.savez_compressed(
        args.out,
        T=T,
        C4=contacts["C4"],
        C6=contacts["C6"],
        C8=contacts["C8"],
        gaps=gaps,
        pdb_path=np.array([args.pdb], dtype=object),
        protein_only=np.array([args.protein_only], dtype=bool),
        **parsed,
    )

    print(f"Saved unified file to: {args.out}")
    print(f"T shape: {T.shape}")
    print(f"C4 shape: {contacts['C4'].shape}")
    print(f"C6 shape: {contacts['C6'].shape}")
    print(f"C8 shape: {contacts['C8'].shape}")
    print(f"Residues: {len(tensor_to_pose)}")


if __name__ == "__main__":
    init_pyrosetta()
    main()