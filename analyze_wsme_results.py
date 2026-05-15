import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def extract_protein_id_from_npz_path(npz_path: str) -> str:
    base = os.path.basename(npz_path)
    name, _ = os.path.splitext(base)
    if name.endswith("_full"):
        name = name[:-5]
    return name


def load_contact_counts(npz_dir: str, protein_ids: list[str], cutoff: str = "C8", min_seq_sep: int = 3) -> dict[str, int]:
    counts = {}

    for protein_id in protein_ids:
        npz_path = os.path.join(npz_dir, f"{protein_id}_full.npz")
        if not os.path.exists(npz_path):
            counts[protein_id] = np.nan
            continue

        try:
            data = np.load(npz_path, allow_pickle=True)
            if cutoff not in data:
                counts[protein_id] = np.nan
                continue

            C = data[cutoff].astype(np.float64)
            N = C.shape[0]

            # remove local contacts exactly as in your fitting scripts
            for i in range(N):
                for j in range(N):
                    if abs(i - j) < min_seq_sep:
                        C[i, j] = 0.0

            # count unique residue-residue contacts
            contact_count = int(np.sum(np.triu(C, k=1) > 0))
            counts[protein_id] = contact_count

        except Exception:
            counts[protein_id] = np.nan

    return counts


def safe_corr(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 2:
        return np.nan
    if np.allclose(np.std(x), 0) or np.allclose(np.std(y), 0):
        return np.nan
    return float(np.corrcoef(x, y)[0, 1])


def fit_line(x: np.ndarray, y: np.ndarray):
    if len(x) < 2:
        return None
    if np.allclose(np.std(x), 0):
        return None
    coeffs = np.polyfit(x, y, 1)
    return coeffs


def save_summary_txt(df_valid: pd.DataFrame, out_txt: str, corr_dg: float, r2_dg: float, corr_contacts_exp: float, corr_contacts_pred: float):
    rmse = float(np.sqrt(np.mean(df_valid["sq_error"])))
    mae = float(np.mean(np.abs(df_valid["error"])))

    with open(out_txt, "w", encoding="utf-8") as f:
        f.write("WSME RESULT SUMMARY\n")
        f.write("===================\n\n")
        f.write(f"n_valid = {len(df_valid)}\n")
        f.write(f"RMSE = {rmse:.6f}\n")
        f.write(f"MAE = {mae:.6f}\n")
        f.write(f"Pearson(predicted vs experimental dG) = {corr_dg:.6f}\n")
        f.write(f"R^2(predicted vs experimental dG) = {r2_dg:.6f}\n")
        f.write(f"Pearson(contacts vs experimental dG) = {corr_contacts_exp:.6f}\n")
        f.write(f"Pearson(contacts vs predicted dG) = {corr_contacts_pred:.6f}\n\n")

        f.write("BEST/WORST ABSOLUTE ERRORS\n")
        f.write("=========================\n")
        df_sorted = df_valid.copy()
        df_sorted["abs_error"] = df_sorted["error"].abs()
        for _, row in df_sorted.sort_values("abs_error", ascending=False).iterrows():
            f.write(
                f"{row['protein_id']:15s} "
                f"target={row['target_dG']:.3f} "
                f"pred={row['predicted_dG']:.3f} "
                f"err={row['error']:.3f} "
                f"Eact={row['E_activation']:.3f}\n"
            )


def make_scatter_dg(df_valid: pd.DataFrame, out_png: str):
    x = df_valid["target_dG"].to_numpy(dtype=float)
    y = df_valid["predicted_dG"].to_numpy(dtype=float)

    corr = safe_corr(x, y)
    r2 = corr ** 2 if np.isfinite(corr) else np.nan

    plt.figure(figsize=(7, 6))
    plt.scatter(x, y)

    xy_min = min(np.min(x), np.min(y))
    xy_max = max(np.max(x), np.max(y))
    plt.plot([xy_min, xy_max], [xy_min, xy_max], linestyle="--")

    line = fit_line(x, y)
    if line is not None:
        a, b = line
        xx = np.linspace(xy_min, xy_max, 200)
        yy = a * xx + b
        plt.plot(xx, yy)

    plt.xlabel("Experimental ΔG (kcal/mol)")
    plt.ylabel("Predicted ΔG (kcal/mol)")
    plt.title(f"Predicted vs Experimental ΔG\nPearson = {corr:.3f}, R² = {r2:.3f}")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()

    return corr, r2


def make_contacts_plot(df_valid: pd.DataFrame, out_png: str):
    df_contacts = df_valid.dropna(subset=["n_contacts"]).copy()
    if len(df_contacts) == 0:
        return np.nan, np.nan

    x = df_contacts["n_contacts"].to_numpy(dtype=float)
    y_exp = df_contacts["target_dG"].to_numpy(dtype=float)
    y_pred = df_contacts["predicted_dG"].to_numpy(dtype=float)

    corr_exp = safe_corr(x, y_exp)
    corr_pred = safe_corr(x, y_pred)

    plt.figure(figsize=(7, 6))
    plt.scatter(x, y_exp, label="Experimental ΔG")
    plt.scatter(x, y_pred, label="Predicted ΔG")

    line_exp = fit_line(x, y_exp)
    if line_exp is not None:
        a, b = line_exp
        xx = np.linspace(np.min(x), np.max(x), 200)
        plt.plot(xx, a * xx + b)

    line_pred = fit_line(x, y_pred)
    if line_pred is not None:
        a, b = line_pred
        xx = np.linspace(np.min(x), np.max(x), 200)
        plt.plot(xx, a * xx + b)

    plt.xlabel("Number of contacts")
    plt.ylabel("ΔG (kcal/mol)")
    plt.title(
        "ΔG vs Number of Contacts\n"
        f"corr(exp) = {corr_exp:.3f}, corr(pred) = {corr_pred:.3f}"
    )
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()

    return corr_exp, corr_pred


def make_barrier_plot(df_valid: pd.DataFrame, out_png: str):
    df_plot = df_valid.sort_values("E_activation", ascending=False).copy()

    plt.figure(figsize=(10, 6))
    plt.bar(df_plot["protein_id"], df_plot["E_activation"])
    plt.xticks(rotation=90)
    plt.xlabel("Protein")
    plt.ylabel("Activation barrier Eact (kcal/mol)")
    plt.title("Activation Barrier Heights")
    plt.grid(True, axis="y")
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()


def make_error_plot(df_valid: pd.DataFrame, out_png: str):
    df_plot = df_valid.copy()
    df_plot["abs_error"] = df_plot["error"].abs()
    df_plot = df_plot.sort_values("abs_error", ascending=False)

    plt.figure(figsize=(10, 6))
    plt.bar(df_plot["protein_id"], df_plot["error"])
    plt.xticks(rotation=90)
    plt.xlabel("Protein")
    plt.ylabel("Prediction error (pred - exp), kcal/mol")
    plt.title("Per-protein ΔG Errors")
    plt.grid(True, axis="y")
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", required=True, help="CSV with fitting results")
    parser.add_argument("--npz_dir", default=None, help="Folder with *_full.npz files for contact counting")
    parser.add_argument("--cutoff", default="C8", choices=["C4", "C6", "C8"])
    parser.add_argument("--min_seq_sep", type=int, default=3)
    parser.add_argument("--out_dir", default="wsme_analysis")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    df = pd.read_csv(args.csv)

    required_cols = ["protein_id", "target_dG", "predicted_dG", "error", "sq_error", "E_activation"]
    for col in required_cols:
        if col not in df.columns:
            raise ValueError(f"CSV is missing required column: {col}")

    df_valid = df.dropna(subset=["target_dG", "predicted_dG", "error", "sq_error", "E_activation"]).copy()

    if args.npz_dir is not None:
        contact_counts = load_contact_counts(
            npz_dir=args.npz_dir,
            protein_ids=df["protein_id"].astype(str).tolist(),
            cutoff=args.cutoff,
            min_seq_sep=args.min_seq_sep,
        )
        df["n_contacts"] = df["protein_id"].map(contact_counts)
        df_valid["n_contacts"] = df_valid["protein_id"].map(contact_counts)
    else:
        df["n_contacts"] = np.nan
        df_valid["n_contacts"] = np.nan

    corr_dg, r2_dg = make_scatter_dg(
        df_valid,
        os.path.join(args.out_dir, "predicted_vs_experimental_dG.png"),
    )

    corr_contacts_exp, corr_contacts_pred = make_contacts_plot(
        df_valid,
        os.path.join(args.out_dir, "dG_vs_contacts.png"),
    )

    make_barrier_plot(
        df_valid,
        os.path.join(args.out_dir, "activation_barriers.png"),
    )

    make_error_plot(
        df_valid,
        os.path.join(args.out_dir, "per_protein_errors.png"),
    )

    save_summary_txt(
        df_valid,
        os.path.join(args.out_dir, "summary.txt"),
        corr_dg=corr_dg,
        r2_dg=r2_dg,
        corr_contacts_exp=corr_contacts_exp,
        corr_contacts_pred=corr_contacts_pred,
    )

    df.to_csv(os.path.join(args.out_dir, "results_with_contacts.csv"), index=False)

    print("Done.")
    print("Saved files to:", args.out_dir)
    print(" - predicted_vs_experimental_dG.png")
    print(" - dG_vs_contacts.png")
    print(" - activation_barriers.png")
    print(" - per_protein_errors.png")
    print(" - summary.txt")
    print(" - results_with_contacts.csv")


if __name__ == "__main__":
    main()