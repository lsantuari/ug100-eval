#!/usr/bin/env python3
"""
Create a side-by-side publication-ready scatterplot of per-gene TP/FP differences
for SNVs and Indels between two hap.py results.

Output:
  - High-quality PNG, PDF, and EPS figures
  - Shared main title and consistent axes
"""

import pandas as pd
import matplotlib.pyplot as plt
from adjustText import adjust_text
import argparse

def plot_one(ax, df, label1, label2, variant_type, top_n=10):
    """Plot a single panel (SNV or Indel)"""
    df["ΔTP"] = df[f"TP_{label1}"] - df[f"TP_{label2}"]
    df["ΔFP"] = df[f"FP_{label1}"] - df[f"FP_{label2}"]
    df["abs_diff"] = (df["ΔTP"].abs() + df["ΔFP"].abs())
    outliers = df.nlargest(top_n, "abs_diff")

    ax.scatter(df["ΔTP"], df["ΔFP"], s=25, color="steelblue", alpha=0.6, edgecolor="none")
    ax.axhline(0, color="gray", linestyle="--", linewidth=0.8)
    ax.axvline(0, color="gray", linestyle="--", linewidth=0.8)

    # Label outliers
    texts = [ax.text(row["ΔTP"], row["ΔFP"], row["Gene"], fontsize=8, color="black")
             for _, row in outliers.iterrows()]
    adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", lw=0.4, color="gray"))

    ax.set_title(variant_type, fontsize=12, pad=8)
    ax.set_xlabel(r"$\Delta$TP (" + f"{label1} − {label2})", fontsize=11)
    ax.set_ylabel(r"$\Delta$FP (" + f"{label1} − {label2})", fontsize=11)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    pad_x = max(1, df["ΔTP"].abs().max() * 0.1)
    pad_y = max(1, df["ΔFP"].abs().max() * 0.1)
    ax.set_xlim(df["ΔTP"].min() - pad_x, df["ΔTP"].max() + pad_x)
    ax.set_ylim(df["ΔFP"].min() - pad_y, df["ΔFP"].max() + pad_y)


def main():
    ap = argparse.ArgumentParser(description="Combined SNV/Indel ΔTP vs ΔFP scatterplot")
    ap.add_argument("--snv-csv", required=True, help="CSV with SNV per-gene results")
    ap.add_argument("--indel-csv", required=True, help="CSV with Indel per-gene results")
    ap.add_argument("--label1", required=True)
    ap.add_argument("--label2", required=True)
    ap.add_argument("--out-prefix", default="gene_diff_combined")
    ap.add_argument("--label-top-n", type=int, default=10)
    ap.add_argument("--main-title", default="Per-gene ΔTP and ΔFP differences")
    args = ap.parse_args()

    # Load data
    df_snv = pd.read_csv(args.snv_csv)
    df_indel = pd.read_csv(args.indel_csv)

    # --- Global style ---
    plt.rcParams.update({
        "font.size": 11,
        "axes.linewidth": 1.0,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "font.family": "serif",
        "font.serif": ["Times New Roman", "DejaVu Serif", "Times"],
    })

    # --- Create figure with 2 side-by-side panels ---
    fig, axes = plt.subplots(1, 2, figsize=(9, 4.5), sharex=False, sharey=False)

    plot_one(axes[0], df_snv, args.label1, args.label2, "SNVs", args.label_top_n)
    plot_one(axes[1], df_indel, args.label1, args.label2, "Indels", args.label_top_n)

    # Main title
    fig.suptitle(args.main_title, fontsize=13, weight="bold", y=1.03)

    plt.tight_layout()

    # --- Save outputs ---
    for ext in ["png", "pdf", "eps"]:
        outfile = f"{args.out_prefix}_TP_FP_diff.{ext}"
        plt.savefig(outfile, dpi=300, bbox_inches="tight", format=ext)
        print(f"[✔] Saved: {outfile}")

    plt.close()


if __name__ == "__main__":
    main()
