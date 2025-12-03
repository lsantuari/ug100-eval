#!/usr/bin/env python3

import pysam
import pybedtools
import csv
import tempfile
from collections import defaultdict
import argparse
import math

def write_csv_tables(merged, vt, l1, l2, out_prefix):
    """
    Saves both filtered and full CSV tables.
    """
    tp_key, fp_key = f"{vt}_TP", f"{vt}_FP"
    
    # Full CSV (all genes)
    full_csv_file = f"{out_prefix}_{vt}_full.csv"
    with open(full_csv_file, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Gene", f"TP_{l1}", f"FP_{l1}", f"TP_{l2}", f"FP_{l2}"])
        for gene, vals in sorted(merged.items()):
            v1, v2 = vals[l1], vals[l2]
            writer.writerow([gene, v1[tp_key], v1[fp_key], v2[tp_key], v2[fp_key]])
    print(f"[✔] Saved full CSV: {full_csv_file}")

    # Filtered CSV (genes with differing TP/FP)
    filtered_csv_file = f"{out_prefix}_{vt}_filtered.csv"
    with open(filtered_csv_file, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Gene", f"TP_{l1}", f"FP_{l1}", f"TP_{l2}", f"FP_{l2}"])
        for gene, vals in sorted(merged.items()):
            v1, v2 = vals[l1], vals[l2]
            if v1[tp_key] != v2[tp_key] or v1[fp_key] != v2[fp_key]:
                writer.writerow([gene, v1[tp_key], v1[fp_key], v2[tp_key], v2[fp_key]])
    print(f"[✔] Saved filtered CSV: {filtered_csv_file}")


def vcf_to_bed(vcf_path, label):
    vcf = pysam.VariantFile(vcf_path)
    second_sample = list(vcf.header.samples)[1]
    bed_file = tempfile.NamedTemporaryFile(prefix=f"{label}_", suffix=".bed", delete=False)

    with open(bed_file.name, "w") as out:
        for rec in vcf:
            bd = rec.samples[second_sample].get("BD")
            if not bd:
                continue
            bd = bd.upper()
            if bd not in ("TP", "FP"):
                continue

            if len(rec.ref) == 1 and all(len(a) == 1 for a in rec.alts):
                var_type = "SNV"
            else:
                max_len = max(abs(len(a) - len(rec.ref)) for a in rec.alts)
                if max_len > 50:
                    continue
                var_type = "Indel"

            start = rec.pos - 1
            end = rec.pos + len(rec.ref) - 1
            out.write(f"{rec.chrom}\t{start}\t{end}\t{var_type}\t{bd}\n")

    return pybedtools.BedTool(bed_file.name)


def count_by_gene(vcf_path, gene_bed, label):
    vbed = vcf_to_bed(vcf_path, label)
    intersected = vbed.intersect(gene_bed, wa=True, wb=True)
    counts = defaultdict(lambda: {"SNV_TP": 0, "SNV_FP": 0, "Indel_TP": 0, "Indel_FP": 0})
    for f in intersected:
        var_type, bd, gene = f[3], f[4], f[8]
        key = f"{var_type}_{bd}"
        counts[gene][key] += 1
    return {label: counts}


def merge_counts(c1, c2, l1, l2):
    genes = set(c1[l1].keys()) | set(c2[l2].keys())
    merged = {}
    for g in genes:
        merged[g] = {
            l1: c1[l1].get(g, {"SNV_TP": 0, "SNV_FP": 0, "Indel_TP": 0, "Indel_FP": 0}),
            l2: c2[l2].get(g, {"SNV_TP": 0, "SNV_FP": 0, "Indel_TP": 0, "Indel_FP": 0}),
        }
    return merged


def filter_differences(merged, vt, l1, l2):
    tp_key, fp_key = f"{vt}_TP", f"{vt}_FP"
    filtered = {}
    for gene, vals in merged.items():
        v1, v2 = vals[l1], vals[l2]
        if v1[tp_key] != v2[tp_key] or v1[fp_key] != v2[fp_key]:
            filtered[gene] = vals
    return filtered


def write_multicol_minipage_table(merged, vt, l1, l2, out_prefix, ncols=4):
    """
    Writes a LaTeX table and a CSV file.
    """
    tp_key, fp_key = f"{vt}_TP", f"{vt}_FP"
    rows = []
    total = {l1: {"TP": 0, "FP": 0}, l2: {"TP": 0, "FP": 0}}

    for gene, vals in sorted(merged.items()):
        v1, v2 = vals[l1], vals[l2]
        t1, f1, t2, f2 = v1[tp_key], v1[fp_key], v2[tp_key], v2[fp_key]
        rows.append([gene, t1, f1, t2, f2])
        total[l1]["TP"] += t1
        total[l1]["FP"] += f1
        total[l2]["TP"] += t2
        total[l2]["FP"] += f2

    if not rows:
        print(f"[⚠] No differing {vt} counts between {l1} and {l2}")
        return

    # Add totals row
    rows.append([
        "Total",
        total[l1]["TP"], total[l1]["FP"],
        total[l2]["TP"], total[l2]["FP"],
    ])

    # Save CSV
    csv_file = f"{out_prefix}_{vt}.csv"
    with open(csv_file, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Gene", f"TP_{l1}", f"FP_{l1}", f"TP_{l2}", f"FP_{l2}"])
        writer.writerows(rows)

    # LaTeX table
    tex_file = f"{out_prefix}_{vt}.tex"
    nrows_per_col = math.ceil(len(rows)/ncols)

    with open(tex_file, "w") as f:
        f.write("% Auto-generated multi-column landscape table with minipages\n")
        f.write("\\begin{landscape}\n")
        f.write("\\begin{table}[htbp]\n\\centering\n")
        f.write("\\resizebox{\\textwidth}{!}{%\n")
        
        for col in range(ncols):
            start = col * nrows_per_col
            end = start + nrows_per_col
            chunk = rows[start:end]
            f.write("\\begin{minipage}{0.32\\textwidth}\n")
            f.write("\\begin{tabular}{lrrrr}\n\\hline\n")
            f.write(f" & \\multicolumn{{2}}{{c}}{{{l1}}} & \\multicolumn{{2}}{{c}}{{{l2}}} \\\\\n")
            f.write("Gene & TP & FP & TP & FP \\\\\n\\hline\n")
            for r in chunk:
                f.write(f"{r[0]} & {r[1]} & {r[2]} & {r[3]} & {r[4]} \\\\\n")
            f.write("\\hline\n\\end{tabular}\n")
            f.write("\\end{minipage}%\n")
            if col < ncols-1:
                f.write("\\hspace{0.5cm}%\n")  # space between columns
        f.write("}\n")  # end resizebox
        f.write(f"\\caption{{Per-gene {vt} differences between {l1} and {l2}}}\n")
        f.write(f"\\label{{tab:{out_prefix.lower()}_{vt.lower()}}}\n")
        f.write("\\end{table}\n\\end{landscape}\n")

    print(f"[✔] Saved CSV: {csv_file}, LaTeX: {tex_file}")


def main():
    ap = argparse.ArgumentParser(description="Generate landscape tables with side-by-side columns")
    ap.add_argument("--vcf1", required=True)
    ap.add_argument("--vcf2", required=True)
    ap.add_argument("--bed", required=True)
    ap.add_argument("--label1", default="VCF1")
    ap.add_argument("--label2", default="VCF2")
    ap.add_argument("--out-prefix", default="results")
    ap.add_argument("--ncols", type=int, default=2, help="Number of side-by-side columns")
    args = ap.parse_args()

    gene_bed = pybedtools.BedTool(args.bed)
    c1 = count_by_gene(args.vcf1, gene_bed, args.label1)
    c2 = count_by_gene(args.vcf2, gene_bed, args.label2)
    merged = merge_counts(c1, c2, args.label1, args.label2)

    for vt in ["SNV", "Indel"]:
        write_csv_tables(merged, vt, args.label1, args.label2, args.out_prefix)
        filtered = filter_differences(merged, vt, args.label1, args.label2)
        write_multicol_minipage_table(filtered, vt, args.label1, args.label2, args.out_prefix, args.ncols)

if __name__ == "__main__":
    main()
