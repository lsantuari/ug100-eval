This analysis generates the table with the per-gene true positives (TP) and false positives (FP) of the ultima and illumina dataset for the CMRG HG002 benchmark.

## Setup
Install miniforge:
```
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh
```

Install and activate the conda environment
```
mamba env create -f env/cmrg-env.yaml

mamba activate cmrg-env
```

## Analysis
The analysis needs the BED file [GRCh38_CMRG_benchmark_gene_coordinates.bed](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh38/SupplementaryFiles/GRCh38_CMRG_benchmark_gene_coordinates.bed)

To generate the table
```
python generate_cmrg_table.py \
    --vcf1 data/HG002_CMRG_Ultima_30X_all_stratifications_v3.6.vcf.gz \
    --vcf2 data/HG002_CMRG_Illumina_happy_all_stratifications_v3.6.vcf.gz \
    --bed data/GRCh38_CMRG_benchmark_gene_coordinates.bed \
    --label1 ultima --label2 illumina
```
This generates the files `results_SNV_full.csv` and `results_Indel_full.csv` needed for the next step.

To generate the plots:
```
python plot_tp_fp_diff.py \
    --snv-csv results_SNV_full.csv \
    --indel-csv results_Indel_full.csv \
    --label1 ultima \
    --label2 illumina \
    --out-prefix CMRG \
    --label-top-n 10 \
    --main-title "CMRG HG002: per-gene ΔFP versus ΔTP"
```