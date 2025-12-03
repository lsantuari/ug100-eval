In this analysis we compute the overlap of ClinVar small variants (SNVs and Indels) with a unique CLNSIG value (Benign, Likely_begign, Uncertain_significance, Likely_pathogenic and Pathogenic) with respect to the GA4GH homopolymer tracks.

Install miniforge:
```
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh
```

Install and activate the conda environment
```
mamba env create -f env/renv.yaml
mamba activate clinvar_overlap
```

Input data:

1. Download the [clinvar_20250623](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/weekly/clinvar_20250623.vcf.gz) with its [index](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/weekly/clinvar_20250623.vcf.gz.tbi);

2. Download the GA4GH genomic tracks:
    - [GRCh38_SimpleRepeat_homopolymer_4to6_slop5](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.6/GRCh38@all/LowComplexity/GRCh38_SimpleRepeat_homopolymer_4to6_slop5.bed.gz);
    - [GRCh38_SimpleRepeat_homopolymer_7to11_slop5](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.6/GRCh38@all/LowComplexity/GRCh38_SimpleRepeat_homopolymer_7to11_slop5.bed.gz);
    - [GRCh38_SimpleRepeat_homopolymer_gt11_slop5](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.6/GRCh38@all/LowComplexity/GRCh38_SimpleRepeat_homopolymer_gt11_slop5.bed.gz);
    - [GRCh38_SimpleRepeat_homopolymer_gt20_slop5](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.6/GRCh38@all/LowComplexity/GRCh38_SimpleRepeat_homopolymer_gt20_slop5.bed.gz);

3. Adjust the file paths in `config.yml` accordingly;

4. Run the script in the `clinvar_overlap` conda environment with: `Rscript find_overlap.R`