# Evaluation of the Ultima Genomics UG 100™ sequencing platform for clinical services
The steep reduction in the cost of genome sequencing started with the introduction of Solexa Sequencing-By-Synthesis technology in 2006 has plateaued recently due to technical limitations in the use of closed flowcells and to the cost of reagents. The Ultima Genomics UG 100™ is the first sequencing machine to lower the cost of human genome sequencing to 80\$. However, technical limitations in resolving long homopolymer regions undermine the application of this technology to short variant calling in clinical settings. Here, we evaluate the ability of UG 100™ to identify short variants in relation to its suitability for clinical accreditation, by comparing it with the Illumina NovaSeq 6000 Systems platform. We focus specifically on the small variant calling performance in long homopolymer regions, both genome-wide and in relation to a set of medically-relevant genes that are challenging to sequence. Our analysis aims at supporting clinicians in determining whether the UG 100™ platform is well-suited for their studies, and to guide clinical sequencing centers in evaluating the adoption of this emerging technology.

This is the repository for the code to reproduce the analysis of the paper:

[**Evaluation of the Ultima Genomics UG 100™ sequencing platform for clinical services**](https://www.biorxiv.org/content/10.64898/2025.12.11.693628v1)  

Authors: Santuari Luca, Kolpakov Ilya, Anthony Blin, Ioannis Xenarios and Cédric Howald

*bioRxiv* (2025)  

This work was carried out at the [Health 2030 Genome Center](https://www.health2030genome.ch/), [Campus Biotech](https://campusbiotech.ch/), Geneva, Switzerland.


## Setup
As a package manager, we use [conda](https://anaconda.org/channels/anaconda/packages/conda/overview).

Install [miniforge](https://github.com/conda-forge/miniforge) with:
```
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh
```

Create and activate the conda environment as follows. For instance for the environment in `analyses/clinvar_overlap`, to run the R script: `analyses/clinvar_overlap/find_overlap.R`:
```
mamba env create -f env/renv.yaml
mamba activate clinvar_overlap
```

Figures and tables can be reproduced:
- Figures 1,2,3 and 5 with the Jupyter notebooks:
    - GIAB pilot study (Figure 1 and 2): `notebooks/HG2030GC_GIAB.ipynb`;
    - CMRG HG002 benchmark (Figure 3): `notebooks/CMRG.ipynb`;
    - UG GIAB dataset (Figure 5): `notebooks/HG2030GC_GIAB.ipynb`;
- Figure 4 with the scripts in 'analyses/cmrg_table';
- Tables 1 and 2 can be reproduced with the R script `analyses/clinvar_overlap/find_overlap.R` and the Jupyter notebook `analyses/clinvar_overlap/generate_latex_table.ipynb`.

Note: ChatGPT was used to assist with streamlining the code.

Feel free to open an issue if needed.

## Citation

If you use this repository, please cite our preprint:

```bibtex
@article {Santuari2025.12.11.693628,
	author = {Santuari, Luca and Kolpakov, Ilya and Blin, Anthony and Xenarios, Ioannis and Howald, Cedric},
	title = {Evaluation of the Ultima Genomics UG 100{\texttrademark} sequencing platform for clinical services},
	elocation-id = {2025.12.11.693628},
	year = {2025},
	doi = {10.64898/2025.12.11.693628},
	publisher = {Cold Spring Harbor Laboratory},
	abstract = {The steep reduction in the cost of genome sequencing started with the introduction of Solexa Sequencing-By-Synthesis technology in 2006 has plateaued recently due to technical limitations in the use of closed flowcells and to the cost of reagents. The Ultima Genomics UG 100{\texttrademark} is the first sequencing machine to lower the cost of human genome sequencing to 80$. However, technical limitations in resolving long homopolymer regions undermine the application of this technology to short variant calling in clinical settings. Here, we evaluate the ability of UG 100{\texttrademark} to identify short variants in relation to its suitability for clinical accreditation, by comparing it with the Illumina NovaSeq 6000 Systems platform. We focus specifically on the small variant calling performance in long homopolymer regions, both genome-wide and in relation to a set of medically-relevant genes that are challenging to sequence. Our analysis aims at supporting clinicians in determining whether the UG 100{\texttrademark} platform is well-suited for their studies, and to guide clinical sequencing centers in evaluating the adoption of this emerging technology.Competing Interest StatementThe authors have declared no competing interest.},
	URL = {https://www.biorxiv.org/content/early/2025/12/12/2025.12.11.693628},
	eprint = {https://www.biorxiv.org/content/early/2025/12/12/2025.12.11.693628.full.pdf},
	journal = {bioRxiv}
}
