This analysis runs hap.py on the Ultima and Illumina dataset of the GIAB pilot study and on the UG GIAB dataset, using the GIAB NIST v.4.2.1 benchmark.

## Setup
Install miniforge:
```
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"

bash Miniforge3-$(uname)-$(uname -m).sh
```

Install and activate the conda environments
```
mamba env create -f env/hydra-env.yaml

mamba env create -f env/happy-env.yaml

mamba activate hydra-env
```

## Analysis
- Run the GIAB pilot study analysis with:

`python run_happy.py config=config-ultima`

- Run the Illumina dataset analysis with:

`python run_happy.py config=config-illumina`

- Run the UG GIAB dataset analysis with:

`python run_happy.py config=config-ug-giab`
