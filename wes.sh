#!/bin/bash
#SBATCH --job-name=wes_run
#SBATCH -n 64
#SBATCH --mem=128G       # total memory need
#SBATCH --time=72:00:00


##Export up the sentieon license
export SENTIEON_LICENSE=/michorlab/aashna/WES/sentieon/Harvard_DFCI_Iris.lic

##Activate the environment
CONDA_PREFIX=/michorlab/jacobg/miniconda3
export CONDA_ROOT=/michorlab/jacobg/miniconda3
export PATH=/michorlab/jacobg/miniconda3/bin:$PATH
source activate wes
module load samtools
module load bcftools

##Link the Reference files
ln -s /michorlab/aashna/multiomics/wes/reference_files

##Run the snakemake command
snakemake --unlock -s wes/wes.snakefile -j 32 -k
snakemake --rerun-incomplete  -s wes/wes.snakefile -j 32 -k
