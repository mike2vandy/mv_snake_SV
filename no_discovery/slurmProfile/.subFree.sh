#! /bin/bash -l
#SBATCH -p msismall
#SBATCH --time=96:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mwvandew@ncsu.edu

cd $SLURM_SUBMIT_DIR

conda activate snakemake 

snakemake -s meltMain.3.smk \
	--use-conda \
	--keep-going \
        --rerun-incomplete \
        --profile slurmProfile 
        
