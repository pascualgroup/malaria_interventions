#!/bin/bash
#SBATCH --job-name=prelim
#SBATCH --time=36:00:00
#SBATCH --output=slurm_output/prelim_%A.out
#SBATCH --error=slurm_output/prelim_%A.err
##SBATCH --array=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=1 
#SBATCH --mem-per-cpu=58000
#SBATCH --partition=broadwl
##SBATCH --mail-type=END
##SBATCH --mail-user=pilosofs@uchicago.edu

# Load modules
module load R
module load gcc/6.1

# Initialize
base_folder='/scratch/midway2/pilosofs/PLOS_Biol/' # Where all the source files are and where the experiment folders will be created

cd $base_folder # Go to base folder
Rscript preliminary_results_PLOS_biol_Midway.R