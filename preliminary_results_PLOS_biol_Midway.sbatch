#!/bin/bash
#SBATCH --job-name=prelim
#SBATCH --time=XXX
#SBATCH --output=slurm_output/prelim%A_%a.out
#SBATCH --error=slurm_output/prelim%A_%a.err
##SBATCH --array=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=1 
#SBATCH --mem-per-cpu=XXX
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