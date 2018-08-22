#!/bin/bash
#SBATCH --job-name=XXX
#SBATCH --time=00:30:00
#SBATCH --output=slurm_output/XXX_%A_%a.out
#SBATCH --error=slurm_output/XXX_%A_%a.err
#SBATCH --array=1-3
#SBATCH --tasks=1
#SBATCH --cpus-per-task=1 
#SBATCH --mem-per-cpu=16000
#SBATCH --partition=broadwl
##SBATCH --mail-type=END
##SBATCH --mail-user=pilosofs@uchicago.edu

# Load modules
module load R/3.3.2
module load gcc/4.9

# Initialize
PS='27' # the parameter space
scenario='S' # Can be: 'S', 'G' or 'N'
exp='001' # The experiment ID
base_folder='/scratch/midway2/pilosofs/malaria_interventions/' # Where all the source files are and where the experiment folders will be created


# Set up
cd $base_folder # Go to base folder
mkdir -p $base_folder'/Results/'$PS'_'$scenario

prog='get_data_Midway.R'
Rscript $prog $PS $scenario $exp
