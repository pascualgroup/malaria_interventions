#!/bin/bash
#SBATCH --job-name=PS_39
#SBATCH --time=06:00:00
#SBATCH --output=slurm_output/PS_39_%A_%a.out
#SBATCH --error=slurm_output/PS_39_%A_%a.err
#SBATCH --array=1-50
#SBATCH --tasks=1
#SBATCH --cpus-per-task=1 
#SBATCH --mem-per-cpu=32000
#SBATCH --partition=broadwl
##SBATCH --mail-type=END
##SBATCH --mail-user=pilosofs@uchicago.edu

# Load modules
module load R

# Initialize
PS='39' # the parameter space
scenario='S' # Can be: 'S', 'G' or 'N'
#exp='001' # The experiment ID
base_folder='/scratch/midway2/pilosofs/malaria_interventions/' # Where all the source files are and where the experiment folders will be created


# Set up
cd $base_folder # Go to base folder
mkdir -p $base_folder'/Results/'$PS'_'$scenario

prog='get_data_Midway.R'
Rscript $prog $PS $scenario '001'
Rscript $prog $PS $scenario '002'
Rscript $prog $PS $scenario '003'
Rscript $prog $PS $scenario '004'
