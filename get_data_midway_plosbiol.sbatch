#!/bin/bash
#SBATCH --job-name=PSXXX
#SBATCH --time=XXX
#SBATCH --output=slurm_output/PSXXX_%A_%a.out
#SBATCH --error=slurm_output/PSXXX_%A_%a.err
#SBATCH --array=XXX
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
PS='XXX' # the parameter space
scenario='XXX' # Can be: 'S', 'G' or 'N'
#exp='001' # The experiment ID
cutoff_prob=XXX
base_folder='/scratch/midway2/pilosofs/PLOS_Biol/' # Where all the source files are and where the experiment folders will be created


# Set up
prog='get_data_Midway.R'
cd $base_folder # Go to base folder
mkdir -p $base_folder'/Results/'$PS'_'$scenario
cp Infomap_v01926 $base_folder'/Results/'$PS'_'$scenario
cp $prog $base_folder'/Results/'$PS'_'$scenario
cp functions.R $base_folder'/Results/'$PS'_'$scenario
cd $base_folder'/Results/'$PS'_'$scenario


# Run R script
Rscript $prog $PS $scenario '001' $cutoff_prob 'make_networks'


# Run Infomap
infomap_name='PS'$PS'_'$scenario'_E001_R'$SLURM_ARRAY_TASK_ID'_'$cutoff_prob'_Infomap_multilayer'
./Infomap_v01926 $infomap_name'.txt' . -i multilayer -d -N 10 --rawdir --two-level --tree --expanded

# Get infomap results
Rscript $prog $PS $scenario '001' $cutoff_prob 'read_infomap_results'

# remove files
rm functions.R
rm $prog
rm Infomap_v01926
