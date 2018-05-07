#!/bin/bash
#SBATCH --job-name=test_01
#SBATCH --time=05:00:00
#SBATCH --output=slurm_output/test_01_%A_%a.out
#SBATCH --error=slurm_output/test_01_%A_%a.err
#SBATCH --array=1-2
#SBATCH --tasks=1
#SBATCH --cpus-per-task=1 
#SBATCH --mem-per-cpu=10000
#SBATCH --partition=broadwl
##SBATCH --mail-type=END
##SBATCH --mail-user=pilosofs@uchicago.edu

# Load modules
module load python/2.7
module load gcc/4.9

# Initialize
run=$SLURM_ARRAY_TASK_ID
exp='test_01'
base_folder='/home/pilosofs/varmodel2/' # Where all the source files are and where the experiment folders will be created
work_folder=$base_folder'/'$exp'/run_'$run
ABM_file=varmodel2.zip

# Set up the experiment
cd $base_folder # Go to base folder
rm -rf $work_folder # Remove folder if exists
mkdir -p $exp # create the experiment folder, if does not exist
mkdir $work_folder
cp $exp'.py' $work_folder # Copy the experiment parameter file
cp $ABM_file $work_folder # Copy the ABM code
cd $work_folder
unzip $ABM_file # Unzip ABM code

# Run the ABM
./build.py -p $exp'.py' -d $exp'_'$run
./$exp'_'$run/bin/varmodel2

# Rename sqlite
mv $exp'.sqlite' $exp'_'$run'.sqlite'

# Copy sqlite to the sqlite folder
cp $exp'_'$run'.sqlite' $base_folder'/sqlite'

# Remove files
rm -r sqlite3
rm -r src
rm build.py
rm generate_managers.py
rm generate_parameters.py
rm README.md
