#!/bin/bash
#SBATCH --job-name=PS01
#SBATCH --time=00:15:00
#SBATCH --output=slurm_output/PS01_%A_%a.out
#SBATCH --error=slurm_output/PS01_%A_%a.err
#SBATCH --array=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=1 
#SBATCH --mem-per-cpu=3000
#SBATCH --partition=broadwl
##SBATCH --mail-type=END
##SBATCH --mail-user=pilosofs@uchicago.edu

# Load modules
module load python/2.7
module load gcc/4.9

# Initialize
PS='01'
scenario='S'
run=$SLURM_ARRAY_TASK_ID # The run in the parameter space
exp='XXX' # The experiment ID
base_folder='/home/pilosofs/malaria_interventions/' # Where all the source files are and where the experiment folders will be created
work_folder=$base_folder'PS'$PS'/run_'$run
ABM_file=varmodel2.zip
CHECKPOINT='create'

# Set up the experiment
cd $base_folder # Go to base folder
rm -rf $work_folder # Remove folder if exists
mkdir -p 'PS'$PS # create the paramter space folder, if does not exist
mkdir $work_folder
cp 'PS'$PS'_'$scenario'.py' $work_folder # Copy the parameter space parameter file
cp $ABM_file $work_folder # Copy the ABM code
cd $work_folder
unzip $ABM_file # Unzip ABM code

# Run the ABM
./build.py -p 'PS'$PS'_'$scenario'.py' -d 'PS'$PS'_'$scenario'_R'$run
./'PS'$PS'_'$scenario'_R'$run/bin/varmodel2

# Rename sqlite files
mv 'PS'$PS'_'$scenario'.sqlite' 'PS'$PS'_'$scenario'_R'$run'.sqlite'
if [ "$CHECKPOINT" = "create" ]; then
	mv 'PS'$PS'_'$scenario'_CP.sqlite' 'PS'$PS'_'$scenario'_R'$run'_CP.sqlite'
fi

# Copy sqlite files to the sqlite folder (include checkpoint files)
cp *.sqlite $base_folder'/sqlite'

# Remove unnecesary files
rm -r sqlite3
rm -r src
rm build.py
rm generate_managers.py
rm generate_parameters.py
rm README.md


#if [ "$CHECKPOINT" = "load" ]; then
	#cd $base_folder'/sqlite'
	#cp $exp'_checkpoint.sqlite' $exp'_'$run'_checkpoint.sqlite'
#fi
