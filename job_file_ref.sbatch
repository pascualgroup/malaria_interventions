#!/bin/bash
#SBATCH --job-name=XXX
#SBATCH --time=20:00:00
#SBATCH --output=slurm_output/XXX_%A_%a.out
#SBATCH --error=slurm_output/XXX_%A_%a.err
#SBATCH --array=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=1 
#SBATCH --mem-per-cpu=20000
#SBATCH --partition=broadwl
##SBATCH --mail-type=END
##SBATCH --mail-user=pilosofs@uchicago.edu

# Load modules
module load python/2.7
module load gcc/4.9

# Initialize
PS = 'XXX' # the parameter space
scenario = 'XXX' # Can be: 'S', 'G' or 'N'
run=$SLURM_ARRAY_TASK_ID # The run in the parameter space
exp='XXX' # The experiment ID
base_folder='/home/pilosofs/malaria_interventions/' # Where all the source files are and where the experiment folders will be created
$base_name='PS'$PS'_'$scenario'_E'$exp'_R'$run
work_folder=$base_folder'PS'$PS'_'$scenario'_E'$exp'/run_'$run
ABM_file=varmodel2.zip
CHECKPOINT='none' # Can be: 'create', 'load', or 'none'

# Set up the experiment
cd $base_folder # Go to base folder
rm -rf $work_folder # Remove folder if exists
mkdir -p 'PS'$PS'_'$scenario'_E'$exp # create the paramter space folder, if does not exist
mkdir $work_folder
cp 'PS'$PS'_'$scenario'_E'$exp'.py' $work_folder # Copy the parameter space parameter file
cp $ABM_file $work_folder # Copy the ABM code
cd $work_folder
unzip -qq $ABM_file # Unzip ABM code

# Compile ABM
echo Paramter file is: 'PS'$PS'_'$scenario'_E'$exp'.py'
./build.py -p 'PS'$PS'_'$scenario'_E'$exp'.py' -d 'PS'$PS'_'$scenario'_R'$run'_E'$exp

if [ "$CHECKPOINT" = "load" ]; then
	cp $base_folder'sqlite/PS'$PS'_'$scenario'_R'$run'_CP.sqlite' $work_folder'/PS'$PS'_'$scenario'_R'$run'_E'$exp
	# Run the ABM
	cd 'PS'$PS'_'$scenario'_R'$run'_E'$exp
	./bin/varmodel2
	# Rename sqlite files
	mv 'PS'$PS'_'$scenario'_E'$exp'.sqlite' 'PS'$PS'_'$scenario'_E'$exp'_R'$run'.sqlite'
	cp 'PS'$PS'_'$scenario'_E'$exp'_R'$run'.sqlite' $base_folder'/sqlite'

fi


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
