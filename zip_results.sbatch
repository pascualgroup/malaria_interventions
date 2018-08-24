#!/bin/bash
#SBATCH --job-name=zipresults
#SBATCH --time=04:00:00
#SBATCH --output=slurm_output/zipresults_%A_%a.out
#SBATCH --error=slurm_output/zipresults_%A_%a.err
#SBATCH --array=27-33
#SBATCH --tasks=1
#SBATCH --cpus-per-task=1 
#SBATCH --mem-per-cpu=16000
#SBATCH --partition=broadwl
##SBATCH --mail-type=END
##SBATCH --mail-user=pilosofs@uchicago.edu

# Initialize
PS=$SLURM_ARRAY_TASK_ID # the parameter space

cd '/scratch/midway2/pilosofs/malaria_interventions/Results'
for scenario in 'S' 'N' 'G'
do
  if [[ -d $PS'_'$scenario ]];
  then
  zip -r -m $PS'_'$scenario'_summary_general.zip' $PS'_'$scenario'/' -i '*summary_general*.*'
  zip -r -m $PS'_'$scenario'_sampled_infections.zip' $PS'_'$scenario'/' -i '*sampled_infections*.*'
  zip -r -m $PS'_'$scenario'_node_list.zip' $PS'_'$scenario'/' -i '*node_list*.*'
  zip -r -m $PS'_'$scenario'_layer_summary.zip' $PS'_'$scenario'/' -i '*layer_summary*.*'
  zip -r -m $PS'_'$scenario'_network_info.zip' $PS'_'$scenario'/' -i '*network_info*.*'
  zip -r -m $PS'_'$scenario'_intralayer.zip' $PS'_'$scenario'/' -i '*intralayer*.*'
  zip -r -m $PS'_'$scenario'_interlayer.zip' $PS'_'$scenario'/' -i '*interlayer*.*'
  zip -r -m $PS'_'$scenario'_Infomap_multilayer.zip' $PS'_'$scenario'/' -i '*Infomap_multilayer*.*'
  fi
done
