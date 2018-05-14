library(tidyverse)
setwd('~/Documents/malaria_interventions')

set_parameter <- function(param_data, parameter, value){
  param_data$value[param_data$param==parameter] <- value
  return(subset(param_data, param==parameter))
}


# Get reference parameter file --------------------------------------------
reference <- readLines('parameter_file_ref.py')
lines <- 1:length(reference)
parameters <- map(lines, function(l){
  tmp <- str_sub(reference[l], 1, str_locate(reference[l], '=')[1]-1)
  if(!is.na(tmp)){return(tmp)}
}) %>% combine()

param_values <- map(lines, function(l){
  tmp <- str_trim(str_sub(reference[l], str_locate(reference[l], '=')[1]+1, 1000), side = 'both')
  if(!is.na(tmp)){return(tmp)}
}) %>% combine()

param_data <- data.frame(param=str_trim(parameters), value=param_values, stringsAsFactors = F)


# Define parameters for a given parameter set -----------------------------
parameter_space <- '01'
scenario <- 'S'
experiment <- '01'
BITING_RATE_MEAN <- 1
N_GENES_INITIAL <- 1200 # The initial gene pool

# The biting rate distribution for seasonality and/or IRS is generated in
# Mathematica and exported to a file. This distribution consists of the number
# of adult mosquitos. This is done per experiment.
mathematica_file <- 'mosquito_population_seasonality.csv' 
biting_rate_mathematica <- read_csv(mathematica_file, col_names = c('day','num_mosquitos'))
DAILY_BITING_RATE_DISTRIBUTION <- biting_rate_mathematica$num_mosquitos
# Or, if using a distribution not from a file
DAILY_BITING_RATE_DISTRIBUTION <- paste(rep(0.03,360),collapse=',') # This is to set a fixed biting rate. Mostly for testing.
T_END <- 2880*2

# Set parameters ----------------------------------------------------------
# Genetic diversity
param_data[param_data$param=='N_GENES_INITIAL',] <- set_parameter(param_data, 'N_GENES_INITIAL', N_GENES_INITIAL)
param_data[param_data$param=='N_ALLELES_INITIAL',] <- set_parameter(param_data, 'N_ALLELES_INITIAL', paste('N_LOCI*[',N_GENES_INITIAL/10,']',sep=''))
# Biting rates
param_data[param_data$param=='BITING_RATE_MEAN',] <- set_parameter(param_data, 'BITING_RATE_MEAN', paste('[',BITING_RATE_MEAN,']',sep=''))
param_data[param_data$param=='DAILY_BITING_RATE_DISTRIBUTION',] <- set_parameter(param_data, 'DAILY_BITING_RATE_DISTRIBUTION', paste('[',paste(DAILY_BITING_RATE_DISTRIBUTION, collapse=','),']',sep=''))
# Timing
param_data[param_data$param=='T_END',] <- set_parameter(param_data, 'T_END', T_END)
param_data[param_data$param=='VERIFICATION_PERIOD',] <- set_parameter(param_data, 'VERIFICATION_PERIOD', T_END)


# Create checkpoints ------------------------------------------------------
# This section prepares the parameter file to create a checkpoint.
param_data[param_data$param=='SAMPLE_DB_FILENAME',] <- set_parameter(param_data, 'SAMPLE_DB_FILENAME', paste('\'\"PS',parameter_space,'_',scenario,'.sqlite\"\'',sep='')) # The run ID will be added while running the job (in the sbatch execution).
param_data[param_data$param=='SAVE_TO_CHECKPOINT',] <- set_parameter(param_data, 'SAVE_TO_CHECKPOINT', 'True')
param_data[param_data$param=='CHECKPOINT_SAVE_PERIOD',] <- set_parameter(param_data, 'CHECKPOINT_SAVE_PERIOD', T_END)
param_data[param_data$param=='CHECKPOINT_SAVE_FILENAME',] <- set_parameter(param_data, 'CHECKPOINT_SAVE_FILENAME', paste('\'\"PS',parameter_space,'_',scenario,'_CP.sqlite\"\'',sep=''))
param_data[param_data$param=='LOAD_FROM_CHECKPOINT',] <- set_parameter(param_data, 'LOAD_FROM_CHECKPOINT', 'False')
param_data[param_data$param=='CHECKPOINT_LOAD_FILENAME',] <- set_parameter(param_data, 'CHECKPOINT_LOAD_FILENAME', '\'\"\"\'')
# Name of output parameter file
output_file=paste('PS',parameter_space,'_',scenario,'.py',sep = '')

# Load from a checkpoint --------------------------------------------------
run <- 1 # This is the particular run of the parameter space which should be loaded
# This section prepares a parameter file to load a checkpoint and run an experiment
param_data[param_data$param=='SAMPLE_DB_FILENAME',] <- set_parameter(param_data, 'SAMPLE_DB_FILENAME', paste('\'\"PS',parameter_space,'_',scenario,'_E',experiment,'.sqlite\"\'',sep='')) # The run ID will be added while running the job (in the sbatch execution).
param_data[param_data$param=='SAVE_TO_CHECKPOINT',] <- set_parameter(param_data, 'SAVE_TO_CHECKPOINT', 'False')
param_data[param_data$param=='CHECKPOINT_SAVE_PERIOD',] <- set_parameter(param_data, 'CHECKPOINT_SAVE_PERIOD', 0)
param_data[param_data$param=='CHECKPOINT_SAVE_FILENAME',] <- set_parameter(param_data, 'CHECKPOINT_SAVE_FILENAME', '\'\"\"\'')
param_data[param_data$param=='LOAD_FROM_CHECKPOINT',] <- set_parameter(param_data, 'LOAD_FROM_CHECKPOINT', 'True')
param_data[param_data$param=='CHECKPOINT_LOAD_FILENAME',] <- set_parameter(param_data, 'CHECKPOINT_LOAD_FILENAME', paste('\'\"PS',parameter_space,'_',scenario,'_R',run,'_CP.sqlite\"\'',sep=''))
param_data[param_data$param=='T_BURNIN',] <- set_parameter(param_data, 'T_BURNIN', 2880)
# Name of output parameter file
output_file=paste('PS',parameter_space,'_',scenario,'_E',experiment,'.py',sep = '')

# Write to a new parameter file --------------------------------------------
param_data$output <- paste(param_data$param,param_data$value,sep='=')
write_lines(param_data$output, output_file)


# Write Midway job file ---------------------------------------------------
job_lines <- readLines('job_file_ref.sbatch')
wall_time <- '00:15:00'
SLURM_ARRAY_RANGE <- '1'
memory <- '3000'
CP_state <- 'load' # Can be 'create', 'load', or 'none'
if (CP_state=='create'){
  output_file=paste('PS',parameter_space,'_',scenario,'.sbatch',sep = '')
  job_lines[47] <- "cd 'PS'$PS'_'$scenario'_R'$run'"
}
if (CP_state=='load'){
  output_file=paste('PS',parameter_space,'_',scenario,'_E',experiment,'.sbatch',sep = '')
  job_lines[24] <- "work_folder=$base_folder'PS'$PS'_'$scenario'_E'$exp'/run_'$run"
  job_lines[31] <- "mkdir -p 'PS'$PS'_'$scenario'_E'$exp # create the paramter space folder, if does not exist"
  job_lines[33] <- "cp 'PS'$PS'_'$scenario'_E'$exp'.py' $work_folder # Copy the parameter space parameter file"
  job_lines[40] <- "./build.py -p 'PS'$PS'_'$scenario'_E'$exp'.py' -d 'PS'$PS'_'$scenario'_R'$run'_E'$exp"
  job_lines[47] <- "cd 'PS'$PS'_'$scenario'_R'$run'_E'$exp"
}

job_lines[2] <- paste('#SBATCH --job-name=PS',parameter_space,sep='')
job_lines[3] <- paste('#SBATCH --time=',wall_time,sep='')
job_lines[4] <- paste('#SBATCH --output=slurm_output/PS',parameter_space,'_%A_%a.out',sep='')
job_lines[5] <- paste('#SBATCH --error=slurm_output/PS',parameter_space,'_%A_%a.err',sep='')
job_lines[6] <- paste('#SBATCH --array=',SLURM_ARRAY_RANGE,sep='')
job_lines[9] <- paste('#SBATCH --mem-per-cpu=',memory,sep='')
job_lines[19] <- paste("PS='",parameter_space,"'",sep='')
job_lines[20] <- paste("scenario='",scenario,"'",sep='')
job_lines[22] <- paste("exp='",experiment,"'",sep='')
job_lines[26] <- paste("CHECKPOINT='",CP_state,"'",sep='')
write_lines(job_lines, output_file)

# plots -------------------------------------------------------------------
seasonality <- read_csv('mosquito_population_seasonality.csv', col_names = c('day','num_mosquitos'))
seasonality$grp <- 'S'
irs01 <- read_csv('mosquito_population_IRS01.csv', col_names = c('day','num_mosquitos'))
irs01$grp <- 'IRS01'
irs02 <- read_csv('mosquito_population_IRS02.csv', col_names = c('day','num_mosquitos'))
irs02$grp <- 'IRS02'
x <- rbind(seasonality,irs01,irs02)
x %>% ggplot(aes(day, num_mosquitos, color=grp))+geom_line()


# Checkpoints -------------------------------------------------------------

