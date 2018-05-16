library(tidyverse)
setwd('~/Documents/malaria_interventions')

# Function to set a parameter
set_parameter <- function(param_data, parameter, value){
  param_data$value[param_data$param==parameter] <- value
  return(subset(param_data, param==parameter))
}

# Function to get reference parameter file
get_parameter_reference <- function(parameter_file_ref='parameter_file_ref.py'){
  reference <- readLines(parameter_file_ref)
  lines <- 1:length(reference)
  parameters <- map(lines, function(l){
    tmp <- str_sub(reference[l], 1, str_locate(reference[l], '=')[1]-1)
    if(!is.na(tmp)){return(tmp)}
  }) %>% combine()
  
  param_values <- map(lines, function(l){
    tmp <- str_trim(str_sub(reference[l], str_locate(reference[l], '=')[1]+1, 1000), side = 'both')
    if(!is.na(tmp)){return(tmp)}
  }) %>% combine()
  
  return(data.frame(param=str_trim(parameters), value=param_values, stringsAsFactors = F))
}

loadExperiments_GoogleSheets <- function(workBookName='malaria_interventions_design',sheetID=4){
  require(googlesheets)
  GS <- gs_title(workBookName)
  col_types <- GS %>% gs_read(ws=1, col_names=T)
  col_t <- unname(as.list(col_types[1,]))
  experiments <- GS %>% gs_read(ws=sheetID, col_names=T, col_types=col_t)
  print(experiments)
  return(experiments)
}

# Get data design ---------------------------------------------------------
design <- loadExperiments_GoogleSheets()

# Function to create the necesary files and pipeline for a single run of an experiment.
# Each run has its own random seed across experiments.
create_run <- function(design_ID, run, RANDOM_SEED){
  # Regime
  parameter_space <- design$PS[design_ID]
  scenario <- design$Scenario[design_ID]
  experiment <- design$Experiment[design_ID] # Use 00 for checkpoints and control
  base_name <- paste('PS',parameter_space,'_',scenario,'_E',experiment,'_R',run,sep='')
  param_data <- get_parameter_reference()
  
  # General parameters
  param_data[param_data$param=='RANDOM_SEED',] <- set_parameter(param_data, 'RANDOM_SEED', RANDOM_SEED)
  T_END <- design$T_END[design_ID]
  param_data[param_data$param=='T_END',] <- set_parameter(param_data, 'T_END', T_END)
  param_data[param_data$param=='VERIFICATION_ON',] <- set_parameter(param_data, 'VERIFICATION_ON', 'False')
  param_data[param_data$param=='VERIFICATION_PERIOD',] <- set_parameter(param_data, 'VERIFICATION_PERIOD', T_END)
  
  # Biting rates
  BITING_RATE_MEAN <- design$BITING_RATE_MEAN[design_ID]
  mathematica_file <- design$DAILY_BITING_RATE_DISTRIBUTION[design_ID]
  biting_rate_mathematica <- read_csv(mathematica_file, col_names = c('day','num_mosquitos'))
  DAILY_BITING_RATE_DISTRIBUTION <- biting_rate_mathematica$num_mosquitos
  param_data[param_data$param=='BITING_RATE_MEAN',] <- set_parameter(param_data, 'BITING_RATE_MEAN', paste('[',BITING_RATE_MEAN,']',sep=''))
  param_data[param_data$param=='DAILY_BITING_RATE_DISTRIBUTION',] <- set_parameter(param_data, 'DAILY_BITING_RATE_DISTRIBUTION', paste('[',paste(DAILY_BITING_RATE_DISTRIBUTION, collapse=','),']',sep=''))
  
  # Genetic diversity
  N_GENES_INITIAL <- design$N_GENES_INITIAL[design_ID]
  param_data[param_data$param=='N_GENES_INITIAL',] <- set_parameter(param_data, 'N_GENES_INITIAL', N_GENES_INITIAL)
  N_ALLELES_INITIAL <- design$N_ALLELES_INITIAL[design_ID]
  param_data[param_data$param=='N_ALLELES_INITIAL',] <- set_parameter(param_data, 'N_ALLELES_INITIAL', paste('N_LOCI*[',N_ALLELES_INITIAL,']',sep=''))

  # Checkpoints
  T_BURNIN <- design$T_BURNIN[design_ID]
  if (experiment == '00'){
    param_data[param_data$param=='SAMPLE_DB_FILENAME',] <- set_parameter(param_data, 'SAMPLE_DB_FILENAME', paste('\'\"',base_name,'.sqlite\"\'',sep='')) # The run ID will be added while running the job (in the sbatch execution).
    param_data[param_data$param=='SAVE_TO_CHECKPOINT',] <- set_parameter(param_data, 'SAVE_TO_CHECKPOINT', 'True')
    param_data[param_data$param=='CHECKPOINT_SAVE_PERIOD',] <- set_parameter(param_data, 'CHECKPOINT_SAVE_PERIOD', T_END) # The save period should be the T_END
    param_data[param_data$param=='CHECKPOINT_SAVE_FILENAME',] <- set_parameter(param_data, 'CHECKPOINT_SAVE_FILENAME', paste('\'\"',base_name,'_CP.sqlite\"\'',sep=''))
    param_data[param_data$param=='LOAD_FROM_CHECKPOINT',] <- set_parameter(param_data, 'LOAD_FROM_CHECKPOINT', 'False')
    param_data[param_data$param=='CHECKPOINT_LOAD_FILENAME',] <- set_parameter(param_data, 'CHECKPOINT_LOAD_FILENAME', '\'\"\"\'')
    param_data[param_data$param=='T_BURNIN',] <- set_parameter(param_data, 'T_BURNIN', T_BURNIN)
  }
  if (experiment != '00'){ # This section prepares a parameter file to load a checkpoint and run an experiment
    param_data[param_data$param=='SAMPLE_DB_FILENAME',] <- set_parameter(param_data, 'SAMPLE_DB_FILENAME', paste('\'\"',base_name,'.sqlite\"\'',sep='')) # The run ID will be added while running the job (in the sbatch execution).
    param_data[param_data$param=='SAVE_TO_CHECKPOINT',] <- set_parameter(param_data, 'SAVE_TO_CHECKPOINT', 'False')
    param_data[param_data$param=='CHECKPOINT_SAVE_PERIOD',] <- set_parameter(param_data, 'CHECKPOINT_SAVE_PERIOD', 0)
    param_data[param_data$param=='CHECKPOINT_SAVE_FILENAME',] <- set_parameter(param_data, 'CHECKPOINT_SAVE_FILENAME', '\'\"\"\'')
    param_data[param_data$param=='LOAD_FROM_CHECKPOINT',] <- set_parameter(param_data, 'LOAD_FROM_CHECKPOINT', 'True')
    param_data[param_data$param=='CHECKPOINT_LOAD_FILENAME',] <- set_parameter(param_data, 'CHECKPOINT_LOAD_FILENAME', paste('\'\"PS',parameter_space,'_',scenario,'_E00','_R',run,'_CP.sqlite\"\'',sep=''))
    param_data[param_data$param=='T_BURNIN',] <- set_parameter(param_data, 'T_BURNIN', T_BURNIN) # The burnin value should be the value where the checkpoint was taken
  }
  
  # Write parameter file
  output_file=paste(base_name,'.py',sep = '')
  param_data$output <- paste(param_data$param,param_data$value,sep='=')
  write_lines(param_data$output, output_file)
}

# This creates several run WITH THE SAME RANDOM SEED for each experiment
for (RUN in 1:2){
  RANDOM_SEED <- round(runif(n = 1, min=1, max=10^7),0)
  for (id in 4:6){
    create_run(design_ID = id, RUN, RANDOM_SEED)
  }
}

SLURM_ARRAY_RANGE <- '1-2'
for (design_ID in 4:6){
  parameter_space <- design$PS[design_ID]
  scenario <- design$Scenario[design_ID]
  experiment <- design$Experiment[design_ID] # Use 00 for checkpoints and control
  base_name <- paste('PS',parameter_space,scenario,'E',experiment,sep='')
  # Write the job file for the exepriment
  job_lines <- readLines('job_file_ref.sbatch')
  wall_time <-  design$wall_time[design_ID]
  mem_per_cpu <-  design$mem_per_cpu[design_ID]
  job_lines[2] <- paste('#SBATCH --job-name=',base_name,sep='')
  job_lines[3] <- paste('#SBATCH --time=',wall_time,sep='')
  job_lines[4] <- paste('#SBATCH --output=slurm_output/',base_name,'_%A_%a.out',sep='')
  job_lines[5] <- paste('#SBATCH --error=slurm_output/',base_name,'_%A_%a.err',sep='')
  job_lines[6] <- paste('#SBATCH --array=',SLURM_ARRAY_RANGE,sep='')
  job_lines[9] <- paste('#SBATCH --mem-per-cpu=',mem_per_cpu,sep='')
  job_lines[19] <- paste("PS='",parameter_space,"'",sep='')
  job_lines[20] <- paste("scenario='",scenario,"'",sep='')
  job_lines[21] <- paste("exp='",experiment,"'",sep='')
  if (experiment == '00'){
    job_lines[26] <- "CHECKPOINT='create'"
  }
  if (experiment != '00'){
    job_lines[26] <- "CHECKPOINT='load'"
  }
  output_file=paste(base_name,'.sbatch',sep = '')
  write_lines(job_lines, output_file)
}



















# Define parameters for a given parameter set -----------------------------
parameter_space <- '01'
scenario <- 'S'
experiment <- '00' # Use 00 for checkpoints and control
run <- 1
base_name <- paste('PS',parameter_space,'_',scenario,'_E',experiment,'_R',run,sep='')

RANDOM_SEED <- round(runif(n = 1, min=1, max=10^7),0)
BITING_RATE_MEAN <- 1
N_GENES_INITIAL <- 1200 # The initial gene pool

# The biting rate distribution for seasonality and/or IRS is generated in
# Mathematica and exported to a file. This distribution consists of the number
# of adult mosquitos. This is done per experiment.
mathematica_file <- 'mosquito_population_seasonality.csv' 
biting_rate_mathematica <- read_csv(mathematica_file, col_names = c('day','num_mosquitos'))
DAILY_BITING_RATE_DISTRIBUTION <- biting_rate_mathematica$num_mosquitos
# Or, if using a distribution not from a file
DAILY_BITING_RATE_DISTRIBUTION <- paste(rep(0.01,360),collapse=',') # This is to set a fixed biting rate. Mostly for testing.
T_END <- 5760

# Set parameters ----------------------------------------------------------
# General
param_data[param_data$param=='RANDOM_SEED',] <- set_parameter(param_data, 'RANDOM_SEED', RANDOM_SEED)
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
param_data[param_data$param=='SAMPLE_DB_FILENAME',] <- set_parameter(param_data, 'SAMPLE_DB_FILENAME', paste('\'\"',base_name,'.sqlite\"\'',sep='')) # The run ID will be added while running the job (in the sbatch execution).
param_data[param_data$param=='SAVE_TO_CHECKPOINT',] <- set_parameter(param_data, 'SAVE_TO_CHECKPOINT', 'True')
param_data[param_data$param=='CHECKPOINT_SAVE_PERIOD',] <- set_parameter(param_data, 'CHECKPOINT_SAVE_PERIOD', T_END) # The save period should be the T_END
param_data[param_data$param=='CHECKPOINT_SAVE_FILENAME',] <- set_parameter(param_data, 'CHECKPOINT_SAVE_FILENAME', paste('\'\"',base_name,'_CP.sqlite\"\'',sep=''))
param_data[param_data$param=='LOAD_FROM_CHECKPOINT',] <- set_parameter(param_data, 'LOAD_FROM_CHECKPOINT', 'False')
param_data[param_data$param=='CHECKPOINT_LOAD_FILENAME',] <- set_parameter(param_data, 'CHECKPOINT_LOAD_FILENAME', '\'\"\"\'')
param_data[param_data$param=='T_BURNIN',] <- set_parameter(param_data, 'T_BURNIN', 0)

# Load from a checkpoint --------------------------------------------------
# This section prepares a parameter file to load a checkpoint and run an experiment
param_data[param_data$param=='SAMPLE_DB_FILENAME',] <- set_parameter(param_data, 'SAMPLE_DB_FILENAME', paste('\'\"',base_name,'.sqlite\"\'',sep='')) # The run ID will be added while running the job (in the sbatch execution).
param_data[param_data$param=='SAVE_TO_CHECKPOINT',] <- set_parameter(param_data, 'SAVE_TO_CHECKPOINT', 'False')
param_data[param_data$param=='CHECKPOINT_SAVE_PERIOD',] <- set_parameter(param_data, 'CHECKPOINT_SAVE_PERIOD', 0)
param_data[param_data$param=='CHECKPOINT_SAVE_FILENAME',] <- set_parameter(param_data, 'CHECKPOINT_SAVE_FILENAME', '\'\"\"\'')
param_data[param_data$param=='LOAD_FROM_CHECKPOINT',] <- set_parameter(param_data, 'LOAD_FROM_CHECKPOINT', 'True')
param_data[param_data$param=='CHECKPOINT_LOAD_FILENAME',] <- set_parameter(param_data, 'CHECKPOINT_LOAD_FILENAME', paste('\'\"PS',parameter_space,'_',scenario,'_E00','_R',run,'_CP.sqlite\"\'',sep=''))
param_data[param_data$param=='T_BURNIN',] <- set_parameter(param_data, 'T_BURNIN', 2880) # The burnin value should be the value where the checkpoint was taken


# Write to a new parameter file --------------------------------------------
output_file=paste(base_name,'.py',sep = '')
param_data$output <- paste(param_data$param,param_data$value,sep='=')
write_lines(param_data$output, output_file)


# Write Midway job file ---------------------------------------------------
job_lines <- readLines('job_file_ref.sbatch')
wall_time <- '00:15:00'
SLURM_ARRAY_RANGE <- '1'
memory <- '3000'
CP_state <- 'load' # Can be 'create', 'load', or 'none'
job_lines[2] <- paste('#SBATCH --job-name=',parameter_space,scenario,experiment,'R',run,sep='')
job_lines[3] <- paste('#SBATCH --time=',wall_time,sep='')
job_lines[4] <- paste('#SBATCH --output=slurm_output/PS',parameter_space,scenario,experiment,'R',run,'_%A_%a.out',sep='')
job_lines[5] <- paste('#SBATCH --error=slurm_output/PS',parameter_space,scenario,experiment,'R',run,'_%A_%a.err',sep='')
job_lines[6] <- paste('#SBATCH --array=',SLURM_ARRAY_RANGE,sep='')
job_lines[9] <- paste('#SBATCH --mem-per-cpu=',memory,sep='')
job_lines[19] <- paste("PS='",parameter_space,"'",sep='')
job_lines[20] <- paste("scenario='",scenario,"'",sep='')
job_lines[21] <- paste("exp='",experiment,"'",sep='')
job_lines[26] <- paste("CHECKPOINT='",CP_state,"'",sep='')
output_file=paste(base_name,'.sbatch',sep = '')
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

