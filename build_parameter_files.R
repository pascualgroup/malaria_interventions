library(tidyverse)
setwd('~/Documents/malaria_interventions')


# General functions -------------------------------------------------------

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
    tmp <- str_trim(str_sub(reference[l], str_locate(reference[l], '=')[1]+1, 10^6), side = 'both')
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


# Setting general parameters ----------------------------------------------

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


# IRS ---------------------------------------------------------------------

# This function sets the parameters for a single IRS. It is possible to run it several times, once for each IRS scheme.
set_IRS <- function(design_ID, run, IRS_START_TIME, IRS_input, IRS_IMMIGRATION){
  # Regime
  parameter_space <- design$PS[design_ID]
  scenario <- design$Scenario[design_ID]
  experiment <- design$Experiment[design_ID] # Use 00 for checkpoints and control
  base_name <- paste('PS',parameter_space,'_',scenario,'_E',experiment,'_R',run,sep='')
  output_file <- paste(base_name,'.py',sep = '')
  param_data <- get_parameter_reference(output_file)
  
  # Turn on IRS
  param_data[param_data$param=='IRS_ON',] <- set_parameter(param_data, 'IRS_ON', 'True')
  
  # Add the start time to the vector IRS_START_TIMES
  x <- substr(param_data[param_data$param=='IRS_START_TIMES',]$value, 1, nchar(param_data[param_data$param=='IRS_START_TIMES',]$value)-1)
  x <- paste(x,',',IRS_START_TIME,']',sep='')
  x <- str_replace(x, '\\[,','\\[')
  param_data[param_data$param=='IRS_START_TIMES',] <- set_parameter(param_data, 'IRS_START_TIMES', x)
  # Add the IRS scheme to the vector BITING_RATE_FACTORS
  IRS_scheme <- read_csv(IRS_input, col_names = c('day','num_mosquitos'))
  x <- substr(param_data[param_data$param=='BITING_RATE_FACTORS',]$value, 1, nchar(param_data[param_data$param=='BITING_RATE_FACTORS',]$value)-1)
  tmp1 <- paste(IRS_scheme$num_mosquitos, collapse=',')
  tmp2 <- paste('[',tmp1,']',sep='')
  x <- paste(x,',',tmp2,']',sep='')
  x <- str_replace(x, '\\[,','\\[')
  param_data[param_data$param=='BITING_RATE_FACTORS',] <- set_parameter(param_data, 'BITING_RATE_FACTORS', x)
  # Add to the IRS_IMMIGRATION_RATE_FACTORS
  x <- substr(param_data[param_data$param=='IRS_IMMIGRATION_RATE_FACTORS',]$value, 1, nchar(param_data[param_data$param=='IRS_IMMIGRATION_RATE_FACTORS',]$value)-1)
  x <- paste(x,',',IRS_IMMIGRATION,']',sep='')
  x <- str_replace(x, '\\[,','\\[')
  param_data[param_data$param=='IRS_IMMIGRATION_RATE_FACTORS',] <- set_parameter(param_data, 'IRS_IMMIGRATION_RATE_FACTORS', x)
  
  # Write parameter file
  param_data$output <- paste(param_data$param,param_data$value,sep='=')
  write_lines(param_data$output, output_file)
}


# MDA ---------------------------------------------------------------------

set_MDA <- function(design_ID, run){
  # Regime
  parameter_space <- design$PS[design_ID]
  scenario <- design$Scenario[design_ID]
  experiment <- design$Experiment[design_ID] # Use 00 for checkpoints and control
  base_name <- paste('PS',parameter_space,'_',scenario,'_E',experiment,'_R',run,sep='')
  output_file <- paste(base_name,'.py',sep = '')
  param_data <- get_parameter_reference(output_file)
  
  # Turn on MDA
  param_data[param_data$param=='MDA_ON',] <- set_parameter(param_data, 'MDA_ON', 'True')
  
  # Set parameters
  param_data[param_data$param=='MDA_START_TIMES',] <- set_parameter(param_data, 'MDA_START_TIMES', paste('[',design$MDA_START[design_ID],']',sep=''))
  param_data[param_data$param=='HOST_FAIL_RATE',] <- set_parameter(param_data, 'HOST_FAIL_RATE', paste('[',design$HOST_FAIL_RATE[design_ID],']',sep=''))
  param_data[param_data$param=='DRUG_EFF_DURATION',] <- set_parameter(param_data, 'DRUG_EFF_DURATION', paste('[',design$DRUG_EFF_DURATION[design_ID],']',sep=''))
  param_data[param_data$param=='MDA_IMMIGRATION_RATE_FACTORS',] <- set_parameter(param_data, 'MDA_IMMIGRATION_RATE_FACTORS', paste('[',design$MDA_IMMIGRATION_RATE_FACTORS[design_ID],']',sep=''))
  
  # Write parameter file
  param_data$output <- paste(param_data$param,param_data$value,sep='=')
  write_lines(param_data$output, output_file)
}

# Create parameter and job files -------------------------------------------
design <- loadExperiments_GoogleSheets() # Get data design 

# This creates several run WITH THE SAME RANDOM SEED for each experiment
for (RUN in 1:1){
  RANDOM_SEED <- round(runif(n = 1, min=1, max=10^7),0)
  RANDOM_SEED <- 9073975
  for (design_ID in 8:10){
    # Create parameter file with main parameters
    create_run(design_ID, RUN, RANDOM_SEED)
    # Set IRS
    if (!is.na(design$IRS_START_TIMES[design_ID])){
      IRS_scheme <- data.frame(IRS_START_TIME=str_split(design$IRS_START_TIMES[design_ID], ',')[[1]],
                               IRS_input=str_split(design$IRS_input[design_ID], ',')[[1]],
                               IRS_IMMIGRATION=str_split(design$IRS_IMMIGRATION[design_ID], ',')[[1]],
                               stringsAsFactors = F)
      for (i in 1:nrow(IRS_scheme)){
        set_IRS(design_ID, RUN, IRS_scheme$IRS_START_TIME[i], IRS_scheme$IRS_input[i], IRS_scheme$IRS_IMMIGRATION[i])
      }
    }
    # Set MDA
    if (!is.na(design$MDA_START[design_ID])){
      set_MDA(design_ID, RUN)
    }
  }
}

SLURM_ARRAY_RANGE <- '1'
for (design_ID in 8:10){
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



# plots -------------------------------------------------------------------
seasonality <- read_csv('mosquito_population_seasonality.csv', col_names = c('day','num_mosquitos'))
seasonality$grp <- 'S'
irs01 <- read_csv('mosquito_population_IRS01.csv', col_names = c('day','num_mosquitos'))
irs01$grp <- 'IRS01'
irs02 <- read_csv('mosquito_population_IRS02.csv', col_names = c('day','num_mosquitos'))
irs02$grp <- 'IRS02'
x <- rbind(seasonality,irs01,irs02)
x %>% ggplot(aes(day, num_mosquitos, color=grp))+geom_line()

