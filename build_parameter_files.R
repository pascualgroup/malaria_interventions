library(tidyverse)
library(sqldf)
setwd('~/Documents/malaria_interventions')


# General functions -------------------------------------------------------

# Function to set a parameter
set_parameter <- function(param_data, parameter, value){
  param_data$value[param_data$param==parameter] <- value
  return(subset(param_data, param==parameter))
}

# Function to get reference parameter file
get_parameter_reference <- function(parameter_file_ref='~/Documents/malaria_interventions/parameter_file_ref.py'){
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

# A function to remove sqlite and parameter files, or any other file for
# specific combinations of parameter space, scenario and experiment. Will remove
# across all runs.
clear_previous_files <- function(parameter_space=NULL, scenario=NULL, experiment=NULL, exclude_sqlite=T, exclude_CP=T, exclute_control=T){
  files <- list.files(path = '~/Documents/malaria_interventions_data/', full.names = T)
  if (!is.null(parameter_space)){
    files <- files[str_detect(files,paste('PS',parameter_space,sep=''))]
  }
  if (!is.null(scenario)){
    files <- files[str_detect(files,paste('_',scenario,'_',sep='')) | str_detect(files,paste(scenario,'E',sep=''))]
  }
  if (!is.null(experiment)){
    files <- files[str_detect(files,paste('E',experiment,sep=''))]
  }
  if(exclude_CP){
    files <- files[!str_detect(files,'E000')]
  }
  if(exclute_control){
    files <- files[!str_detect(files,'E001')]
  }
  if(exclude_sqlite){
    files <- files[!str_detect(files,'\\.sqlite')]
  }
  file.remove(files)
}


# Functions to set parameters ---------------------------------------------

# This function obtains the duration of infection from the selection mode
# counterpart simulations to calculate the var switching rate in the neutral
# scenario. it makes more sense to take the duration from the control
# experiment, otherwise interventions can affect it, so I set exepriment=001 as default.
set_transition_rate_neutral <- function(parameter_space, experiment='001', N_GENES_PER_STRAIN = 60){
  sqlite_files <- list.files(path = '~/Documents/malaria_interventions_data/', pattern=paste('PS',parameter_space,'_S_E',experiment,'_R',sep=''), full.names = T)
  sqlite_files <- sqlite_files[str_detect(sqlite_files,'\\.sqlite')] # only sqlite files
  sqlite_files <- sqlite_files[!str_detect(sqlite_files,'_CP')] # no checkpoint files
  
  doi <- map(sqlite_files, function(f){
    db <- dbConnect(SQLite(), dbname = f)
    sampled_duration <- dbGetQuery(db, 'SELECT duration FROM sampled_duration')
    # print(head(sampled_duration))
    return(sampled_duration)
  }) %>% bind_rows()
  
  mean_doi <- mean(doi$duration)
  TRANSITION_RATE_NOT_IMMUNE <- 1/(mean_doi/N_GENES_PER_STRAIN)
  # setwd('~/Documents/malaria_interventions/')
  return(TRANSITION_RATE_NOT_IMMUNE)
}


#A function to set the parameters for generazlied immunity scenario. Unlike the
#set_transition_rate_neutral, which returns a single value, this function
#returns a list which describes a fit of a curve of a particualr run, so it
#cannot be aggregated across runs.
set_generalized_immunity <- function(parameter_space, experiment='001', run){
require('rPython')
  sqlite_file <- paste('/home/shai/Documents/malaria_interventions_data/','PS',parameter_space,'_S_E',experiment,'_R',run,'.sqlite',sep='')
  pyFile <- readLines('/home/shai/Documents/malaria_interventions/generalized_immunity_fitting.py')
  pyFile[11] <- paste('path=','"',sqlite_file,'"',sep='')
  writeLines(pyFile,'generalized_immunity_fitting_experiment.py')
  python.load('generalized_immunity_fitting_experiment.py')
  GENERAL_IMMUNITY_PARAMS <- python.get('generalImmunityParams')
  GENERAL_IMMUNITY_PARAMS <- paste('[',paste(GENERAL_IMMUNITY_PARAMS, collapse = ','),']',sep='')
  N_INFECTIONS_FOR_GENERAL_IMMUNITY <- python.get('infectionTimesToImmune')
  CLEARANCE_RATE_IMMUNE <- python.get('clearanceRateConstantImmune')
  file.remove('generalized_immunity_fitting_experiment.py')
  return(list(GENERAL_IMMUNITY_PARAMS,N_INFECTIONS_FOR_GENERAL_IMMUNITY,CLEARANCE_RATE_IMMUNE))
}

# This function sets the parameters for a single IRS. It is possible to run it several times, once for each IRS scheme.
set_IRS <- function(design_ID, run, IRS_START_TIME, IRS_input, IRS_IMMIGRATION, experimental_design){
  # Regime
  parameter_space <- experimental_design$PS[design_ID]
  scenario <- experimental_design$scenario[design_ID]
  experiment <- experimental_design$exp[design_ID] # Use 000 for checkpoints and 001 for control
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

set_MDA <- function(design_ID, run, experimental_design){
  # Regime
  parameter_space <- experimental_design$PS[design_ID]
  scenario <- experimental_design$scenario[design_ID]
  experiment <- experimental_design$exp[design_ID] # Use 000 for checkpoints and 001 for control
  base_name <- paste('PS',parameter_space,'_',scenario,'_E',experiment,'_R',run,sep='')
  output_file <- paste(base_name,'.py',sep = '')
  param_data <- get_parameter_reference(output_file)
  
  # Turn on MDA
  param_data[param_data$param=='MDA_ON',] <- set_parameter(param_data, 'MDA_ON', 'True')
  
  # Set parameters
  param_data[param_data$param=='MDA_START_TIMES',] <- set_parameter(param_data, 'MDA_START_TIMES', paste('[',experimental_design$MDA_START[design_ID],']',sep=''))
  param_data[param_data$param=='HOST_FAIL_RATE',] <- set_parameter(param_data, 'HOST_FAIL_RATE', paste('[',experimental_design$HOST_FAIL_RATE[design_ID],']',sep=''))
  param_data[param_data$param=='DRUG_EFF_DURATION',] <- set_parameter(param_data, 'DRUG_EFF_DURATION', paste('[',experimental_design$DRUG_EFF_DURATION[design_ID],']',sep=''))
  param_data[param_data$param=='MDA_IMMIGRATION_RATE_FACTORS',] <- set_parameter(param_data, 'MDA_IMMIGRATION_RATE_FACTORS', paste('[',experimental_design$MDA_IMMIGRATION_RATE_FACTORS[design_ID],']',sep=''))
  
  # Write parameter file
  param_data$output <- paste(param_data$param,param_data$value,sep='=')
  write_lines(param_data$output, output_file)
}


# Functions to generate files ---------------------------------------------

get_random_seed <- function(PS, scenario, run_range){
  seeds <- c()
  for (run in run_range){
    x <- readLines(paste('PS',PS,'_',scenario,'_','E000_R',run,'.py',sep=''))
    seeds <- c(seeds, parse_number(x[1]))
  }
  return(seeds)
}

# Function to create the necesary files and pipeline for a single run of an experiment.
# Each run has its own random seed across experiments.
create_run <- function(design_ID, run, RANDOM_SEED, experimental_design){
  
  collapse_with_commas <- function(x){
    paste('[',paste(x,collapse=','),']',sep='')
  }
  
  
  # Regime
  parameter_space <- experimental_design$PS[design_ID]
  scenario <- experimental_design$scenario[design_ID]
  experiment <- experimental_design$exp[design_ID] # Use 000 for checkpoints and 001 for control
  base_name <- paste('PS',parameter_space,'_',scenario,'_E',experiment,'_R',run,sep='')
  param_data <- get_parameter_reference()
  
  # General parameters
  param_data[param_data$param=='RANDOM_SEED',] <- set_parameter(param_data, 'RANDOM_SEED', RANDOM_SEED)
  T_END <- experimental_design$T_END[design_ID]
  param_data[param_data$param=='T_END',] <- set_parameter(param_data, 'T_END', T_END)
  param_data[param_data$param=='VERIFICATION_ON',] <- set_parameter(param_data, 'VERIFICATION_ON', 'False')
  param_data[param_data$param=='VERIFICATION_PERIOD',] <- set_parameter(param_data, 'VERIFICATION_PERIOD', T_END)
  
  # Scenario
  if(scenario=='N'){
    param_data[param_data$param=='SELECTION_MODE',] <- set_parameter(param_data, 'SELECTION_MODE', "\'NEUTRALITY\'")
    TRANSITION_RATE_NOT_IMMUNE <- set_transition_rate_neutral(parameter_space, experiment='001')
    param_data[param_data$param=='TRANSITION_RATE_NOT_IMMUNE',] <- set_parameter(param_data, 'TRANSITION_RATE_NOT_IMMUNE', TRANSITION_RATE_NOT_IMMUNE)
  }
  if(scenario=='G'){
    param_data[param_data$param=='SELECTION_MODE',] <- set_parameter(param_data, 'SELECTION_MODE', "\'GENERAL_IMMUNITY\'")
    x <- set_generalized_immunity(parameter_space, run = run)
    param_data[param_data$param=='GENERAL_IMMUNITY_PARAMS',] <- set_parameter(param_data, 'GENERAL_IMMUNITY_PARAMS', x[[1]])
    param_data[param_data$param=='N_INFECTIONS_FOR_GENERAL_IMMUNITY',] <- set_parameter(param_data, 'N_INFECTIONS_FOR_GENERAL_IMMUNITY', x[[2]])
    param_data[param_data$param=='CLEARANCE_RATE_IMMUNE',] <- set_parameter(param_data, 'CLEARANCE_RATE_IMMUNE', x[[3]])
  }
  
  # Populations. This duplicates the following parameters as the number of populations.
  N_POPULATIONS <- experimental_design$N_POPULATIONS[design_ID]
  param_data[param_data$param=='N_POPULATIONS',] <- set_parameter(param_data, 'N_POPULATIONS', N_POPULATIONS)
  # This creates the distance matrix, assuming equal distances between all populations (parameter populations_dist)
  DISTANCE_MAT <- matrix(experimental_design$populations_dist[design_ID], ncol=N_POPULATIONS, nrow=N_POPULATIONS)
  diag(DISTANCE_MAT) <- 1
  DISTANCE_MAT_vectorized <- c()
  for (i in 1:N_POPULATIONS){
    DISTANCE_MAT_vectorized <- paste(DISTANCE_MAT_vectorized, collapse_with_commas(DISTANCE_MAT[i,]), sep=',')
  }
  DISTANCE_MAT <- paste(str_replace(DISTANCE_MAT_vectorized,',','['),']',sep='')
  param_data[param_data$param=='DISTANCE_MAT',] <- set_parameter(param_data, 'DISTANCE_MAT', DISTANCE_MAT)
  # Here are some parameters with fixed values
  param_data[param_data$param=='N_HOSTS',] <- set_parameter(param_data, 'N_HOSTS', collapse_with_commas(rep(10000,N_POPULATIONS)))
  param_data[param_data$param=='N_INITIAL_INFECTIONS',] <- set_parameter(param_data, 'N_INITIAL_INFECTIONS', collapse_with_commas(rep(20,N_POPULATIONS)))
  param_data[param_data$param=='BITING_RATE_RELATIVE_AMPLITUDE',] <- set_parameter(param_data, 'BITING_RATE_RELATIVE_AMPLITUDE', collapse_with_commas(rep(0,N_POPULATIONS)))
  param_data[param_data$param=='BITING_RATE_PEAK_PHASE',] <- set_parameter(param_data, 'BITING_RATE_PEAK_PHASE', collapse_with_commas(rep(0,N_POPULATIONS)))
  param_data[param_data$param=='IMMIGRATION_RATE',] <- set_parameter(param_data, 'IMMIGRATION_RATE', collapse_with_commas(rep(1,N_POPULATIONS)))
  
  # Biting rates
  BITING_RATE_MEAN <- experimental_design$BITING_RATE_MEAN[design_ID]
  param_data[param_data$param=='BITING_RATE_MEAN',] <- set_parameter(param_data, 'BITING_RATE_MEAN', paste('[',paste(rep(BITING_RATE_MEAN,N_POPULATIONS),collapse = ','),']',sep=''))
  mathematica_file <- experimental_design$DAILY_BITING_RATE_DISTRIBUTION[design_ID]
  biting_rate_mathematica <- read_csv(mathematica_file, col_names = c('day','num_mosquitos'))
  DAILY_BITING_RATE_DISTRIBUTION <- biting_rate_mathematica$num_mosquitos
  param_data[param_data$param=='DAILY_BITING_RATE_DISTRIBUTION',] <- set_parameter(param_data, 'DAILY_BITING_RATE_DISTRIBUTION', paste('[',paste(DAILY_BITING_RATE_DISTRIBUTION, collapse=','),']',sep=''))
  
  # Genetic diversity
  N_GENES_INITIAL <- experimental_design$N_GENES_INITIAL[design_ID]
  param_data[param_data$param=='N_GENES_INITIAL',] <- set_parameter(param_data, 'N_GENES_INITIAL', N_GENES_INITIAL)
  N_ALLELES_INITIAL <- experimental_design$N_ALLELES_INITIAL[design_ID]
  param_data[param_data$param=='N_ALLELES_INITIAL',] <- set_parameter(param_data, 'N_ALLELES_INITIAL', paste('N_LOCI*[',N_ALLELES_INITIAL,']',sep=''))
  
  # Checkpoints
  T_BURNIN <- experimental_design$T_BURNIN[design_ID]
  if (experiment == '000'){
    param_data[param_data$param=='SAMPLE_DB_FILENAME',] <- set_parameter(param_data, 'SAMPLE_DB_FILENAME', paste('\'\"',base_name,'.sqlite\"\'',sep='')) # The run ID will be added while running the job (in the sbatch execution).
    param_data[param_data$param=='SAVE_TO_CHECKPOINT',] <- set_parameter(param_data, 'SAVE_TO_CHECKPOINT', 'True')
    param_data[param_data$param=='CHECKPOINT_SAVE_PERIOD',] <- set_parameter(param_data, 'CHECKPOINT_SAVE_PERIOD', T_END) # The save period should be the T_END
    param_data[param_data$param=='CHECKPOINT_SAVE_FILENAME',] <- set_parameter(param_data, 'CHECKPOINT_SAVE_FILENAME', paste('\'\"',base_name,'_CP.sqlite\"\'',sep=''))
    param_data[param_data$param=='LOAD_FROM_CHECKPOINT',] <- set_parameter(param_data, 'LOAD_FROM_CHECKPOINT', 'False')
    param_data[param_data$param=='CHECKPOINT_LOAD_FILENAME',] <- set_parameter(param_data, 'CHECKPOINT_LOAD_FILENAME', '\'\"\"\'')
    param_data[param_data$param=='T_BURNIN',] <- set_parameter(param_data, 'T_BURNIN', T_BURNIN)
  }
  if (experiment != '000'){ # This section prepares a parameter file to load a checkpoint and run an experiment
    param_data[param_data$param=='SAMPLE_DB_FILENAME',] <- set_parameter(param_data, 'SAMPLE_DB_FILENAME', paste('\'\"',base_name,'.sqlite\"\'',sep='')) # The run ID will be added while running the job (in the sbatch execution).
    param_data[param_data$param=='SAVE_TO_CHECKPOINT',] <- set_parameter(param_data, 'SAVE_TO_CHECKPOINT', 'False')
    param_data[param_data$param=='CHECKPOINT_SAVE_PERIOD',] <- set_parameter(param_data, 'CHECKPOINT_SAVE_PERIOD', 0)
    param_data[param_data$param=='CHECKPOINT_SAVE_FILENAME',] <- set_parameter(param_data, 'CHECKPOINT_SAVE_FILENAME', '\'\"\"\'')
    param_data[param_data$param=='LOAD_FROM_CHECKPOINT',] <- set_parameter(param_data, 'LOAD_FROM_CHECKPOINT', 'True')
    param_data[param_data$param=='CHECKPOINT_LOAD_FILENAME',] <- set_parameter(param_data, 'CHECKPOINT_LOAD_FILENAME', paste('\'\"PS',parameter_space,'_',scenario,'_E000','_R',run,'_CP.sqlite\"\'',sep=''))
    param_data[param_data$param=='T_BURNIN',] <- set_parameter(param_data, 'T_BURNIN', T_BURNIN) # The burnin value should be the value where the checkpoint was taken
  }
  
  # Write parameter file
  output_file=paste(base_name,'.py',sep = '')
  param_data$output <- paste(param_data$param,param_data$value,sep='=')
  write_lines(param_data$output, output_file)
}

# This function generates parameter files and corresponding job files for
# several experiments and runs. It keeps the random seed for each RUN the same.
# Also possible to provide a random seed.
# row_range is the row numbers in the design data frame
generate_files <- function(row_range, run_range, random_seed=NULL, experimental_design){
  # Generate a parameter file for each combination of experiment and run
  for (RUN in run_range){
    # Keep random seed the same across runs
    RANDOM_SEED <- ifelse(is.null(random_seed),round(runif(n = 1, min=1, max=10^7),0),random_seed[RUN])
    for (design_ID in row_range){
      # Create parameter file with main parameters
      create_run(design_ID, RUN, RANDOM_SEED, experimental_design)
      # Set IRS
      if (!is.na(experimental_design$IRS_START_TIMES[design_ID])){
        IRS_scheme <- data.frame(IRS_START_TIME=str_split(experimental_design$IRS_START_TIMES[design_ID], ',')[[1]],
                                 IRS_input=str_split(experimental_design$IRS_input[design_ID], ',')[[1]],
                                 IRS_IMMIGRATION=str_split(experimental_design$IRS_IMMIGRATION[design_ID], ',')[[1]],
                                 stringsAsFactors = F)
        for (i in 1:nrow(IRS_scheme)){
          set_IRS(design_ID, RUN, IRS_scheme$IRS_START_TIME[i], IRS_scheme$IRS_input[i], IRS_scheme$IRS_IMMIGRATION[i], experimental_design)
        }
      }
      # Set MDA
      if (!is.na(experimental_design$MDA_START[design_ID])){
        set_MDA(design_ID, RUN, experimental_design)
      }
    }
  }
  
  # Generate batch file for each experiment
  SLURM_ARRAY_RANGE <- paste("\'",min(run_range),'-',max(run_range),"\'",sep='')
  for (design_ID in row_range){
    parameter_space <- experimental_design$PS[design_ID]
    scenario <- experimental_design$scenario[design_ID]
    experiment <- experimental_design$exp[design_ID] # Use 000 for checkpoints and 001 for control
    base_name <- paste('PS',parameter_space,scenario,'E',experiment,sep='')
    # Write the job file for the exepriment
    job_lines <- readLines('~/Documents/malaria_interventions/job_file_ref.sbatch')
    wall_time <-  experimental_design$wall_time[design_ID]
    mem_per_cpu <-  experimental_design$mem_per_cpu[design_ID]
    job_lines[2] <- paste('#SBATCH --job-name=',paste(parameter_space,scenario,'E',experiment,sep=''),sep='')
    job_lines[3] <- paste('#SBATCH --time=',wall_time,sep='')
    job_lines[4] <- paste('#SBATCH --output=slurm_output/',base_name,'_%A_%a.out',sep='')
    job_lines[5] <- paste('#SBATCH --error=slurm_output/',base_name,'_%A_%a.err',sep='')
    job_lines[6] <- paste('#SBATCH --array=',SLURM_ARRAY_RANGE,sep='')
    job_lines[9] <- paste('#SBATCH --mem-per-cpu=',mem_per_cpu,sep='')
    job_lines[19] <- paste("PS='",parameter_space,"'",sep='')
    job_lines[20] <- paste("scenario='",scenario,"'",sep='')
    job_lines[21] <- paste("exp='",experiment,"'",sep='')
    if (experiment == '000'){
      job_lines[26] <- "CHECKPOINT='create'"
    }
    if (experiment != '000'){
      job_lines[26] <- "CHECKPOINT='load'"
    }
    output_file=paste(base_name,'.sbatch',sep = '')
    write_lines(job_lines, output_file)
  }
  
  # Printout of generated file names
  cat('Generated files:\n')
  for (e in row_range){
    cat(paste('PS',experimental_design$PS[e],experimental_design$scenario[e],'E',experimental_design$exp[e],'.sbatch','\n',sep=''))
    for (r in run_range){
      cat(paste('--- PS',experimental_design$PS[e],'_',experimental_design$scenario[e],'_E',experimental_design$exp[e],'_R',r,'\n',sep=''))
    }
  }
  
  cat('Run this on Midway: \n')
  for (e in row_range){
    cat(paste('sbatch PS',experimental_design$PS[e],experimental_design$scenario[e],'E',experimental_design$exp[e],'.sbatch','\n',sep=''))
  }
}


create_intervention_scheme_IRS <- function(PS_benchmark, scenario_benchmark, IRS_START_TIMES=NULL, immigration_range, length_range, coverage_range, write_to_file=T, design_ref=NULL){
  # This fixes the 1. in the file name produced in Mathematica
  files <- list.files(path = '~/Documents/malaria_interventions_data/', full.names = T, pattern = '1\\._')
  if(length(files)>0){sapply(files, function(x){file.rename(from = x, to = str_replace(x, '\\.', ''))})}
  
  # PS_benchmark is the parameter space for which intervention is tested. it should already have the checkpoint (000) and control (001) experiments in the google sheet.
  # IRS_START_TIMES can also be something like: '29000,35000'
  if(is.null(design_ref)){
    design_ref <- loadExperiments_GoogleSheets() # Get data design
  }
  design_ref <- subset(design_ref, PS==PS_benchmark & scenario==scenario_benchmark)
  reference_row <- which(grepl(PS_benchmark, design_ref$PS))[2] # reference_row is the row in the design data set which is the reference for this set of IRS exepriments. Should be the line of the control experiment (not the checkpoint)
  
  scheme <- expand.grid(PS=design_ref$PS[reference_row],
                        scenario=design_ref$scenario[reference_row],
                        BITING_RATE_MEAN=design_ref$BITING_RATE_MEAN[reference_row],
                        DAILY_BITING_RATE_DISTRIBUTION=design_ref$DAILY_BITING_RATE_DISTRIBUTION[reference_row],
                        N_GENES_INITIAL=design_ref$N_GENES_INITIAL[reference_row],
                        N_LOCI=design_ref$N_LOCI[reference_row],
                        N_ALLELES_INITIAL=design_ref$N_ALLELES_INITIAL[reference_row],
                        N_POPULATIONS=design_ref$N_POPULATIONS[reference_row],
                        populations_dist=design_ref$populations_dist[reference_row],
                        T_BURNIN=design_ref$T_BURNIN[reference_row],
                        T_END=design_ref$T_END[reference_row],
                        IRS_START_TIMES=IRS_START_TIMES,
                        IRS_IMMIGRATION=immigration_range,
                        IRS_length=length_range,
                        IRS_coverage=coverage_range,
                        MDA_START=NA, # This so to not have an MDA (for function generate_files)
                        wall_time=design_ref$wall_time[reference_row],
                        mem_per_cpu=design_ref$mem_per_cpu[reference_row],
                        stringsAsFactors = F)
  scheme$exp <- sprintf('%0.3d', 2:(nrow(scheme)+1)) # Start with 002 because 000 and 001 are checkpoint and control
  # Important!!! The files for the IRS in the next line are generated separately in Mathematica, and should already exist.
  scheme$IRS_input <- paste('IRS_120_300_',scheme$IRS_coverage,'_',scheme$IRS_length,'.csv',sep='')
  if (write_to_file){
    write_csv(scheme, paste('PS',scheme$PS[1],'_',scheme$scenario[1],'_IRS_scheme.csv',sep=''))
  }
  return(scheme)
}


# Create parameter and job files -------------------------------------------

setwd('~/Documents/malaria_interventions_data/')

# Clear previous files if necessary
clear_previous_files(parameter_space='12', scenario = 'S', exclude_sqlite = F, exclude_CP = F, exclute_control = F)
# Create the reference experiments (checkpoint and control)
design <- loadExperiments_GoogleSheets() # Get data design 
generate_files(row_range = 23:24, run_range = 1:2, experimental_design = design, random_seed = c(9198087,3012346))

# Create the corresponding IRS experiments  
PS <- sprintf('%0.2d', 12)
for (ps in PS){
  design_irs <- create_intervention_scheme_IRS(PS_benchmark = ps, scenario_benchmark = 'S', IRS_START_TIMES = '29160', immigration_range=c(0), length_range=c(720,1800,3600), coverage_range=0.9, write_to_file = T, design_ref=design)
  generate_files(row_range = 1:nrow(design_irs), run_range = 1:2, random_seed = get_random_seed(ps, 'S', 1:2), design_irs)
}

for (ps in PS){
  print(get_random_seed(ps, 'S', 1:2))
}


# Generate command to run experiment jobs
paste('for i in ', paste(sprintf('%0.3d', 101:221), collapse=' '),'; do sbatch PS04SE$i.sbatch; done', sep='')

# First run the checkpoints
for i in 06 07 08 09 10 11; do sbatch 'PS'$i'SE000.sbatch'; done;

# Then run control and interventions
exp <- sprintf('%0.3d', 1:4)
PS <- sprintf('%0.2d', 12)
jobids <- c('47323910') # 1 job id per PS
for (ps in PS){
  for (e in exp){
    cat(paste('sbatch -d afterok:',jobids[which(PS==ps)],' PS',ps,'SE',e,'.sbatch',sep=''));cat('\n')
  }
}
    
for x in 001 002 003 004; do for i in 06 07 08; do sbatch 'PS'$i'SE'$x'.sbatch'; done; done


design_irs <- create_intervention_scheme_IRS(PS_benchmark = '09', scenario_benchmark = 'S', IRS_START_TIMES = '29160', immigration_range=c(0,0.001), length_range=c(1800,3600), coverage_range=0.9, write_to_file = T, design_ref=design)



# plots -------------------------------------------------------------------
setwd('~/Documents/malaria_interventions/')
seasonality <- read_csv('mosquito_population_seasonality.csv', col_names = c('day','num_mosquitos'))
seasonality$grp <- 'S'
seasonality2 <- NULL
for (i in 1:10){
  seasonality2 <- rbind(seasonality2,seasonality)
}
seasonality=seasonality2
seasonality$day <- 1:3600
irs01 <- read_csv('IRS01.csv', col_names = c('day','num_mosquitos'))
irs01$grp <- 'IRS01'
x <- rbind(seasonality,irs01)
x %>% ggplot(aes(day, num_mosquitos, color=grp))+geom_line()

