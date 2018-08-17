library(tidyverse)
library(magrittr)
library(sqldf)
setwd('~/Documents/malaria_interventions')


# General functions -------------------------------------------------------

collapse_with_commas <- function(x, use_brackets=T){
  if (use_brackets){
    paste('[',paste(x,collapse=','),']',sep='')
  } else {
    paste(x,collapse=',')
  }
}

sink.reset <- function(){
  for(i in seq_len(sink.number())){
    sink(NULL)
  }
}

# Function to set a parameter
set_parameter <- function(param_data, parameter, value){
  param_data$value[param_data$param==parameter] <- value
  return(subset(param_data, param==parameter))
}

# Function to get reference parameter file
get_parameter_reference <- function(parameter_file_ref='~/Documents/malaria_interventions/parameter_file_ref.py'){
  reference <- readLines(parameter_file_ref)
  for (l in 1:length(reference)){
      if (str_detect(reference[l],'#')){
        reference[l] <- str_sub(reference[l],1,str_locate(reference[l],'#')[1]-1)
      }
  }
  
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
clear_previous_files <- function(parameter_space=NULL, scenario=NULL, experiment=NULL, run = NULL, exclude_sqlite=T, exclude_CP=T, exclude_control=T, test=F){
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
  if (!is.null(run_range)){
    files <- files[str_detect(files,paste('_R',run,sep=''))]
  }
  if(exclude_CP){
    files <- files[!str_detect(files,'E000')]
  }
  if(exclude_control){
    files <- files[!str_detect(files,'E001')]
  }
  if(exclude_sqlite){
    files <- files[!str_detect(files,'\\.sqlite')]
  }
  if (test){
    print('test mode, not actually removing')
    print(files)
  } else {
    print(files)
    file.remove(files)
  }
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
#cannot be aggregated across runs. It makes most sense to fir the curve to an
#experiment without interevntions and in stable state, so 001 is by default.
set_generalized_immunity <- function(parameter_space, experiment='001', run){
  require('rPython')
  sqlite_file <- paste('/home/shai/Documents/malaria_interventions_data/','PS',parameter_space,'_S_E',experiment,'_R',run,'.sqlite',sep='')
  pyFile <- readLines('/home/shai/Documents/malaria_interventions/generalized_immunity_fitting.py')
  pyFile[11] <- paste('path=','"',sqlite_file,'"',sep='')
  writeLines(pyFile,'generalized_immunity_fitting_experiment.py')
  
  res <- try(python.load('generalized_immunity_fitting_experiment.py'))
  if(inherits(res, "try-error"))
  {
    return(NA)
  } else {
    GENERAL_IMMUNITY_PARAMS <- python.get('generalImmunityParams')
    GENERAL_IMMUNITY_PARAMS <- paste('[',paste(GENERAL_IMMUNITY_PARAMS, collapse = ','),']',sep='')
    N_INFECTIONS_FOR_GENERAL_IMMUNITY <- python.get('infectionTimesToImmune')
    CLEARANCE_RATE_IMMUNE <- python.get('clearanceRateConstantImmune')
    file.remove('generalized_immunity_fitting_experiment.py')
    return(list(GENERAL_IMMUNITY_PARAMS,N_INFECTIONS_FOR_GENERAL_IMMUNITY,CLEARANCE_RATE_IMMUNE))
  }
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

get_random_seed <- function(PS, scenario, experiment='000', run_range){
  seeds <- c()
  for (run in run_range){
    x <- try(readLines(paste('PS',PS,'_',scenario,'_','E',experiment,'_R',run,'.py',sep='')))
    if (inherits(x, "try-error")){
      print('Cannot find parameter file to extract random number. Most probably a generalized-immunity experiment/run that does not exist.')
      seeds <- c(seeds, NA)
    } else {
      seeds <- c(seeds, parse_number(x[1]))  
    }
  }
  return(seeds)
}

# Function to create the necesary files and pipeline for a single run of an experiment.
# Each run has its own random seed across experiments.
create_run <- function(design_ID, run, RANDOM_SEED, experimental_design, biting_rate_mathematica=NULL){
  
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
    try_fit <- set_generalized_immunity(parameter_space, run = run)
    if (is.na(try_fit)[1]){
      print(paste('ERROR (create_run): Could not fit function for generalized immunity, skipping and NOT PRODUCING FILE ',base_name,'.py',sep = ''))
      return(NULL)
    }
    param_data[param_data$param=='GENERAL_IMMUNITY_PARAMS',] <- set_parameter(param_data, 'GENERAL_IMMUNITY_PARAMS', try_fit[[1]])
    param_data[param_data$param=='N_INFECTIONS_FOR_GENERAL_IMMUNITY',] <- set_parameter(param_data, 'N_INFECTIONS_FOR_GENERAL_IMMUNITY', try_fit[[2]])
    param_data[param_data$param=='CLEARANCE_RATE_IMMUNE',] <- set_parameter(param_data, 'CLEARANCE_RATE_IMMUNE', try_fit[[3]])
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
  # It is possible to provide a specific file for biting rates for each
  # experiments, but usually that is not necessary so no need to re-read this
  # when producing batches of files. Just provide it in the call to the
  # function.
  if (is.null(biting_rate_mathematica)){
    mathematica_file <- experimental_design$DAILY_BITING_RATE_DISTRIBUTION[design_ID]
    biting_rate_mathematica <- read_csv(mathematica_file, col_names = c('day','num_mosquitos'))
  }
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



# A function to generate sbatch files on demand
generate_sbatch <- function(ps, scen, experiment, runs, unzip_py_files){
  base_name <- paste('PS',ps,scen,'E',experiment,sep='')
  SLURM_ARRAY_RANGE <- collapse_with_commas(runs, use_brackets = F)
  
  # This is to get the correct walltime and memory for Midway
  if (experiment == '000'){
    experimental_design <- subset(design, PS==ps & scenario==scen & exp=='000')
  }
  if (experiment != '000'){
    experimental_design <- subset(design, PS==ps & scenario==scen & exp=='001')
  }  
  
  # Write the job file for the exepriment
  job_lines <- readLines('~/Documents/malaria_interventions/job_file_ref.sbatch')
  wall_time <-  experimental_design$wall_time
  mem_per_cpu <-  experimental_design$mem_per_cpu
  job_lines[2] <- paste('#SBATCH --job-name=',paste(ps,scen,'E',experiment,sep=''),sep='')
  job_lines[3] <- paste('#SBATCH --time=',wall_time,sep='')
  job_lines[4] <- paste('#SBATCH --output=slurm_output/',base_name,'_%A_%a.out',sep='')
  job_lines[5] <- paste('#SBATCH --error=slurm_output/',base_name,'_%A_%a.err',sep='')
  job_lines[6] <- paste('#SBATCH --array=',SLURM_ARRAY_RANGE,sep='')
  job_lines[9] <- paste('#SBATCH --mem-per-cpu=',mem_per_cpu,sep='')
  job_lines[19] <- paste("PS='",ps,"'",sep='')
  job_lines[20] <- paste("scenario='",scen,"'",sep='')
  job_lines[21] <- paste("exp='",experiment,"'",sep='')
  if (experiment == '000'){
    job_lines[26] <- "CHECKPOINT='create'"
  }
  if (experiment != '000'){
    job_lines[26] <- "CHECKPOINT='load'"
  }  
  
  output_file=paste(base_name,'.sbatch',sep = '')
  write_lines(job_lines, output_file)
  cat(paste('sbatch ',base_name,'.sbatch',sep=''));cat('\n')
  
  if (unzip_py_files){
    unzip(paste(ps,'_',scen,'_py.zip',sep=''), files = paste('PS',ps,'_',scen,'_','E',experiment,'_R',runs,'.py',sep=''))
  }
}


# This function generates parameter files and corresponding job files (to run on
# Midway), for several experiments and runs. It keeps the random seed for each
# RUN across experiments AND PARAMETER SPACES the same. Also possible to provide
# a random seed. row_range is the row numbers in the design data frame.
generate_files <- function(row_range, run_range, random_seed=NULL, experimental_design, biting_rate_file='mosquito_population_seasonality.csv'){
  
  if (!is.null(biting_rate_file)){
    biting_rate_mathematica <- read_csv(biting_rate_file, col_names = c('day','num_mosquitos'))
  }
  
  # Using this data structure makes it easier to control specific combinations
  # of runs in experimets, for eaxmple for parameter files that cannot be
  # produced because of ill fitting of curves in the generlized immunity.
  cases <- expand.grid(run = run_range, design_ID=row_range, file_created=T)
  
  # Random seeds are equal for each run across experiments AND PARAMETER SPACES
  # (e.g., the same seed for run 1 in experiments '000', '001', '002' in
  # parameter spaces 01 and 02).
  if (is.null(random_seed)){
    cases$seed <- rep(round(runif(n = length(run_range), min=1, max=10^7),0),length(row_range))
  } else {
    cases$seed <- rep(random_seed,length(row_range))
  }
  
  cases$param_file <- NA
  
  for (idx in 1:nrow(cases)){
    RANDOM_SEED <- cases$seed[idx]
    RUN <- cases$run[idx]
    design_ID <- cases$design_ID[idx]
    print(paste('Run: ',RUN,' | Row: ',design_ID,sep=''))
    
    # Create parameter file with main parameters. If the scenario is generalized
    # immunity then the fit for some runs may have not convereged (functions
    # 'create_run' and 'set_generalized_immunity'). Parameter files cannot be
    # produced in these cases and so the function skips to the next case.
    try_run <- create_run(design_ID, RUN, RANDOM_SEED, experimental_design, biting_rate_mathematica)
    if (is.null(try_run)){
      cases$file_created[idx] <- F
      next
    }
    
    # Set IRS
    if (!is.na(experimental_design$IRS_START_TIMES[design_ID])){
      IRS_scheme <- data.frame(IRS_START_TIME=str_split(experimental_design$IRS_START_TIMES[design_ID], ',')[[1]],
                               IRS_input=str_split(experimental_design$IRS_input[design_ID], ',')[[1]],
                               IRS_IMMIGRATION=str_split(experimental_design$IRS_IMMIGRATION[design_ID], ',')[[1]],
                               stringsAsFactors = F)
      for (i in 1:nrow(IRS_scheme)){
        set_IRS(design_ID = design_ID, run = RUN, IRS_START_TIME = IRS_scheme$IRS_START_TIME[i], IRS_IMMIGRATION = IRS_scheme$IRS_IMMIGRATION[i], IRS_input = IRS_scheme$IRS_input[i], experimental_design)
      }
    }
    # Set MDA
    if (!is.na(experimental_design$MDA_START[design_ID])){
      set_MDA(design_ID, RUN, experimental_design)
    }
    
    cases$param_file[idx] <- paste('PS',experimental_design$PS[design_ID],'_',experimental_design$scenario[design_ID],'_E',experimental_design$exp[design_ID],'_R',RUN,sep='')
  }
  
  # Generate batch file for each experiment
  for (design_ID in row_range){
    parameter_space <- experimental_design$PS[design_ID]
    scenario <- experimental_design$scenario[design_ID]
    experiment <- experimental_design$exp[design_ID] # Use 000 for checkpoints and 001 for control
    base_name <- paste('PS',parameter_space,scenario,'E',experiment,sep='')
    
    SLURM_ARRAY_RANGE <- cases$file_created[which(cases$design_ID==design_ID)]
    if (any(SLURM_ARRAY_RANGE==F)){  # If some runs are missing (e.g., runs for which the fit has not converged in generalized immunity)
      SLURM_ARRAY_RANGE <-  paste("\'",collapse_with_commas(run_range[SLURM_ARRAY_RANGE],F),"\'",sep='')
    } else {
      SLURM_ARRAY_RANGE <- paste("\'",min(run_range),'-',max(run_range),"\'",sep='')
    }
    
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
  
  print(cases)
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
clear_previous_files(run = 6, exclude_sqlite = F, exclude_CP = F, exclude_control = F, test = T)
# Get data design 
design <- loadExperiments_GoogleSheets() 

# Create the reference experiments (checkpoint and control)
ps_range <- sprintf('%0.2d', 27:39)
exp_range <- sprintf('%0.3d', 0:4)
run_range <- 11:50
work_scenario <- 'S'
design_subset <- subset(design, PS %in% ps_range & scenario==work_scenario)
generate_files(row_range = 1:nrow(design_subset), run_range = run_range, experimental_design = design_subset)

# # If checkpoints already exist, can create the reference experiments (control only)
for (ps in ps_range){
  print(ps)
  design_control <- subset(design, PS==ps & scenario == work_scenario & exp=='001')
  # clear_previous_files(parameter_space = ps, scenario = work_scenario, exclude_sqlite = F, exclude_CP = T, exclude_control = F, test = F)
  seeds <- get_random_seed(PS = ps, scenario = work_scenario, run_range = run_range)
  generate_files(row_range = 1, run_range = run_range, experimental_design = design_control, random_seed = seeds)
}

# Create the corresponding IRS experiments  
for (ps in ps_range){
  design_irs <- create_intervention_scheme_IRS(PS_benchmark = ps, scenario_benchmark = work_scenario, IRS_START_TIMES = '29160', immigration_range=c(0), length_range=c(720,1800,3600), coverage_range=0.9, write_to_file = F, design_ref=design)
  generate_files(row_range = 1:nrow(design_irs), run_range = run_range, random_seed = get_random_seed(ps, work_scenario, run_range = run_range), design_irs)
}

# ZIP all the PY and sbatch files
unlink('files_to_run.txt')
sink('files_to_run.txt', append = T)
# Add sbatch files
files <- list.files(path = '~/Documents/malaria_interventions_data/', pattern = 'sbatch', full.names = F) 
files <- files[str_detect(string = files, paste(work_scenario,'E',sep=''))]
files <- files[str_sub(files,3,4) %in% ps_range]
files <- files[str_sub(files,7,9) %in% exp_range]
for (i in 1:length(files)){
  cat(files[i]);cat('\n')
}
# Add py files
files <- list.files(path = '~/Documents/malaria_interventions_data/', pattern = '.py', full.names = F) 
files <- files[str_detect(string = files, paste('_',work_scenario,sep=''))]
files <- files[str_sub(files,3,4) %in% ps_range]
files <- files[str_sub(files,9,11) %in% exp_range]
files <- files[parse_number(str_sub(files,14,15)) %in% run_range]
for (i in 1:length(files)){
  cat(files[i]);cat('\n')
}
sink.reset()
# Create zip
unlink('files_to_run.zip')
system('zip files_to_run.zip -@ < files_to_run.txt')

# Copy the file to Midway and unzip it

# First run the checkpoints
paste("for i in ",paste(ps_range,collapse=' '),"; do sbatch 'PS'$i'",work_scenario,"E000.sbatch'; done;",sep='')

# Then run control and interventions
## Run in Midway terminal:
cat("rm job_ids.txt; sacct -u pilosofs --format=jobid,jobname --starttime 2018-08-07T14:48:00 --name=");cat(paste(paste(ps_range,work_scenario,'E000',sep=''),collapse = ','));cat(" >> 'job_ids.txt'")
## Copy file from Midway and run:
jobids <- read.table('job_ids.txt', header = F, skip=2) 
jobids <- na.omit(unique(parse_number(jobids$V1))) # 1 job id per PS
length(jobids)==length(ps_range)
## Copy the output of the following loop and paste in Midway
for (ps in ps_range){
  for (e in exp_range){
    cat(paste('sbatch -d afterok:',jobids[which(ps_range==ps)],' PS',ps,work_scenario,'E',e,'.sbatch',sep=''));cat('\n')
  }
}

# Or, if checkpoints are already finished:
unlink('jobs_to_run.sh')
sink('jobs_to_run.sh')
for (ps in ps_range){
  for (e in exp_range[-1]){
    cat(paste('sbatch PS',ps,work_scenario,'E',e,'.sbatch',sep=''));cat('\n')
  }
}
sink.reset()

# Zip files to reduce the clutter -----------------------------------------

system('rm *.sbatch')

scenario_range <- c('G')
exp_range <- sprintf('%0.3d', 0:4)
ps_range <- sprintf('%0.2d', 27:39)

for (ps in ps_range){
  for (scenario in scenario_range){
    unlink('files_tmp.txt')
    sink('files_tmp.txt', append = T)
    # Add py files
    files <- list.files(path = '~/Documents/malaria_interventions_data/', pattern = '.py', full.names = F) 
    files <- files[str_sub(files,3,4) %in% ps]
    files <- files[str_sub(files,6,6) %in% scenario]
    files <- files[str_sub(files,9,11) %in% exp_range]
    for (i in 1:length(files)){
      cat(files[i]);cat('\n')
    }
    sink.reset()
    # Create zip
    # unlink(paste(ps,'_sbatch_py.zip',sep=''))
    if (file.size('files_tmp.txt')>10){
      system(paste('zip ',ps,'_',scenario,'_py.zip -@ < files_tmp.txt -mu',sep=''))
    }
  }
}


# Verify files ------------------------------------------------------------

# --- sqlite files ---
files <- list.files(path = '~/Documents/malaria_interventions_data/', pattern = '\\.sqlite', full.names = F)
# Also consider adding following information: extracting random seeds; existence of corresponding parameter files
files_sqlite <- tibble(file_sqlite=files,
                    PS = sapply(str_split(files,'_'),function (x) parse_number(x[1])),
                    scenario=sapply(str_split(files,'_'),function (x) x[2]),
                    experiment=sapply(str_split(files,'_'),function (x) str_sub(x[3],2,4)),
                    run= sapply(str_split(files,'_'),function (x) parse_number(x[4])),
                    CP=sapply(str_split(files,'_'),function (x) str_detect(x[5],'CP')),
                    size=round(file.size(files)/1024^2,2)
)
files_sqlite$PS <- sprintf('%0.2d', files_sqlite$PS)

files_sqlite$run_time <- unlist(lapply(files_sqlite$file_sqlite[is.na(files_sqlite$CP)], function (f) {
  print(f)
  db <- dbConnect(SQLite(), dbname = paste('~/Documents/malaria_interventions_data/',f,sep=''))
  x <- dbGetQuery(db, 'SELECT max(time) FROM summary')
  dbDisconnect(db)
  x
}))
unique(files_sqlite$run_time)
files_sqlite %>% filter(is.na(run_time))


files_sqlite %<>% mutate(scenario=factor(scenario, levels=c('S','N','G')))  %>% arrange(CP, scenario, PS, experiment, run)         
files_sqlite %>% filter(is.na(CP) & as.numeric(PS)>=27) %>% group_by(PS,scenario,experiment) %>% 
  summarise(runs_completed=length(run)) %>% print(n = Inf)
files_sqlite %>% group_by(scenario, PS, experiment) %>% summarise(x=sum(runs_completed)) %>% print(n = Inf)
files_sqlite %>% filter(scenario=='G') %>% print(n = Inf)

#--- PY files ---
files <- list.files(path = '~/Documents/malaria_interventions_data/', pattern = 'py.zip', full.names = T)
files <- map(files, function(f){
  unzip(f, list=T)
}) %>% bind_rows()

files_py <- tibble(file_py=files$Name,
                   PS = sapply(str_split(files$Name,'_'),function (x) parse_number(x[1])),
                   scenario=sapply(str_split(files$Name,'_'),function (x) x[2]),
                   experiment=sapply(str_split(files$Name,'_'),function (x) str_sub(x[3],2,4)),
                   run= sapply(str_split(files$Name,'_'),function (x) parse_number(x[4])),
                   size=round(files$Length,2),
                   date=files$Date
)
files_py$PS <- sprintf('%0.2d', files_py$PS)

files_py %<>% filter(as.numeric(PS)>=27) %>% mutate(scenario=factor(scenario, levels=c('S','N','G'))) %>% arrange(scenario, PS, experiment, run)
files_py %<>% filter(as.numeric(PS)>=27) %>% group_by(PS,scenario,experiment) %>% summarise(runs_count=length(run))
files_py %>% filter(scenario=='G') %>% filter(as.numeric(PS)>27) %>% distinct(PS,run)

# --- CP sqlite files ---
# In Midway run: ls *.sqlite >> CP_files.txt and copy the file to local computer
files <- read.table('CP_files.txt',header = F, stringsAsFactors = F)$V1
files_CP <- tibble(file_CP=files,
          PS = sapply(str_split(files,'_'),function (x) parse_number(x[1])),
          scenario=sapply(str_split(files,'_'),function (x) x[2]),
          experiment=sapply(str_split(files,'_'),function (x) str_sub(x[3],2,4)),
          run= sapply(str_split(files,'_'),function (x) parse_number(x[4])),
          CP=sapply(str_split(files,'_'),function (x) str_detect(x[5],'CP'))
)
files_CP$PS <- sprintf('%0.2d', files_CP$PS)

missing_py_files <- full_join(files_CP, files_py, by = c("PS", "scenario", "experiment", "run")) %>% 
  filter(as.numeric(PS)>=27) %>% 
  select(scenario, PS, experiment, run, file_CP, file_py) %>% 
  mutate(scenario=factor(scenario, levels=c('S','N','G'))) %>% 
  arrange(PS, scenario, experiment, run)

missing_py_files %>% filter(is.na(file_py))

#--- Check seeds ---
for(i in 1:nrow(files_df)){
  print(i)
  unzip(paste(files_df$PS[i],'_',files_df$scenario[i],'_py.zip',sep=''), files = files_df$file_py[i])
  files_df$seed[i] <- get_random_seed(PS = files_df$PS[i], scenario = as.character(files_df$scenario[i]), experiment = files_df$experiment[i], run_range = files_df$run[i])
  unlink(files_df$file_py[i])
}

seed_check <- files_df %>% distinct(scenario, PS, run, seed) %>% group_by(scenario, PS) %>% summarise(length(run))


# Generate sbatch files to create the missing checkpoints -----------------

ps_scen_combinations <- as.tibble(files_CP_missing) %>% distinct(scenario,PS)

for (i in 1:nrow(ps_scen_combinations)){
  ps <- ps_scen_combinations$PS[i]
  scen <- ps_scen_combinations$scenario[i]

  unzip(paste(ps,'_',scen,'_py.zip',sep=''), files = subset(files_CP_missing, scenario==scen & PS==ps)$file_py)
  experiment <- '000' # Use 000 for checkpoints and 001 for control
  base_name <- paste('PS',ps,scen,'E',experiment,sep='')

  SLURM_ARRAY_RANGE <- collapse_with_commas(subset(files_CP_missing, scenario==scen & PS==ps)$run, use_brackets = F)

  experimental_design <- subset(design, PS==ps & scenario==scen & exp=='000')
  # Write the job file for the exepriment
  job_lines <- readLines('~/Documents/malaria_interventions/job_file_ref.sbatch')
  wall_time <-  experimental_design$wall_time
  mem_per_cpu <-  experimental_design$mem_per_cpu
  job_lines[2] <- paste('#SBATCH --job-name=',paste(ps,scen,'E',experiment,sep=''),sep='')
  job_lines[3] <- paste('#SBATCH --time=',wall_time,sep='')
  job_lines[4] <- paste('#SBATCH --output=slurm_output/',base_name,'_%A_%a.out',sep='')
  job_lines[5] <- paste('#SBATCH --error=slurm_output/',base_name,'_%A_%a.err',sep='')
  job_lines[6] <- paste('#SBATCH --array=',SLURM_ARRAY_RANGE,sep='')
  job_lines[9] <- paste('#SBATCH --mem-per-cpu=',mem_per_cpu,sep='')
  job_lines[19] <- paste("PS='",ps,"'",sep='')
  job_lines[20] <- paste("scenario='",scen,"'",sep='')
  job_lines[21] <- paste("exp='",experiment,"'",sep='')
  job_lines[26] <- "CHECKPOINT='create'"
  output_file=paste(base_name,'.sbatch',sep = '')
  write_lines(job_lines, output_file)
  cat(paste('sbatch ',base_name,'.sbatch',sep=''));cat('\n')
}


# Add a line in the parameter files ---------------------------------------
# This adds a line to the end of each py file
ps_range <- sprintf('%0.2d', 28:39)
work_scenario <- 'G'
exp_range <- sprintf('%0.3d', 2:4)
run_range <- 1:10
for (ps in ps_range){
  for (e in exp_range){
    for (r in run_range){
      file <- paste('PS',ps,'_',work_scenario,'_E',e,'_R',r,'.py',sep='')
      print(file)
      unzip(paste(ps,'_',work_scenario,'_py.zip',sep=''), file)
      if(file.exists(file)) {write('POOLSIZE_BOUNCE_BACK_AFTER_INTERVENTION=False',file, append=TRUE)}
    }
  }
}

for (ps in ps_range){
  generate_sbatch(ps, scen = 'N', experiment = '003', runs = 1:10, unzip_py_files = F)
}

# Or generate sbatch for particular runs, as in the generalized immunity
cases <- files_py %>% filter(scenario=='G') %>% filter(as.numeric(PS)>27) %>% distinct(PS,run)
for (ps in ps_range){
  generate_sbatch(ps, scen = 'G', experiment = '002', runs = subset(cases, PS==ps)$run, unzip_py_files = F)
  generate_sbatch(ps, scen = 'G', experiment = '003', runs = subset(cases, PS==ps)$run, unzip_py_files = F)
  generate_sbatch(ps, scen = 'G', experiment = '004', runs = subset(cases, PS==ps)$run, unzip_py_files = F)
}

# General stuff -----------------------------------------------------------

# Extract the random seeds of all experiments and runs
random_seeds <- c()
for (ps in ps_range){
  for (e in exp_range){
    for(r in 1:5){
      x <- data.frame(ps,e,r,seed=get_random_seed(PS = ps, scenario = 'N', run_range = r))
      random_seeds <- rbind(random_seeds, x)
    }
  }
}
y <- random_seeds %>% group_by(ps,r,seed) %>% summarise(s=length(seed))


for (ps in sprintf('%0.2d', 27:31)){
  for (e in sprintf('%0.3d', 1:4)){
    cat(paste('sbatch PS',ps,'NE',e,'.sbatch',sep=''));cat('\n')
  }
}



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

