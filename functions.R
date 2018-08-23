# General functions -------------------------------------------------------
on_Midway <- function(){
  if (Sys.getenv('USER')=='pilosofs'){
    return(TRUE)
  } else {
    return(FALSE)
  }
}
  
prep.packages <- function(package.list) {
  loaded = package.list %in% .packages()
  if ( all(loaded) ) return(invisible())
  
  package.list = package.list[!loaded]
  installed = package.list %in% .packages(TRUE)
  if ( !all(installed) ) install.packages(package.list[!installed], repos="http://cran.rstudio.com/")
  for ( p in package.list )
  {
    print(paste("Loading package:",p))
    library(p,character.only=TRUE)
  }
}

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

chunk2 <- function(x,n) {split(x, cut(seq_along(x), n, labels = FALSE))}

## A function used to get the repertoire names from unique repertoire copy
## names.
splitText <- function(str,after=T,splitchar='\\.'){
  if (after){
    return(sapply(strsplit(str, split=splitchar), tail, 1))
  }
  if (after==F){
    return(sapply(strsplit(str, split=splitchar), head, 1))
  }
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

# Requires package googlesheets
loadExperiments_GoogleSheets <- function(workBookName='malaria_interventions_design',sheetID=4){
  GS <- gs_title(workBookName)
  col_types <- GS %>% gs_read(ws=1, col_names=T)
  col_t <- unname(as.list(col_types[1,]))
  experiments <- GS %>% gs_read(ws=sheetID, col_names=T, col_types=col_t)
  print(experiments)
  return(experiments)
}

# Functions to set parameters ---------------------------------------------

# Function to get reference parameter file and output a data frame of parameters
get_parameter_reference <- function(parameter_file_ref='parameter_file_ref.py'){
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

# Function to set a parameter in the parameters data frame
set_parameter <- function(param_data, parameter, value){
  param_data$value[param_data$param==parameter] <- value
  return(subset(param_data, param==parameter))
}

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

# This function extracts the biting rates from the parameter file, which are in
# daily resolution and returns a vector of the averaged biting rates for a given
# period. (e.g. 30 would be biting rates averaged over a month and 1 will return
# the daily biting rates).
get_biting_rate <- function(parameter_file, sampling_period=30){
  x <- readLines(parameter_file)
  y <- x[grep('BITING_RATE_MEAN',x)[1]]
  BITING_RATE_MEAN <- parse_number(y)
  y <- x[grep('DAILY_BITING_RATE_DISTRIBUTION',x)[1]]
  DAILY_BITING_RATE_DISTRIBUTION <- eval(parse(text=paste('c(',(str_sub(y, str_locate(y, '\\[(.*?)\\]')[1]+1, str_locate(y, '\\[(.*?)\\]')[2]-1)),')',sep='')))
  BITING_RATE <- BITING_RATE_MEAN*DAILY_BITING_RATE_DISTRIBUTION
  BITING_RATE <- chunk2(BITING_RATE, 360/sampling_period)
  sapply(BITING_RATE, mean)
}

# This function creates a data frame with an experimental IRS design. This is a
# convenience function to avoid putting the desing in the Google sheet.
create_intervention_scheme_IRS <- function(PS_benchmark, scenario_benchmark, IRS_START_TIMES=NULL, immigration_range, length_range, coverage_range, poolsize_bounce='False', write_to_file=T, design_ref=NULL){
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
                        POOLSIZE_BOUNCE_BACK_AFTER_INTERVENTION=poolsize_bounce,
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
# Functions for neutral models parameters ---------------------------------

# This function obtains the duration of infection from the selection mode
# counterpart simulations to calculate the var switching rate in the neutral
# scenario. it makes more sense to take the duration from the control experiment
# (001), otherwise interventions can affect it, so I set exepriment=001 as
# default.
set_transition_rate_neutral <- function(parameter_space, run, N_GENES_PER_STRAIN = 60){
  if (on_Midway()){
    sqlite_file <- list.files(path = 'sqlite/', pattern=paste('PS',parameter_space,'_S_E001_R',run,'.sqlite',sep=''), full.names = T)
  } else {
    sqlite_file <- list.files(path = '~/Documents/malaria_interventions_data/', pattern=paste('PS',parameter_space,'_S_E001_R',run,'.sqlite',sep=''), full.names = T)
  }
  db <- dbConnect(SQLite(), dbname = sqlite_file)
  sampled_duration <- dbGetQuery(db, 'SELECT duration FROM sampled_duration')
  mean_doi <- mean(sampled_duration$duration)
  TRANSITION_RATE_NOT_IMMUNE <- 1/(mean_doi/N_GENES_PER_STRAIN)
  # setwd('~/Documents/malaria_interventions/')
  return(TRANSITION_RATE_NOT_IMMUNE)
}

# A function to obtain the transition rate from an existing checkpoint parameter
# file.
get_transition_rate_neutral <- function(PS, run){
  transition_rate <- c()
  x <- try(readLines(paste('PS',PS,'_N_','E000_R',run,'.py',sep='')))
  if (inherits(x, "try-error")){
    print(paste('Cannot find parameter file to extract transition rate! (','PS',PS,'_N_','E000_R',run,'.py)',sep=''))
    transition_rate <- NA
  } else {
    transition_rate <- parse_number(x[10])
  }
  return(transition_rate)
}

#A function to set the parameters for generazlied immunity scenario. Unlike the
#set_transition_rate_neutral, which returns a single value, this function
#returns a list which describes a fit of a curve of a particular run. It makes
#most sense to fit the curve to an experiment without interevntions and in
#stable state, so 001 is by default.
# REQUIRES rPython
set_generalized_immunity <- function(parameter_space, run){
  # Prepare the python file for the curve-fittign code
  if (on_Midway()){
    sqlite_file <- paste('sqlite/','PS',parameter_space,'_S_E001_R',run,'.sqlite',sep='')
    pyFile <- readLines('generalized_immunity_fitting.py')
  } else {
    sqlite_file <- paste('/home/shai/Documents/malaria_interventions_data/','PS',parameter_space,'_S_E001_R',run,'.sqlite',sep='')
    pyFile <- readLines('/home/shai/Documents/malaria_interventions/generalized_immunity_fitting.py')
  }
  pyFile[11] <- paste('path=','"',sqlite_file,'"',sep='')
  writeLines(pyFile,'generalized_immunity_fitting_experiment.py')
  
  # Try to run the code. Sometimes the curve does not fit
  res <- try(python.load('generalized_immunity_fitting_experiment.py'))  
  if(inherits(res, "try-error")){
    return(NA) # If fit is not sucessful return a NA
  } else {
    N_INFECTIONS_FOR_GENERAL_IMMUNITY <- python.get('infectionTimesToImmune')
    GENERAL_IMMUNITY_PARAMS <- python.get('generalImmunityParams')
    GENERAL_IMMUNITY_PARAMS <- paste('[',paste(GENERAL_IMMUNITY_PARAMS, collapse = ','),']',sep='')
    CLEARANCE_RATE_IMMUNE <- python.get('clearanceRateConstantImmune')
    file.remove('generalized_immunity_fitting_experiment.py')
    return(list(N_INFECTIONS_FOR_GENERAL_IMMUNITY=N_INFECTIONS_FOR_GENERAL_IMMUNITY,
                GENERAL_IMMUNITY_PARAMS=GENERAL_IMMUNITY_PARAMS,
                CLEARANCE_RATE_IMMUNE=CLEARANCE_RATE_IMMUNE))
  }
}

# A function to obtain the GI parameters from an existing checkpoint parameter
# file.
get_generalized_immunity <- function(PS, run){
  x <- try(readLines(paste('PS',PS,'_G_','E000_R',run,'.py',sep='')))
  if (inherits(x, "try-error")){
    print(paste('Cannot find parameter file to extract parameters! (','PS',PS,'_G_','E000_R',run,'.py)',sep=''))
    params <- NA
  } else {
    params <- list(N_INFECTIONS_FOR_GENERAL_IMMUNITY=parse_number(x[7]),
                   GENERAL_IMMUNITY_PARAMS=str_sub(x[8],25,str_length(x[8])),
                   CLEARANCE_RATE_IMMUNE=parse_number(x[9]))
  }
  return(params)
}

# This function sets the parameters for a single IRS. It is possible to run it
# several times, once for each IRS scheme. It uses a parameter file which has
# already has the rest of the parameters in place, and just changes the relevant
# parameters for IRS.
set_IRS <- function(design_ID, run, IRS_START_TIME, IRS_input, IRS_IMMIGRATION, POOLSIZE_BOUNCE_BACK_AFTER_INTERVENTION, experimental_design){
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
  # Add POOLSIZE_BOUNCE_BACK_AFTER_INTERVENTION
  param_data[param_data$param=='POOLSIZE_BOUNCE_BACK_AFTER_INTERVENTION',] <- set_parameter(param_data, 'POOLSIZE_BOUNCE_BACK_AFTER_INTERVENTION', POOLSIZE_BOUNCE_BACK_AFTER_INTERVENTION)
  
  
  # Write parameter file
  param_data$output <- paste(param_data$param,param_data$value,sep='=')
  write_lines(param_data$output, output_file)
}

# This function sets the parameters for a single MDA. It uses a parameter file which has
# already has the rest of the parameters in place, and just changes the relevant
# parameters for MDA.
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

# Function to create the necesary files and pipeline for a single run of an experiment.
# Each run has its own random seed across experiments and scenarios.
create_run <- function(design_ID, run, RANDOM_SEED, experimental_design, biting_rate_mathematica=NULL, params_GI=NULL){
  
  # Regime
  parameter_space <- experimental_design$PS[design_ID]
  scenario <- experimental_design$scenario[design_ID]
  experiment <- experimental_design$exp[design_ID] # Use 000 for checkpoints and 001 for control
  base_name <- paste('PS',parameter_space,'_',scenario,'_E',experiment,'_R',run,sep='')
  if (on_Midway()){
    param_data <- get_parameter_reference('parameter_file_ref.py')
  } else {
    param_data <- get_parameter_reference('~/Documents/malaria_interventions/parameter_file_ref.py')
  }
  
  # General parameters
  param_data[param_data$param=='RANDOM_SEED',] <- set_parameter(param_data, 'RANDOM_SEED', RANDOM_SEED)
  T_END <- experimental_design$T_END[design_ID]
  param_data[param_data$param=='T_END',] <- set_parameter(param_data, 'T_END', T_END)
  param_data[param_data$param=='VERIFICATION_ON',] <- set_parameter(param_data, 'VERIFICATION_ON', 'False')
  param_data[param_data$param=='VERIFICATION_PERIOD',] <- set_parameter(param_data, 'VERIFICATION_PERIOD', T_END)
  
  # Scenario
  if(scenario=='N'){
    param_data[param_data$param=='SELECTION_MODE',] <- set_parameter(param_data, 'SELECTION_MODE', "\'NEUTRALITY\'")
    # For each run, there is a need to calculate the transition rate once, for
    # the CP experiment. the rest of the experiments use the same value.
    if (experiment=='000'){ 
      print(paste('Calculating transition rate from PS',parameter_space,'_S_E001_R',run,sep = ''))
      TRANSITION_RATE_NOT_IMMUNE <- set_transition_rate_neutral(parameter_space, run)
    } else {
      print(paste('Getting transition rate from PS',parameter_space,'_N_E000_R',run,sep = ''))
      TRANSITION_RATE_NOT_IMMUNE <- get_transition_rate_neutral(parameter_space, run)
    }
    param_data[param_data$param=='TRANSITION_RATE_NOT_IMMUNE',] <- set_parameter(param_data, 'TRANSITION_RATE_NOT_IMMUNE', TRANSITION_RATE_NOT_IMMUNE)
  }
  if(scenario=='G'){
    param_data[param_data$param=='SELECTION_MODE',] <- set_parameter(param_data, 'SELECTION_MODE', "\'GENERAL_IMMUNITY\'")
    if (is.null(params_GI)){
      # For each run, there is a need to calculate the general immunity parameters once, for
      # the CP experiment. the rest of the experiments use the same values.
      if (experiment=='000'){ 
        print(paste('Calculating generalized immunity parameters from PS',parameter_space,'_S_E001_R',run,sep = ''))
        params_GI <- set_generalized_immunity(parameter_space, run)
      } else {
        print(paste('Getting generalized immunity parameters from PS',parameter_space,'_G_E000_R',run,sep = ''))
        params_GI <- get_generalized_immunity(parameter_space, run)
      }
      if (is.na(params_GI)[1]){
        print(paste('ERROR (create_run): Could not fit function for generalized immunity, skipping and NOT PRODUCING FILE ',base_name,'.py',sep = ''))
        return(NULL)
      }
      param_data[param_data$param=='GENERAL_IMMUNITY_PARAMS',] <- set_parameter(param_data, 'GENERAL_IMMUNITY_PARAMS', params_GI$GENERAL_IMMUNITY_PARAMS)
      param_data[param_data$param=='N_INFECTIONS_FOR_GENERAL_IMMUNITY',] <- set_parameter(param_data, 'N_INFECTIONS_FOR_GENERAL_IMMUNITY', params_GI$N_INFECTIONS_FOR_GENERAL_IMMUNITY)
      param_data[param_data$param=='CLEARANCE_RATE_IMMUNE',] <- set_parameter(param_data, 'CLEARANCE_RATE_IMMUNE', params_GI$CLEARANCE_RATE_IMMUNE)
    } else {
      print(paste('Using external GI parameters for file ',base_name,'.py',sep = ''))
      param_data[param_data$param=='GENERAL_IMMUNITY_PARAMS',] <- set_parameter(param_data, 'GENERAL_IMMUNITY_PARAMS', params_GI$GENERAL_IMMUNITY_PARAMS)
      param_data[param_data$param=='N_INFECTIONS_FOR_GENERAL_IMMUNITY',] <- set_parameter(param_data, 'N_INFECTIONS_FOR_GENERAL_IMMUNITY', params_GI$N_INFECTIONS_FOR_GENERAL_IMMUNITY)
      param_data[param_data$param=='CLEARANCE_RATE_IMMUNE',] <- set_parameter(param_data, 'CLEARANCE_RATE_IMMUNE', params_GI$CLEARANCE_RATE_IMMUNE)
    }
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

# This function generates parameter files and corresponding job files (to run on
# Midway), for several experiments and runs. It keeps the random seed for each
# RUN across experiments AND PARAMETER SPACES the same. Also possible to provide
# a random seed. row_range is the row numbers in the design data frame.
generate_files <- function(row_range, run_range, random_seed=NULL, experimental_design, biting_rate_file='mosquito_population_seasonality.csv', params_GI=NULL){
  
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
    # print(paste('Run: ',RUN,' | Row: ',design_ID,sep=''))
    
    # Create parameter file with main parameters. If the scenario is generalized
    # immunity then the fit for some runs may have not convereged (functions
    # 'create_run' and 'set_generalized_immunity'). Parameter files cannot be
    # produced in these cases and so the function skips to the next case.
    try_run <- create_run(design_ID, RUN, RANDOM_SEED, experimental_design, biting_rate_mathematica, params_GI=params_GI)
    if (is.null(try_run)){
      cases$file_created[idx] <- F
      next
    }
    
    # Set IRS
    if (!is.na(experimental_design$IRS_START_TIMES[design_ID])){
      IRS_scheme <- data.frame(IRS_START_TIME=str_split(experimental_design$IRS_START_TIMES[design_ID], ',')[[1]],
                               IRS_input=str_split(experimental_design$IRS_input[design_ID], ',')[[1]],
                               IRS_IMMIGRATION=str_split(experimental_design$IRS_IMMIGRATION[design_ID], ',')[[1]],
                               POOLSIZE_BOUNCE_BACK_AFTER_INTERVENTION=experimental_design$POOLSIZE_BOUNCE_BACK_AFTER_INTERVENTION[design_ID],
                               stringsAsFactors = F)
      for (i in 1:nrow(IRS_scheme)){
        set_IRS(design_ID = design_ID, 
                run = RUN, 
                IRS_START_TIME = IRS_scheme$IRS_START_TIME[i], 
                IRS_IMMIGRATION = IRS_scheme$IRS_IMMIGRATION[i], 
                IRS_input = IRS_scheme$IRS_input[i], 
                POOLSIZE_BOUNCE_BACK_AFTER_INTERVENTION=IRS_scheme$POOLSIZE_BOUNCE_BACK_AFTER_INTERVENTION[i], 
                experimental_design)
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
    if (on_Midway()){
      job_lines <- readLines('job_file_ref.sbatch')
    } else {
      job_lines <- readLines('~/Documents/malaria_interventions/job_file_ref.sbatch')
    }
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



# Plotting ----------------------------------------------------------------
library(ggplot2)
mytheme <- theme_bw() + theme(
  legend.title  = element_text(colour = "black", size=17),
  # legend.position = "none",
  #	legend.direction = "horizontal",
  legend.key = element_blank(),
  legend.text  = element_text(colour = "black", size=17),
  panel.background = element_blank(),
  # panel.grid.major = element_blank(),
  # panel.grid.minor = element_blank(),
  axis.text = element_text(color='black', family="Helvetica", size=14),
  axis.title = element_text(color='black', family="Helvetica", size=14),
  strip.text.x = element_text(family = "Helvetica", size = 14),
  strip.text.y = element_text(family = "Helvetica", size = 14),
  panel.border = element_rect(colour = "black", size=1.3),
  axis.ticks = element_line(size = 1.3),
  strip.background = element_rect( fill = "transparent", size = 1.3, colour = "black"  ),
  strip.text = element_text(size = 19)
)

gg_color_hue <- function(n, hue_min = 10, hue_max = 280, l = 62, c = 100) {
  hues = seq(hue_min, hue_max, length=n+1)
  hcl(h=hues, l=l, c=c)[1:n]
}

# This function uses the list obtained by get_data() and generates relevant plots
generate_plots <- function(data, time_range=NULL){
  data1 <- data[[1]]
  data2 <- data[[2]]
  if (is.null(time_range)){
    time_range <- c(min(data1$time), max(data1$time))
  }
  plot_variables <- data1 %>% 
    select(-n_infected, -year, -month, -PS, -exp, -run, -scenario) %>% 
    filter(time>time_range[1]&time<time_range[2]) %>% gather(variable, value, -time) %>% 
    ggplot(aes(time, value, color=variable))+
    geom_line()+
    facet_wrap(~variable, scales = 'free')+
    mytheme+theme(legend.position = 'none')
  
  plot_eir <- data1 %>% 
    ggplot(aes(x=month,y=EIR))+
    geom_boxplot()+
    geom_point(stat='summary', fun.y=mean, color='red')+
    stat_summary(fun.y=mean, geom="line")+mytheme
  
  # Plot the age structure of infected hosts
  plot_age_structure <- data2 %>% group_by()%>% ggplot(aes(x=host_age))+geom_histogram() + 
    labs(x='Infected host age (months)') + 
    geom_vline(xintercept = 60) +
    mytheme
  return(list(plot_variables, plot_eir, plot_age_structure))
}


plotLayer <- function(network_object, l, edge_weight_multiply=1, remove.loops=T, ver.col=NULL, coords,...) { 
  g <- network_object$temporal_network[[l]]
  g <- graph.adjacency(g, weighted = T, mode = 'directed')
  if(remove.loops){g <- simplify(g, remove.multiple = F, remove.loops = T)}
  # g <- delete_edges(g, which(E(g)$weight<quantile(E(g)$weight, cutoff_g))) # remove all edges smaller than the cutoff
  # layout <-layout.kamada.kawai(g)
  # plot.new()
  # par(mar=c(0,3,3,3))
  if (is.null(ver.col)){ # All nodes have the same color?
    V(g)$color <- 'deepskyblue2'
  } else {
    if (class(ver.col)=='character') {V(g)$color <- ver.col}
    if (class(ver.col)=='data.frame') {
      V(g)$color <- ver.col$color[match(V(g)$name, ver.col$node_name)]
    }
  }
  if (!is.null(coords)){
    V(g)$repertoire <- splitText(V(g)$name,splitchar = '_', after = F)
    V(g)$x <- coords$x[match(V(g)$repertoire, coords$node_name)]
    V(g)$y <- coords$y[match(V(g)$repertoire, coords$node_name)]
  }
  plot(g, 
       vertex.color=V(g)$color,
       # vertex.label=V(g)$module,
       vertex.size=3,
       vertex.label=NA,
       # edge.arrow.mode='-', 
       edge.arrow.width=0.2,
       edge.arrow.size=0.2,
       edge.curved=0.5, 
       edge.width=E(g)$weight*edge_weight_multiply,
       # asp=0, # This needed if changing margins
       ...)
}

# Data manipulation -------------------------------------------------------

# This function obtains data from an sqlite file and prepares them for further analysis.
# requires sqldf
get_data <- function(parameter_space, scenario, experiment, run, sampling_period=30, host_age_structure=F){
  # Initialize
  base_name <- paste('PS',parameter_space,'_',scenario,'_E',experiment,'_R',run,sep='')
  if (on_Midway()){
    sqlite_file <- paste('sqlite/',base_name,'.sqlite',sep='')
  } else {
    sqlite_file <- paste('~/Documents/malaria_interventions_data/',base_name,'.sqlite',sep='')
  }
  
  if (!file.exists(sqlite_file)) {
    print (paste(sqlite_file, ' does not exist, ignoring and terminating function'))
    return(NULL)
  }
  # parameter_file <- paste(base_name,'.py',sep='') # This may be necessary so I leave it
  
  # Extract data from sqlite. variable names correspond to table names
  db <- dbConnect(SQLite(), dbname = sqlite_file)
  summary_general <- dbGetQuery(db, 'SELECT * FROM summary')
  summary_general$PS <- parameter_space
  summary_general$exp <- experiment
  summary_general$scenario <- scenario
  summary_general$run <- run
  summary_general$year <- gl(n = max(summary_general$time)/360, length = nrow(summary_general), k = 1)
  summary_general$month <- gl(n = 12, k = 1, length = nrow(summary_general),labels = c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'), ordered = F)
  
  summary_alleles <- dbGetQuery(db, 'SELECT * FROM summary_alleles')
  summary_alleles %<>% group_by(time) %>% summarise(n_alleles=sum(n_circulating_alleles))
  
  summary_general <- suppressMessages(left_join(summary_general, summary_alleles))
  
  sampled_infections <- dbGetQuery(db, 'SELECT * FROM sampled_infections')
  sampled_infections$PS <- parameter_space
  sampled_infections$exp <- experiment
  sampled_infections$scenario <- scenario
  sampled_infections$run <- run
  
  # Extract statistics
  if (nrow(summary_general)%%sampling_period==1){# The beginning of the data set in E00 has a timestap of 0 and this line is unnecesary. So remove it
    summary_general <- summary_general[-1,]
  }
  
  ## Prevalence
  summary_general$prevalence <- summary_general$n_infected/10^4
  
  ## Get EIR from the table in the sqlite
  summary_general$EIR <- summary_general$n_infected_bites/10000 # 10000 is the size of the human population
  
  # ##EIR can also be calculated manually as EIR=biting_rate * prevalence
  # biting_rate <- get_biting_rate(parameter_file)
  # summary_general$b <- NA
  # summary_general$b[] <- biting_rate # The [] is for recycling the biting_rate
  # summary_general$EIR1 <- summary_general$prevalence*summary_general$b*sampling_period
  
  # # Annual EIR is given by dividing by the sampling period to get a daily EIR, then multiplying by 360
  # summary_general %>% mutate(eir_y=EIR2/sampling_period*360) %>% group_by(year) %>% summarise(EIR_Year=mean(eir_y))
  
  # MOI#
  meanMOI <- sampled_infections %>% group_by(time, host_id) %>% summarise(MOI=length(strain_id)) %>% group_by(time) %>% summarise(meanMOI=mean(MOI))
  summary_general <- suppressMessages(inner_join(summary_general, meanMOI)) # Add MOI to the summary
  
  
  # Host age structure
  if(host_age_structure) {
    hosts <- dbGetQuery(db, 'SELECT * FROM hosts')
    names(hosts)[1] <- 'host_id'
    hosts$lifespan <- round((hosts$death_time-hosts$birth_time))
    hosts <- subset(hosts, host_id%in%sampled_infections$host_id)
    sampled_infections <- suppressMessages(left_join(sampled_infections, hosts, by='host_id'))
    sampled_infections$host_age <- round((sampled_infections$time-sampled_infections$birth_time)/30)    
  }
  
  return(list(summary_general=as.tibble(summary_general), sampled_infections=as.tibble(sampled_infections)))
}





# This function compares an experiment to control and calculates statistics for the intervention
post_intervention_stats <- function(PS, scenario='S', exp, run, post_intervention_lag=360, control_data=NULL, plot.it=F, design_irs){
  
  plots_out <- list()
  
  # If a data frame with control data is not included (NULL), then load the
  # control experiment
  if (is.null(control_data)){
    control_data <- get_data(parameter_space = PS, scenario = scenario, experiment = '001', run = run)[[1]]
  }
  
  # Calculate variable means (across time) for control
  control_means <- control_data %>% select(-year, -month, -n_infected) %>% 
    gather(variable, value, -time, -exp, -PS, -scenario, -run, -pop_id) %>% 
    group_by(PS, exp, variable) %>% summarise(mean_value=mean(value, na.rm = T))
  
  # Plot control and means
  if (plot.it){
    plots_out$control <- control_data %>% select(-year, -month, -n_infected) %>%
      gather(variable, value, -time, -exp, -PS, -scenario, -run, -pop_id) %>%
      ggplot(aes(x=time, y=value))+
      geom_line()+
      facet_wrap(~variable,scales='free')+
      geom_hline(aes(yintercept=mean_value), data=control_means, color='blue')+
      mytheme
  }
  
  # Load experiment data
  experiment_data <- get_data(parameter_space = PS, scenario = scenario, experiment = exp, run = run)[[1]]
  
  # Plot control vs experiment
  if (plot.it){
    plots_out$control_experiment <- bind_rows(control_data, experiment_data) %>%
      select(-year, -month, -n_infected) %>% 
      gather(variable, value, -time, -exp, -PS, -scenario, -run, -pop_id) %>% 
      ggplot(aes(x=time, y=value, color=exp))+
      geom_vline(xintercept = c(29160,29160+1800), linetype='dashed', color='tomato')+
      geom_line()+
      scale_color_manual(values=c('#504E4E','purple'))+
      scale_x_continuous(breaks=pretty(x=subset(d, time>time_range[1]&time<time_range[2])$time,n=5))+
      geom_hline(aes(yintercept=mean_value), data=control_means, color='tomato')+
      facet_wrap(~variable, scales = 'free')+mytheme
  }
  
  # Calculate absolute difference from control at every timepoint
  x <- control_data %>%
    select(-year, -month, -n_infected) %>% 
    gather(variable, value, -time, -exp, -PS, -scenario, -run, -pop_id)
  y <- experiment_data %>%
    select(-year, -month, -n_infected) %>% 
    gather(variable, value, -time, -exp, -PS, -scenario, -run, -pop_id)
  diff_control <- suppressMessages(inner_join(x,y,by=c('pop_id','PS','scenario','run','time','variable'))) %>% 
    mutate(diff=value.y-value.x, ratio=value.y/value.x) %>% 
    select(variable, time, diff, ratio)
  
  if (plot.it){
    plots_out$control_diff <- diff_control %>%
      ggplot(aes(time, diff))+
      geom_line()+
      facet_wrap(~variable,scales='free')+
      geom_hline(aes(yintercept=mean_value), data=control_means, color='blue')+
      geom_vline(xintercept = c(29160,29160+3600), color='red')+
      mytheme
  }
  
  # Maximum amplitudes after an intervention
  amplitude <- experiment_data %>% left_join(subset(design_irs, select=c(IRS_START_TIMES,IRS_length,exp)), by='exp') %>% 
    mutate(intervention_lift=as.numeric(IRS_START_TIMES)+as.numeric(IRS_length)) %>% 
    filter(time>intervention_lift) %>% 
    select(-year, -month, -n_infected) %>% 
    gather(variable, value, -time, -exp, -PS, -scenario, -run, -pop_id,-IRS_START_TIMES,-IRS_length) %>% 
    group_by(variable) %>% summarise(max_value=max(value))
  
  # New mean from some lag post intervention
  new_mean <- experiment_data %>% left_join(subset(design_irs, select=c(IRS_START_TIMES,IRS_length,exp)), by='exp') %>% 
    mutate(intervention_lift=as.numeric(IRS_START_TIMES)+as.numeric(IRS_length)) %>% 
    filter(time>intervention_lift+post_intervention_lag) %>% 
    select(-year, -month, -n_infected) %>% 
    gather(variable, value, -time, -exp, -PS, -scenario, -run, -pop_id,-IRS_START_TIMES,-IRS_length) %>% 
    group_by(variable) %>% summarise(mean_value=mean(value, na.rm = T))
  
  # Time when extinction occured (last time point with data)
  time_extinct <- experiment_data %>% left_join(subset(design_irs, select=c(IRS_START_TIMES,IRS_length,exp)), by='exp') %>% 
    select(-year, -month, -n_infected) %>% 
    gather(variable, value, -time, -exp, -PS, -scenario, -run, -pop_id,-IRS_START_TIMES,-IRS_length) %>% 
    group_by(variable) %>% summarise(time_extinct=max(time))
  
  summary_stats <- suppressMessages(left_join(time_extinct, new_mean))
  summary_stats <- suppressMessages(left_join(summary_stats, amplitude))
  summary_stats$PS <- PS
  summary_stats$scenario <- scenario
  summary_stats$exp <- exp
  summary_stats$run <- run  
  
  diff_control$PS <- PS
  diff_control$scenario <- scenario
  diff_control$exp <- exp
  diff_control$run <- run 
  
  if (plot.it){
    return(list(plots=plots_out,stats=list(diff_control=diff_control, summary_stats=summary_stats)))
  } else {
    return(list(diff_control=diff_control, summary_stats=summary_stats))  
  }
}



# Functions for generating networks ---------------------------------------

# Calculate edge values in networks of repertories
overlapAlleleAdj<-function(mat){
  newmat<-tcrossprod(mat>0)
  newmat<-newmat/rowSums(mat>0)
  return(newmat)  
}


# A function to build the similarity matrix for a single layer and calculate some summary stats
build_layer <- function(infection_df){
  infection_df %<>% group_by(strain_id) %>%
    mutate(freq = n()/120) %>% # strain frequency (the number of strain copies should be equal to the frequency)
    arrange(strain_id_unique) 
  # Calculate the edges
  similarity_matrix <- table(infection_df$strain_id_unique, infection_df$allele_locus)
  similarity_matrix <- overlapAlleleAdj(similarity_matrix)
  # Some summary
  layer_summary <- with(infection_df, 
                        data.frame(hosts=length(unique(host_id)),
                                   repertoires_unique=length(unique(strain_id)),
                                   repertoires_total=length(unique(strain_id_unique))
                        ))
  return(list(similarity_matrix=similarity_matrix, infections=infection_df, layer_summary=layer_summary))
}


createTemporalNetwork <- function(ps, scenario, exp, run, cutoff_prob=0.9, cutoff_value=NULL, write_files=F, layers_to_include=NULL, sampled_infections=NULL){
  base_name <- paste('PS',ps,'_',scenario,'_E',exp,'_R',run,sep='')
  if (on_Midway()){
    sqlite_file <- paste('sqlite/',base_name,'.sqlite',sep='')
  } else {
    sqlite_file <- paste('~/Documents/malaria_interventions_data/',base_name,'.sqlite',sep='')
  }
  
  # Extract data from sqlite. variable names correspond to table names
  db <- dbConnect(SQLite(), dbname = sqlite_file)
  print('Getting genetic data from sqlite...')
  sampled_strains <- as.tibble(dbGetQuery(db, 'SELECT id, gene_id FROM sampled_strains'))
  names(sampled_strains)[1] <- c('strain_id')
  sampled_alleles <- as.tibble(dbGetQuery(db, 'SELECT * FROM sampled_alleles'))
  names(sampled_alleles)[3] <- c('allele_id')
  sampled_strains <- full_join(sampled_strains, sampled_alleles)
  sampled_strains$allele_locus <- paste(sampled_strains$allele_id,sampled_strains$locus,sep='_') # each allele in a locus is unique
  dbDisconnect(db)
  # Get infection data
  if (is.null(sampled_infections)){
    print('Getting infection data from sqlite...')
    sampled_infections <- get_data(parameter_space = ps, experiment = exp, scenario = scenario, run = run)$sampled_infections
  }
  print('Building data set...')
  # Add layer numbers. Each time point is a layer
  sampled_infections$layer <- group_indices(sampled_infections, time) 
  sampled_infections %<>% arrange(layer,strain_id,host_id) %>% 
    group_by(layer,strain_id) %>% 
    mutate(strain_copy = row_number()) # add unique id for each strain copy within each layer. A strain copy is an instance of a strain in a particualr host
  sampled_infections$strain_id_unique <- paste(sampled_infections$strain_id,sampled_infections$strain_copy,sep='_') # Create the unique strains
  # Integrate the strain composition into the infections table
  if (all(unique(sampled_strains$strain_id)%in%unique(sampled_infections$strain_id))==F || all(unique(sampled_infections$strain_id)%in%unique(sampled_strains$strain_id))==F) {
    stop('There may be a mismatch in repertoires between the sampled_strains and sampled_infections data sets. Revise!')
  }
  sampled_infections %<>% select(time, layer, host_id, strain_id, strain_copy, strain_id_unique) %>% 
    left_join(sampled_strains)
  
  print('Building layers...')
  Layers <- list()
  if (is.null(layers_to_include)){
    layers_to_include <- sort(unique(sampled_infections$layer))
  }
  for (l in layers_to_include){
    cat(paste('[',Sys.time(), '] building layer ',l,'\n',sep=''))
    sampled_infections_layer <- subset(sampled_infections, layer==l) # This is MUCH faster than sampled_infections_layer <- sampled_infections %>% filter(layer==l)
    Layers[[which(layers_to_include==l)]] <- build_layer(sampled_infections_layer)
  }
  # Get just the matrices
  temporal_network <- lapply(Layers, function(x) x$similarity_matrix)
  layer_summary <- do.call(rbind, lapply(Layers, function(x) x$layer_summary))
  layer_summary$layer <- layers_to_include
  sapply(temporal_network, nrow)
  
  print('Applying cutoff...')
  # Apply a cutoff
  ## Get the similarity matrix for all the repertoires and alleles
  x <- xtabs(~strain_id+allele_locus, subset(sampled_strains, strain_id%in%sampled_infections$strain_id))
  similarityMatrix <- overlapAlleleAdj(x)
  dim(similarityMatrix)
  
  if (is.null(cutoff_value)){
    cutoff_value <- quantile(as.vector(similarityMatrix), probs = cutoff_prob)
  }
  print(cutoff_value)
  # hist(as.vector(similarityMatrix));abline(v=cutoff_value,col='red')
  similarityMatrix[similarityMatrix<cutoff_value] <- 0
  for (i in 1:length(temporal_network)){
    # print(i)
    x <- temporal_network[[i]]
    x[x<cutoff_value] <- 0
    temporal_network[[i]] <- x
  }
  
  if(write_files){
    print('Writing similarity matrices to files...')
    # Write similarity matrix without cutoff. This is used to plot the edge weight distributions.
    write.table(similarityMatrix, paste('../',filenameBase,'_similarityMatrix_nocutoff.csv',sep=''), row.names = T, col.names = T, sep=',')
    # write.table(similarityMatrix, paste(filenameBase,'_similarityMatrix.csv',sep=''), row.names = T, col.names = T, sep=',')
    # write.table(similarityMatrix, paste('../',filenameBase,'_similarityMatrix.csv',sep=''), row.names = T, col.names = T, sep=',')
  }
  
  print('Done!')
  return(list(temporal_network=temporal_network, similarityMatrix=similarityMatrix, cutoff_prob=cutoff_prob, cutoff_value=cutoff_value, layer_summary=layer_summary, base_name = base_name, ps=ps, scenario=scenario, experiment=exp, run=run))
}







# Functions for network properties ----------------------------------------
f_01_averageLocalClusteringCoeff <- function(g,GC=F){
  if (GC) {
    gc <- giant.component(g)
    return(mean(transitivity(gc, type = 'local'), na.rm = T))
  } else {
    return(mean(transitivity(g, type = 'local'), na.rm = T))
  }
}

f_02_averageLocalClusteringCoeffWeighted <- function(g,GC=F){
  if (GC) {
    gc <- giant.component(g)
    return(mean(transitivity(gc, type = 'barrat'), na.rm = T))
  } else {
    return(mean(transitivity(g, type = 'barrat'), na.rm = T))
  }
}


f_03_globalClusteringCoeff <- function(g, GC=F){
  if (GC) {
    gc <- giant.component(g)
    return(transitivity(gc, type = 'global'))
  } else {
    return(transitivity(g, type = 'global'))
  }
}

f_04_gdensity <- function(g){
  return(graph.density(g))
}

f_05_proportionSingletons <- function(g){
  sum(degree(g)==0)/length(V(g))
}

f_06_proportionEndpoints <- function(g){
  sum(degree(g)==1)/length(V(g))
}

f_07_meanDegree <- function(g){
  mean(degree(g))
}

f_08_meanDegreeNotSingletons <- function(g){
  mean(degree(g)[degree(g)!=0])
}

f_09_meanStrength <- function(g){
  mean(strength(g))
}

f_10_meanStrengthNotSingletons <- function(g){
  mean(strength(g)[degree(g)!=0])
}

f_11_entropyDegreeDistribution <-  function(g,verbose=F){
  y=degree(g)
  freq=prop.table(table(y))
  if (verbose){print(freq)}
  -sum(freq * log(freq, base = 2))
}

f_12_ratioComponents <- function(g) { 
  cl <- clusters(g) 
  cl$no/length(V(g))
}

giant.component <- function(g) { 
  cl <- clusters(g) 
  induced.subgraph(g, which(cl$membership == which.max(cl$csize)))
}

f_13_averageComponentSize <- function(g) { 
  cl <- clusters(g) 
  mean(cl$csize)
}

f_14_entropyComponentSize <-  function(g,verbose=F){
  cl <- clusters(g)
  y=cl$csize
  freq=prop.table(table(y))
  if (verbose){print(freq)}
  -sum(freq * log(freq, base = 2))
}

f_15_giantConnectedRatio <- function(g){
  cl <- clusters(g)
  max(cl$csize)/length(V(g))
}

f_16_meanEccentricity <- function(g){
  mean(eccentricity(g))
}

f_17_gdiameter <- function(g){
  return(diameter(g))
}

f_18_meanDiameterComponents <- function(g){
  cl <- clusters(g)
  d <- 0
  for (m in 1:(cl$no)){
    d <- d + diameter(induced.subgraph(g, which(cl$membership == m)))
  }
  return(d/cl$no)
}

f_19_globalEfficiency <- function(g){
  d_ij <- shortest.paths(g)
  d_ij <- d_ij[lower.tri(d_ij)]
  d_ij <- d_ij[!is.infinite(d_ij)]
  N=length(V(g))
  1/(N*(N-1))*sum(1/d_ij)
}

f_20_averageClosenessCentrality <- function(g){
  mean(closeness(g, weights = NULL))
}

graphDensityDirected <- function(m){ # density of directed graphs
  E = sum(m!=0)-length(diag(m)) # in our networks the diagonal is 1 so this should be considered.
  N=nrow(m)
  return(E/(N*(N-1)))
}




motifsProportion <- function(g){ # Calculate the proportion of each of the 16 motifs out of the total motifs found
  
  # This code shows the motifs by order
  # par(mfrow=c(4,4), mar=c(.75,.75,.75,.75))
  # for (i in 0:15){ # note that counting starts at 0
  #   g <- graph_from_isomorphism_class(size = 3, number = i)
  #   V(g)$name <- c('B','A','C')
  #   plot(g,
  #        edge.arrow.size = 0.4,
  #        edge.color='black',
  #        main = i + 1)
  #   box(col='red')
  # }
  
  motifs <- graph.motifs(g, size = 3)
  motifs.prop <- motifs/sum(motifs, na.rm = T)
  names(motifs.prop) <- c('M01--A_B_C', #1
                          'M02--A->B_C', #2
                          'M03--A->B<-C', #3
                          'M04--A<->B_C', #4
                          'M05--A->B->C', #5
                          'M06--A<->B<-C',#6
                          'M07--A<-B->C', #7
                          'M08--A->B<-C_A->C', #8
                          'M09--A<-B->C_A<->C', #9
                          'M10--A<->B->C', #10
                          'M11--A<->B<->C', #11
                          'M12--A<-B<-C_A->C', #12
                          'M13--A->B->C_A<->C', #13
                          'M14--A->B<-C_A<->C', #14
                          'M15--A->B<->C_A<->C', #15
                          'M16--A<->B<->C_A<->C') #16
  return(motifs.prop)
}


calculateFeatures <- function(network_object, l, remove.loops=F){
  g <- network_object$temporal_network[[l]]
  g <- graph.adjacency(g, weighted = T, mode = 'directed')
  if(remove.loops){g <- simplify(g, remove.multiple = F, remove.loops = T)}
  
  featureVector <- NULL
  featureVector <- c(featureVector, vcount(g))
  featureVector <- c(featureVector, f_03_globalClusteringCoeff(g))
  featureVector <- c(featureVector, f_04_gdensity(g))
  featureVector <- c(featureVector, f_07_meanDegree(g))
  featureVector <- c(featureVector, f_12_ratioComponents(g))
  featureVector <- c(featureVector, f_15_giantConnectedRatio(g))
  featureVector <- c(featureVector, f_17_gdiameter(g))
  featureVector <- c(featureVector, motifsProportion(g))
  names(featureVector)[1:7] <- c('Num_nodes','GCC','density','mean_degree','ratio_comp','nodes_in_giant_comp','diameter')
  
  # featureVector <- vector(length=32)
  # Diagnostics of transitivity
  # featureVector[1] <- f_01_averageLocalClusteringCoeff(g,F)    # Clustering coefficient averaged across all nodes
  # featureVector[2] <- f_02_averageLocalClusteringCoeffWeighted(g,F)    # Barrat's clustering coefficient averaged across all nodes
  # featureVector[3] <- f_03_globalClusteringCoeff(g,F)
  # Diagnostics of degree/sterngth
  # featureVector[4] <- f_04_gdensity(g)                    # Graph density
  # featureVector[5] <- f_05_proportionSingletons(g)             # Proportion of nodes with degree 0 of all the nodes
  # featureVector[6] <- f_06_proportionEndpoints(g)              # Proportion of nodes with degree 1 of all the nodes
  # featureVector[7] <- f_07_meanDegree(g)                       # Average degree 
  # featureVector[8] <- f_08_meanDegreeNotSingletons(g)
  # featureVector[9] <- f_09_meanStrength(g)                     # Average strength
  # featureVector[10] <- f_10_meanStrengthNotSingletons(g)
  # featureVector[11] <- f_11_entropyDegreeDistribution(g)       # Average measurement of the heterogeneity of the network. See eq. 14 in: da F. Costa, L., Rodrigues, F. A., Travieso, G. & Boas, P. R. V. Characterization of complex networks: A survey of measurements. Adv. Phys. 56, 167–242 (2007).
  # featureVector[12] <- f_12_ratioComponents(g)                 # Number of components relative to networks size
  # featureVector[13] <- f_13_averageComponentSize(g)
  # featureVector[14] <- f_14_entropyComponentSize(g)
  # featureVector[15] <- f_15_giantConnectedRatio(g)             # Proportion of nodes in the giant component
  # #Diagnostics of shortest-paths
  # featureVector[16] <- f_16_meanEccentricity(g)                 # Eccentricity is the maximum shortest distance from a node to all other nodes. This is averaged across all nodes
  # featureVector[17] <- f_17_gdiameter(g)                        # length of the longest geodesic
  # featureVector[18] <- f_18_meanDiameterComponents(g)
  # featureVector[19] <- f_19_globalEfficiency(g)                # See eq. 14 in: da F. Costa, L., Rodrigues, F. A., Travieso, G. & Boas, P. R. V. Characterization of complex networks: A survey of measurements. Adv. Phys. 56, 167–242 (2007).
  # featureVector[20] <- f_20_averageClosenessCentrality(g)      # 
  # # Diagnostics of motifs
  # featureVector[21:32] <- motifsProportion(g)
  return(featureVector)
}


# Infomap -----------------------------------------------------------------


##  A function that gets the layer as a matrix and writes it for infomap as an edge list
# network_object is a list of matrices, each element in the list is a layer.
# requires igraph
matrix_to_infomap <- function(l, nodeList, network_object){
  current_layer <- network_object$temporal_network[[l]]
  if(nrow(current_layer)<2){
    print(paste('Less than 2 repertoires in layer',l,'!!! skipping intralayer edges'))
    return(NULL)
  }
  g <- graph.adjacency(current_layer, mode = 'directed', weighted = T)
  current_layer_el <- as.tibble(as_data_frame(g, what = 'edges'))
  names(current_layer_el) <- c('node_s','node_t','w')
  current_layer_el$layer_s <- l
  current_layer_el$layer_t <- l
  current_layer_el$node_s <- nodeList$nodeID[match(current_layer_el$node_s,nodeList$nodeLabel)]
  current_layer_el$node_t <- nodeList$nodeID[match(current_layer_el$node_t,nodeList$nodeLabel)]
  current_layer_el %<>% select(layer_s, node_s, layer_t, node_t, w) # Re-arrange columns for the infomap input order
  print(paste('[',Sys.time(), '] Created edge list of layer ',l,' for Infomap | ', nrow(current_layer_el),' edges',sep=''))
  return(current_layer_el)
}

## A function that calculates interlayer edges between layer t and t+1
build_interlayer_edges_1step <- function(t, nodeList, network_object){
  strain_copies_t <- rownames(network_object$temporal_network[[t]]) # repertoires at time t
  strain_copies_t1 <- rownames(network_object$temporal_network[[t+1]]) # repertoires at time t+1
  # need minimum of 2 strains in t and t+1 to build a matrix
  if (length(strain_copies_t)<2 | length(strain_copies_t1)<2){
    print(paste('No interlayer edges between layers',t, 'and',t+1,'because there are < 2 repertoires.'))
    return(NULL)
  } 
  # Pull the similarity values between the repertoires from the general similarity matrix, then write back the repertoire copy names
  inter_layer_edges_matrix <- network_object$similarityMatrix[splitText(strain_copies_t, after = F, '_'),splitText(strain_copies_t1, after = F, '_')] 
  rownames(inter_layer_edges_matrix) <- strain_copies_t
  colnames(inter_layer_edges_matrix) <- strain_copies_t1
  if (all(inter_layer_edges_matrix==0)){
    print(paste('All interlayer edges between layers',t, 'and',t+1,' are 0 (due to the cutoff). skipping.'))
    return(NULL)
  }
  # transform to an edge list
  g <- graph.incidence(inter_layer_edges_matrix, directed = T, mode = 'out', weighted = T)
  inter_layer_edges <- as.tibble(as_data_frame(g, what = 'edges'))
  names(inter_layer_edges) <- c('node_s','node_t','w')
  inter_layer_edges$layer_s <- t
  inter_layer_edges$layer_t <- t+1
  inter_layer_edges$node_s <- nodeList$nodeID[match(inter_layer_edges$node_s,nodeList$nodeLabel)]
  inter_layer_edges$node_t <- nodeList$nodeID[match(inter_layer_edges$node_t,nodeList$nodeLabel)]
  inter_layer_edges %<>% select(layer_s, node_s, layer_t, node_t, w)
  print(paste('Built interlayer edges for layers: ',t,'-->',t+1,' | ', nrow(inter_layer_edges),' edges',sep=''))
  return(inter_layer_edges)
}


# This function takes a list of temporal matrices and returns an intralayer and
# interlayer edge lists in a format: [layer_source, node_source, layer_target node_target, weight].
# It also returns the list of node names.
# requires igraph
build_infomap_objects <- function(network_object, write_to_infomap_file=T, return_objects=T){
  temporal_network <- network_object$temporal_network
  base_name <- network_object$base_name

  # Get the node list
  nodeLabel <- sort(unique(unlist(lapply(temporal_network,rownames))))
  nodeList <- data.frame(nodeID=1:length(nodeLabel), nodeLabel)
  
  # Create intralayer edge lists
  layers <- 1:length(temporal_network)
  infomap_intralayer <- lapply(layers, function (x) matrix_to_infomap(x, nodeList = nodeList, network_object = network_object))
  infomap_intralayer <- do.call("rbind", infomap_intralayer)
  
  # Create interlayer edge lists. Only connect consecutive layers (interlayer edges ONLY go from layer l to l+1).
  layers <- layers[-length(layers)]
  infomap_interlayer <- lapply(layers, function (x) build_interlayer_edges_1step(x, nodeList = nodeList, network_object = network_object))
  infomap_interlayer <- do.call("rbind", infomap_interlayer)
  
  if (write_to_infomap_file){
    ## Write file for infomap
    print('Writing Infomap files')
    if (on_Midway()){
      file <- paste('Results/',job_ps,'_',job_scenario,'/',base_name,'_Infomap_multilayer','.txt',sep='')
    } else {
      file <- paste(base_name,'_Infomap_multilayer','.txt',sep='')
    }
    print(paste('Infomap file:',file))
    if (file.exists(file)){unlink(file)}
    sink(file, append = F)
    cat("# A network in a general multiplex format");cat('\n')
    cat(paste("*Vertices",nrow(nodeList)));cat('\n')
    write.table(nodeList, file, append = T,sep=' ', quote = T, row.names = F, col.names = F)
    cat("*Multiplex");cat('\n')
    cat("# layer node layer node [weight]");cat('\n')
    cat("# Intralayer edges");cat('\n')
    write.table(infomap_intralayer, file, sep = ' ', row.names = F, col.names = F, quote = F, append = T)
    cat("# Interlayer edges");cat('\n')
    write.table(infomap_interlayer, file, sep = ' ', row.names = F, col.names = F, quote = F, append = T)
    sink.reset()
  }
  
  if (return_objects){
    return(list(infomap_intralayer=infomap_intralayer, infomap_interlayer=infomap_interlayer, nodeList=nodeList))
  }
}


