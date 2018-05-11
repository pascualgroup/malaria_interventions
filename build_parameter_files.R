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


# Define parameters to set ------------------------------------------------
experiment <- 'CPsave_test'
BITING_RATE_MEAN <- 1
N_GENES_INITIAL <- 1200 # The initial gene pool

# The biting rate distribution for seasonality and/or IRS is generated in
# Mathematica and exported to a file. This distribution consists of the number
# of adult mosquitos. This is done per experiment.
mathematica_file <- 'mosquito_population_seasonality.csv' 
biting_rate_mathematica <- read_csv(mathematica_file, col_names = c('day','num_mosquitos'))
DAILY_BITING_RATE_DISTRIBUTION <- biting_rate_mathematica$num_mosquitos

DAILY_BITING_RATE_DISTRIBUTION <- paste(rep(0.01,360),collapse=',') # This is to set a fixed biting rate. Mostly for testing.

# Set parameters ----------------------------------------------------------
# File names
param_data[param_data$param=='SAMPLE_DB_FILENAME',] <- set_parameter(param_data, 'SAMPLE_DB_FILENAME', paste('\'\"',experiment,'.sqlite\"\'',sep=''))
# Genetic diversity
param_data[param_data$param=='N_GENES_INITIAL',] <- set_parameter(param_data, 'N_GENES_INITIAL', N_GENES_INITIAL)
param_data[param_data$param=='N_ALLELES_INITIAL',] <- set_parameter(param_data, 'N_ALLELES_INITIAL', paste('N_LOCI*[',N_GENES_INITIAL/10,']',sep=''))
# Biting rates
param_data[param_data$param=='BITING_RATE_MEAN',] <- set_parameter(param_data, 'BITING_RATE_MEAN', paste('[',BITING_RATE_MEAN,']',sep=''))
param_data[param_data$param=='DAILY_BITING_RATE_DISTRIBUTION',] <- set_parameter(param_data, 'DAILY_BITING_RATE_DISTRIBUTION', paste('[',paste(DAILY_BITING_RATE_DISTRIBUTION, collapse=','),']',sep=''))
# Timing
param_data[param_data$param=='T_END',] <- set_parameter(param_data, 'T_END', 2880)
param_data[param_data$param=='VERIFICATION_PERIOD',] <- set_parameter(param_data, 'VERIFICATION_PERIOD', 2880)
param_data[param_data$param=='CHECKPOINT_SAVE_PERIOD',] <- set_parameter(param_data, 'CHECKPOINT_SAVE_PERIOD', 2880)
# Checkpoints
param_data[param_data$param=='SAVE_TO_CHECKPOINT',] <- set_parameter(param_data, 'SAVE_TO_CHECKPOINT', 'True')
param_data[param_data$param=='CHECKPOINT_SAVE_FILENAME',] <- set_parameter(param_data, 'CHECKPOINT_SAVE_FILENAME', paste('\'\"',experiment,'_checkpoint.sqlite\"\'',sep=''))




# Write to a new paramter file --------------------------------------------
output_file=paste(experiment,'py',sep = '.')
param_data$output <- paste(param_data$param,param_data$value,sep='=')
write_lines(param_data$output, output_file)


# Write Midway job file ---------------------------------------------------
output_file=paste(experiment,'sbatch',sep = '.')
job_lines <- readLines('job_file_ref.sbatch')
wall_time <- '01:00:00'
SLURM_ARRAY_RANGE <- '1'
memory <- '4000'
job_lines[2] <- paste('#SBATCH --job-name=',experiment,sep='')
job_lines[3] <- paste('#SBATCH --time=',wall_time,sep='')
job_lines[4] <- paste('#SBATCH --output=slurm_output/',experiment,'_%A_%a.out',sep='')
job_lines[5] <- paste('#SBATCH --error=slurm_output/',experiment,'_%A_%a.err',sep='')
job_lines[6] <- paste('#SBATCH --array=',SLURM_ARRAY_RANGE,sep='')
job_lines[9] <- paste('#SBATCH --mem-per-cpu=',memory,sep='')
job_lines[20] <- paste("exp='",experiment,"'",sep='')
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

