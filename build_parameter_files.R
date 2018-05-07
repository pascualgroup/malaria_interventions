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


# Set parameters ----------------------------------------------------------
experiment <- 'test_09'
mathematica_file <- 'mosquito_population_IRS02.csv' # This file contains the number of adult mosquitos for the particular experiment
param_data[param_data$param=='SAMPLE_DB_FILENAME',] <- set_parameter(param_data, 'SAMPLE_DB_FILENAME', paste('\'\"',experiment,'.sqlite\"\'',sep=''))
param_data[param_data$param=='BITING_RATE_MEAN',] <- set_parameter(param_data, 'BITING_RATE_MEAN', '[0.00005]')
biting_rate_mathematica <- read_csv(mathematica_file, col_names = c('day','num_mosquitos'))
biting_rate_distribution <- biting_rate_mathematica$num_mosquitos
param_data[param_data$param=='DAILY_BITING_RATE_DISTRIBUTION',] <- set_parameter(param_data, 'DAILY_BITING_RATE_DISTRIBUTION', paste('[',paste(biting_rate_distribution, collapse=','),']',sep=''))


# Write to a new paramter file --------------------------------------------
output_file=paste(experiment,'py',sep = '.')
param_data$output <- paste(param_data$param,param_data$value,sep='=')
write_lines(param_data$output, output_file)


# Write Midway job file ---------------------------------------------------
output_file=paste(experiment,'sbatch',sep = '.')
job_lines <- readLines('job_file_ref.sbatch')
wall_time <- '20:00:00'
SLURM_ARRAY_RANGE <- '1'
memory <- '16000'
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
