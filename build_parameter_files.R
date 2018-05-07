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
param_data[param_data$param=='SAMPLE_DB_FILENAME',] <- set_parameter(param_data, 'SAMPLE_DB_FILENAME', '\'\"test_05.sqlite\"\'')
param_data[param_data$param=='BITING_RATE_MEAN',] <- set_parameter(param_data, 'BITING_RATE_MEAN', '[0.0002]')
x <- read_csv('mosquito_population_seasonality.csv', col_names = c('day','num_mosquitos'))
param_data[param_data$param=='DAILY_BITING_RATE_DISTRIBUTION',] <- set_parameter(param_data, 'DAILY_BITING_RATE_DISTRIBUTION', paste('[',paste(x$num_mosquitos, collapse=','),']',sep=''))


# Write to a new paramter file --------------------------------------------
output_file='test_05.py'
param_data$output <- paste(param_data$param,param_data$value,sep='=')
write_lines(param_data$output, output_file)

