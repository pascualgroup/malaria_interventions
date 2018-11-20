# Initialize --------------------------------------------------------------
source('functions.R')
library(sqldf)
library(tidyverse)
library(magrittr)
library(igraph)
library(data.table)

if (length(commandArgs(trailingOnly=TRUE))==0) {
  args <- c('03','N','001',0.9)
} else {
  args <- commandArgs(trailingOnly=TRUE)
}
job_ps <- as.character(args[1]) 
job_scenario <- as.character(args[2]) 
job_exp <- as.character(args[3]) 
# job_run <- 1
job_run <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
cutoff_prob <- as.numeric(args[4])

base_name <- paste('PS',job_ps,'_',job_scenario,'_E',job_exp,'_R',job_run,'_',cutoff_prob,sep='')

# Data from sqlite
data <- get_data(parameter_space = job_ps, scenario = job_scenario, experiment = job_exp, run = job_run, host_age_structure = T, use_sqlite = T)
write_csv(data$summary_general, paste(base_name,'_summary_general.csv',sep=''))
write_csv(data$sampled_infections, paste(base_name,'_sampled_infections.csv',sep=''))

# Network objects
# If experiment is not control then take the cutoff value from the control. This
# requires the control to be run first.
if (job_exp=='001'){
  cutoff_value <- NULL
} else {
  x <- readLines(paste('PS',job_ps,'_',job_scenario,'_E001_R',job_run,'_',cutoff_prob,'_network_info.csv',sep=''))
  cutoff_value <- x[6]
}

network <- createTemporalNetwork(ps = job_ps,
                                 scenario = job_scenario,
                                 exp = job_exp, 
                                 run = job_run,
                                 cutoff_prob = cutoff_prob, 
                                 cutoff_value = cutoff_value,
                                 layers_to_include = 1:300,
                                 sampled_infections = data$sampled_infections)

# edges <- c(network$intralayer_edges_no_cutoff, network$interlayer_edges_no_cutoff)
# as.tibble(edges) %>% ggplot(aes(value))+geom_density()+geom_vline(xintercept = network$cutoff_value)

sink(paste(base_name,'_network_info.csv',sep=''), append = F)
cat(network$ps);cat('\n')
cat(network$scenario);cat('\n')
cat(network$experiment);cat('\n')
cat(network$run);cat('\n')
cat(network$cutoff_prob);cat('\n')
cat(network$cutoff_value);cat('\n')
sink.reset()
write_csv(network$layer_summary, paste(base_name,'_layer_summary.csv',sep=''))
write_csv(as.tibble(network$intralayer_edges_no_cutoff), paste(base_name,'_intralayer_no_cutoff.csv',sep=''))
write_csv(as.tibble(network$interlayer_edges_no_cutoff), paste(base_name,'_interlayer_no_cutoff.csv',sep=''))

# Infomap objects
infomap <- build_infomap_objects(network_object = network, 
                                 write_to_infomap_file = T, 
                                 infomap_file_name = paste(base_name,'_Infomap_multilayer.txt',sep=''), 
                                 return_objects = T)
write_csv(infomap$nodeList, paste(base_name,'_node_list.csv',sep=''))





