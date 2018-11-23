# Initialize --------------------------------------------------------------
source('functions.R')
library(sqldf, quietly = T)
library(tidyverse, quietly = T)
library(magrittr, quietly = T)
library(igraph, quietly = T)
library(data.table, quietly = T)

if (length(commandArgs(trailingOnly=TRUE))==0) {
  args <- c('03','N','001',0.9, '')
} else {
  args <- commandArgs(trailingOnly=TRUE)
}
job_ps <- as.character(args[1]) 
job_scenario <- as.character(args[2]) 
job_exp <- as.character(args[3]) 
# job_run <- 1
job_run <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
cutoff_prob <- as.numeric(args[4])
task <- as.character(args[5]) # Can be: make_networks or read_infomap_results

base_name <- paste('PS',job_ps,'_',job_scenario,'_E',job_exp,'_R',job_run,'_',cutoff_prob,sep='')

if (task=='make_networks'){
  
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
  # write_csv(as.tibble(network$intralayer_edges_no_cutoff), paste(base_name,'_intralayer_no_cutoff.csv',sep=''))
  # write_csv(as.tibble(network$interlayer_edges_no_cutoff), paste(base_name,'_interlayer_no_cutoff.csv',sep=''))
  
  # Infomap objects
  infomap <- build_infomap_objects(network_object = network, 
                                   write_to_infomap_file = T, 
                                   infomap_file_name = paste(base_name,'_Infomap_multilayer.txt',sep=''), 
                                   return_objects = T)
  write_csv(infomap$nodeList, paste(base_name,'_node_list.csv',sep=''))
} # End task 'make_networks'

if (task=='read_infomap_results'){
  print('Getting results from Infomap...')
  x <- infomap_readTreeFile(PS = job_ps,scenario = job_scenario,exp = job_exp,run=job_run,cutoff_prob = cutoff_prob)
  
  # Find unique repertoires using usearch. This needs to be done because the ABM
  # does not keep track on strain composition to give strains unique names.
  print('Finding unique repertoires using usearch...')
  modules <- x$modules
  sampled_strains <- x$sampled_strains
  sampled_alleles <- x$sampled_alleles
  
  print(paste('Same strains in module and strain data frames?',setequal(sampled_strains$strain_id, modules$strain_id)))
  
  # Produce allele combinations in repertoires
  sampled_alleles$allele_locus <- paste(sampled_alleles$allele, sampled_alleles$locus,sep='_')
  sampled_alleles %<>% select(gene_id, allele_locus) %>% arrange(gene_id,allele_locus)
  sampled_strains %<>% left_join(sampled_alleles) %>% distinct(strain_id,allele_locus)
  
  # Make a file for usearch. Need to use letters because it does not accept 0 and 1
  strain_allele_mat <- xtabs(~strain_id+allele_locus, sampled_strains)
  strain_allele_mat[strain_allele_mat==0] <- 'A'
  strain_allele_mat[strain_allele_mat==1] <- 'G'
  f <- paste('/scratch/midway2/pilosofs/PLOS_Biol/Results/',job_ps,'_',job_scenario,'/',base_name,'_strain_sequences.fasta',sep='')
  sink(f, append = F)
  for (i in 1:nrow(strain_allele_mat)){
    cat('>');cat(rownames(strain_allele_mat)[i]);cat('\n')
    cat(paste(strain_allele_mat[i,],collapse=''));cat('\n')
  }
  sink.reset()
  # rep_clusters_file <- paste('/scratch/midway2/pilosofs/PLOS_Biol/Results/',job_ps,'_',job_scenario,'/',base_name,'_unique_repertoires.txt',sep='')
  rep_clusters_file <- paste(base_name,'_unique_repertoires.txt',sep='')
  system(paste("/scratch/midway2/pilosofs/PLOS_Biol/usearch10.0.240_i86linux32 -fastx_uniques ",f," -relabel rep -tabbedout ",rep_clusters_file,sep='')) # Find unique sequences
  
  print('Getting results from usearch...')
  unique_repertoires <- read_delim(rep_clusters_file,
                                   delim='\t', 
                                   col_names = c('strain_id','strain_cluster','cluster','member_id','repertoires_in_cluster','first_strain_id'))
  
  unique_repertoires$strain_id <- as.character(unique_repertoires$strain_id)
  sampled_strains$strain_id <- as.character(sampled_strains$strain_id)
  modules$strain_id <- as.character(modules$strain_id)
  
  sampled_strains %<>% left_join(unique_repertoires) %>% select(strain_id,allele_locus,strain_cluster)
  modules %<>% left_join(unique_repertoires) %>% 
    mutate(strain_cluster=parse_number(strain_cluster)) %>% 
    select(layer, module, strain_cluster, strain_id, strain_unique, nodeID, PS, scenario, exp, run, cutoff_prob, -first_strain_id, -repertoires_in_cluster, -member_id, -cluster)
    
  
  print('Writing Infomap results with unique repertoires to files...')
  # write_csv(modules, paste('/scratch/midway2/pilosofs/PLOS_Biol/Results/',job_ps,'_',job_scenario,'/',base_name,'_modules.csv',sep=''))
  # write_csv(sampled_strains, paste('/scratch/midway2/pilosofs/PLOS_Biol/Results/',job_ps,'_',job_scenario,'/',base_name,'_sampled_strains.csv',sep=''))
  # write_csv(sampled_alleles, paste('/scratch/midway2/pilosofs/PLOS_Biol/Results/',job_ps,'_',job_scenario,'/',base_name,'_sampled_alleles.csv',sep=''))
  write_csv(modules, paste(base_name,'_modules.csv',sep=''))
  write_csv(sampled_strains, paste(base_name,'_sampled_strains.csv',sep=''))
  write_csv(sampled_alleles, paste(base_name,'_sampled_alleles.csv',sep=''))
  
  print('Done!')
} # End task read_infomap_results

