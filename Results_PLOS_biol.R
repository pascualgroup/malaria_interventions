library(tidyverse, quietly = T)
library(magrittr, quietly = T)
library(sqldf, quietly = T)
library(igraph, quietly = T)
library(data.table, quietly = T)
library(googlesheets, quietly = T)
library(utils, quietly = T)

# Functions ---------------------------------------------------------------

## @knitr FUNCTIONS
source('~/Documents/malaria_interventions/functions.R')

## @knitr INITIALIZE

# Initialize important variables ------------------------------------------
setwd('/media/Data/PLOS_Biol/')
scenario_cols <- c('red','orange','blue')

my_labels <- as_labeller(c(`04` = 'Low',
                           `05` = 'Medium',
                           `06` = 'High',
                           `S` = 'Selection',
                           `G` = 'Generalized immunity',
                           `N` = 'Complete neutrality'))

# Fig. 2 ------------------------------------------------------------------
# Figure has 3 scenarios in high diversity, showing the following:
# Module examples, Relative persistence, Temporal Diversity, Number of repertoires per module

PS <- '06'; exp <- '001'; cutoff_prob <- 0.85 # These are fixed becaues the figure is for high diversity
experiments <- expand.grid(run=1:10,scenario=c('S','N','G'),stringsAsFactors = F)

get_modularity_results <- function(PS,scenario,run,cutoff_prob,folder='/media/Data/PLOS_Biol/'){
  file <- paste(folder,'Results/',PS,'_',scenario,'/PS',PS,'_',scenario,'_E',exp,'_R',run,'_',cutoff_prob,'_modules.csv',sep='')
  if(file.exists(file)){
    print(paste(PS,scenario,exp,run,cutoff_prob,sep=' | '))
    x <- read_csv(file, col_types = 'iiccciccccd')  
    return(x)
  } else {
    print(paste('File does not exist: ',file,sep=''))
  }
}

module_results <- c()
for (i in 1:nrow(experiments)){
  x <- get_modularity_results(PS,experiments$scenario[i],experiments$run[i],cutoff_prob)
  module_results <- rbind(module_results,x)
}

# Module examples
# Can only be done for 1 run! So selecta nice run :)
module_results %>% 
  filter(run==2) %>%
  distinct(module, layer, PS, scenario) %>% 
  mutate(scenario=factor(scenario, levels=c('S','N','G'))) %>%
  ggplot(aes(x=layer, y=module, color=scenario))+
  geom_point(size=2)+
  scale_color_manual(values = scenario_cols)+
  labs(y= 'module ID', x='Time (months)')+
  facet_wrap(~scenario, scales='free', labeller = my_labels)+mytheme

# Relative persistence
module_persistence <- results_cutoff %>% 
  select(scenario, PS, run, cutoff_prob, layer, module) %>% 
  group_by(scenario, PS,run,cutoff_prob,module) %>% 
  summarise(birth_layer=min(layer), death_layer=max(layer), persistence=death_layer-birth_layer+1) %>% 
  mutate(relative_persistence=persistence/(300-birth_layer+1)) %>% 
  mutate(type='Module') %>% 
  rename(id=module) %>% mutate(id=as.character(id))

strain_persistence <- results_cutoff %>% 
  select(scenario, PS, run, cutoff_prob, layer, strain_cluster) %>% 
  group_by(scenario, PS,run,cutoff_prob,strain_cluster) %>% 
  summarise(birth_layer=min(layer), death_layer=max(layer), persistence=death_layer-birth_layer+1) %>% 
  mutate(relative_persistence=persistence/(300-birth_layer+1)) %>% 
  mutate(type='Repertoire') %>% 
  rename(id=strain_cluster) %>% mutate(id=as.character(id))

persistence_df <- module_persistence %>% bind_rows(strain_persistence) %>% 
  group_by(scenario, PS, cutoff_prob, type) %>% 
  summarise(mean_persistence=mean(persistence),
            mean_relative_persistence=mean(relative_persistence),
            median_persistence=median(persistence),
            median_relative_persistence=median(relative_persistence)
  )
