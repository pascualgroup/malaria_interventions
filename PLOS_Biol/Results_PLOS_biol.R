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

exp <- '001'

scenario_cols <- c('red','orange','blue') # Order is: S, N, G
my_labels <- as_labeller(c(`04` = 'Low',
                           `05` = 'Medium',
                           `06` = 'High',
                           `S` = 'Selection',
                           `G` = 'Generalized immunity',
                           `N` = 'Complete neutrality'))

# Fig. 2 ------------------------------------------------------------------
# Figure has 3 scenarios in high diversity, showing the following:
# Module examples, Relative persistence, Temporal Diversity, mFst.
# Number of repertoires per module is not so useful.

PS <- '06'; cutoff_prob <- 0.85 # These are fixed becaues the figure is for high diversity
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
# Can only be done for 1 run! So select a nice run :)
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
module_persistence <- module_results %>% 
  select(scenario, PS, run, cutoff_prob, layer, module) %>% 
  group_by(scenario, PS,run,cutoff_prob,module) %>% 
  summarise(birth_layer=min(layer), death_layer=max(layer), persistence=death_layer-birth_layer+1) %>% 
  mutate(relative_persistence=persistence/(300-birth_layer+1)) %>% 
  mutate(type='Module') %>% 
  rename(id=module) %>% mutate(id=as.character(id))

strain_persistence <- module_results %>% 
  select(scenario, PS, run, cutoff_prob, layer, strain_cluster) %>% 
  group_by(scenario, PS,run,cutoff_prob,strain_cluster) %>% 
  summarise(birth_layer=min(layer), death_layer=max(layer), persistence=death_layer-birth_layer+1) %>% 
  mutate(relative_persistence=persistence/(300-birth_layer+1)) %>% 
  mutate(type='Repertoire') %>% 
  rename(id=strain_cluster) %>% mutate(id=as.character(id))

persistence_df <- module_persistence %>% bind_rows(strain_persistence)
persistence_df$scenario <- factor(persistence_df$scenario, levels=c('S','N','G'))

persistence_df %>% 
  ggplot()+
  geom_density(data=subset(persistence_df, type=='Module'), aes(relative_persistence, y=..scaled.., fill=scenario))+
  geom_density(data=subset(persistence_df, type=='Repertoire'), aes(relative_persistence, y=..scaled..),fill='gray',alpha=0.4)+
  # geom_rug(aes(x=relative_persistence, y=0, color=scenario), position = position_jitter(height = 0))+
  scale_fill_manual(values = scenario_cols)+
  labs(x='Relative persistence', y='Density (scaled)')+
  facet_grid(~scenario, scales='free_y', labeller = my_labels)+mytheme

persistence_df %>% 
  ggplot()+
  geom_boxplot(aes(x=type, y=relative_persistence, fill=scenario), notch=T)+
  scale_fill_manual(values = scenario_cols)+
  labs(y='Relative persistence')+
  facet_grid(~scenario, scales='free_y', labeller = my_labels)+mytheme

# Statistical analysis for differences in persistence
persistence_df %>% 
  group_by(scenario, type) %>% 
  summarise(mean_relative_persistence=mean(relative_persistence),
            sd_relative_persistence=sd(relative_persistence),
            median_relative_persistence=median(relative_persistence))

for (scen in c('S','G','N')){
  print(scen)
  module_pers <- subset(persistence_df, type=='Module' & scenario==scen)
  rep_pers <- subset(persistence_df, type=='Repertoire' & scenario==scen)
  print(wilcox.test(module_pers$relative_persistence, rep_pers$relative_persistence))
}

# Repertoires per module
# module_results %>% 
#   mutate(scenario=factor(scenario, levels=c('S','N','G'))) %>% 
#   group_by(scenario, run, module) %>% 
#   summarise(repertoires_per_module=length(unique(strain_cluster))) %>% 
#   select(scenario, run, repertoires_per_module) %>% 
#   ggplot(aes(x=scenario, y=repertoires_per_module, fill=scenario))+
#   geom_boxplot()+
#   scale_fill_manual(values=scenario_cols)+mytheme
# 
# x=module_results %>%
#   filter(scenario=='S') %>%
#   group_by(run, module) %>%
#   summarise(repertoires_per_module=length(unique(strain_cluster))) %>%
#   select(run,repertoires_per_module)
# y=module_results %>%
#   filter(scenario=='N') %>%
#   group_by(run, module) %>%
#   summarise(repertoires_per_module=length(unique(strain_cluster))) %>%
#   select(run,repertoires_per_module)
# t.test(x$repertoires_per_module,y$repertoires_per_module)
# wilcox.test(x$repertoires_per_module,y$repertoires_per_module)

# Temporal diversity
module_results_T <- c()
for (i in 1:nrow(experiments)){
  x <- calculate_module_diversity(PS = PS,scenario = experiments$scenario[i],exp = exp, run = experiments$run[i], cutoff_prob = cutoff_prob)
  module_results_T <- rbind(module_results_T, x)
}
module_results_T$scenario <- factor(module_results_T$scenario, levels=c('S','N','G'))


module_results_T %>% 
  ggplot()+
  geom_density(aes(statistic, y=..scaled.., fill=scenario))+
  geom_rug(aes(x=statistic, y=0, color=scenario), position = position_jitter(height = 0))+
  scale_fill_manual(values = scenario_cols)+
  scale_color_manual(values = scenario_cols)+
  labs(x='Temporal diversity', y='Density (scaled)')+
  facet_grid(~scenario, scales='free_y', labeller = my_labels)+mytheme
module_results_T %>% 
  ggplot()+
  geom_boxplot(aes(x=scenario, y=statistic, fill=scenario))+
  scale_fill_manual(values = scenario_cols)+
  labs(y='Temporal diversity')+mytheme

# mFst

