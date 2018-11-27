# Functions ---------------------------------------------------------------

## @knitr FUNCTIONS
source('functions.R')
prep.packages(c('tidyverse','magrittr','sqldf','igraph','data.table','utils'))
scenario_cols <- c('red','orange','blue')



# Edge cutoff -------------------------------------------------------------

# Get results
design_cutoff <- expand.grid(PS=sprintf('%0.2d', 4:6),
                            scenario=c('S','N','G'), 
                            exp='001',
                            run_range=1:3,
                            cutoff_prob=seq(0.35,0.4,0.05),
                            stringsAsFactors = F)
for (run in unique(design_cutoff$run_range)){
  results_cutoff <- c()
  design_cutoff_run <- subset(design_cutoff, run_range==run)
  print(design_cutoff_run)
  for (i in 1:nrow(design_cutoff_run)){
    PS <- design_cutoff_run[i,1]
    scenario <- design_cutoff_run[i,2]
    exp <- design_cutoff_run[i,3]
    # run <- design_cutoff_run[i,4]
    cutoff_prob <- design_cutoff_run[i,5]
    file <- paste('Results/',PS,'_',scenario,'/PS',PS,'_',scenario,'_E',exp,'_R',run,'_',cutoff_prob,'_modules.csv',sep='')
    if(file.exists(file)){
      print(paste(PS,scenario,exp,run,cutoff_prob,sep=' | '))
      x <- read_csv(file, col_types = 'iiccciccccd')  
      results_cutoff <- rbind(results_cutoff, x)
    } else {
      print(paste('File does not exist: ',file,sep=''))
    }
  }
  write_csv(results_cutoff, paste('results_cutoff_',run,'.csv',sep=''))
}

results_cutoff <- c()
for (run in unique(design_cutoff$run_range)){
  # x <- read_csv(paste('results_cutoff_',run,'.csv',sep=''))
  x <- fread(paste('results_cutoff_',run,'.csv',sep=''), colClasses = 'iiiicicc')
  results_cutoff <- rbind(results_cutoff, x)
}
results_cutoff <- as.tibble(results_cutoff)

results_cutoff %<>% mutate(scenario=factor(scenario, levels=c('S','N','G')))

# Examples for modules
pdf('module_examples.pdf', width = 16, height = 10)
results_cutoff %>% 
  distinct(module,layer,PS, scenario,cutoff_prob) %>% 
  filter(PS=='06') %>% 
  ggplot(aes(x=layer, y=module, color=scenario))+
  geom_point(size=2)+
  scale_color_manual(values = scenario_cols)+
  labs(y= 'module ID', x='Time (months)', title='Structure example')+
  facet_grid(cutoff_prob~scenario, scales='free')+
  mytheme
dev.off()

# Module and repertoire persistence
module_persistence <- results_cutoff %>% 
  mutate(type='module') %>% 
  group_by(PS,scenario,cutoff_prob,type,module) %>% 
  summarise(birth_layer=min(layer), death_layer=max(layer), persistence=death_layer-birth_layer+1) %>% 
  rename(id=module) %>% mutate(id=as.character(id)) %>% 
  mutate(relative_persistence=persistence/(300-birth_layer+1)) 
strain_persistence <- results_cutoff %>% 
  mutate(type='repertoire') %>% 
  group_by(PS,scenario,cutoff_prob,type,strain_cluster) %>% 
  summarise(birth_layer=min(layer), death_layer=max(layer), persistence=death_layer-birth_layer+1) %>% 
  rename(id=strain_cluster) %>% 
  mutate(relative_persistence=persistence/(300-birth_layer+1)) 

x <- module_persistence %>% 
  bind_rows(strain_persistence) %>%
  group_by(PS,scenario,cutoff_prob,type) %>% 
  summarise(mean_persistence=mean(persistence),
            median_persistence=median(persistence),
            mean_relative_persistence=mean(relative_persistence),
            median_relative_persistence=median(relative_persistence))


# Plot mean and median persistence or relative_persistence
my_labels <- as_labeller(c(`04` = 'Low',
                           `05` = 'Medium',
                           `06` = 'High',
                           `S` = 'Selection',
                           `G` = 'Generalized immunity',
                           `N` = 'Complete neutrality'))

pdf('persistence.pdf', width = 16, height = 10)
x %>% 
  filter(type=='module') %>% 
  ggplot(aes(x=cutoff_prob, color=scenario))+
  geom_point(aes(y=mean_persistence))+geom_line(aes(y=mean_persistence))+
  geom_point(aes(y=median_persistence), shape=15)+geom_line(aes(y=median_persistence), linetype='dashed')+
  facet_grid(PS~scenario, scales='free', labeller = my_labels)+
  scale_x_continuous(breaks = seq(0.05,0.95,0.1))+
  labs(x='Cut off %', y=' Mean or median persistence', title='MODULE persistence')+
  scale_color_manual(values=scenario_cols)+mytheme

x %>% 
  filter(type=='module') %>% 
  ggplot(aes(x=cutoff_prob, color=scenario))+
  geom_point(aes(y=mean_relative_persistence))+geom_line(aes(y=mean_relative_persistence))+
  geom_point(aes(y=median_relative_persistence), shape=15)+geom_line(aes(y=median_relative_persistence), linetype='dashed')+
  facet_grid(PS~scenario, scales='free', labeller = my_labels)+
  scale_x_continuous(breaks = seq(0.05,0.95,0.1))+
  labs(x='Cut off %', y=' Mean or median relative persistence', title='MODULE relative persistence')+
  scale_color_manual(values=scenario_cols)+mytheme

x %>% 
  filter(type=='repertoire') %>% 
  ggplot(aes(x=cutoff_prob, color=scenario))+
  geom_point(aes(y=mean_persistence))+geom_line(aes(y=mean_persistence))+
  geom_point(aes(y=median_persistence), shape=15)+geom_line(aes(y=median_persistence), linetype='dashed')+
  facet_grid(PS~scenario, scales='free', labeller = my_labels)+
  scale_x_continuous(breaks = seq(0.05,0.95,0.1))+
  labs(x='Cut off %', y=' Mean or median persistence', title='REPERTOIRE persistence')+
  scale_color_manual(values=scenario_cols)+mytheme

x %>% 
  filter(type=='repertoire') %>% 
  ggplot(aes(x=cutoff_prob, color=scenario))+
  geom_point(aes(y=mean_relative_persistence))+geom_line(aes(y=mean_relative_persistence))+
  geom_point(aes(y=median_relative_persistence), shape=15)+geom_line(aes(y=median_relative_persistence), linetype='dashed')+
  facet_grid(PS~scenario, scales='free', labeller = my_labels)+
  labs(x='Cut off %', y=' Mean or median relative persistence', title='REPERTOIRE relative persistence')+
  scale_color_manual(values=scenario_cols)+mytheme

dev.off()

# Modules per layer
pdf('modules_per_layer.pdf', width = 16, height = 10)
results_cutoff %>% 
  group_by(PS, scenario, cutoff_prob, layer) %>% 
  summarise(modules_per_layer=length(unique(module))) %>% 
  group_by(PS, scenario, cutoff_prob) %>% 
  summarise(mean_modules_per_layer=mean(modules_per_layer),
            median_modules_per_layer=median(modules_per_layer)) %>% 
  ggplot(aes(x=cutoff_prob, color=scenario))+
  geom_point(aes(y=mean_modules_per_layer))+geom_line(aes(y=mean_modules_per_layer))+
  geom_point(aes(y=median_modules_per_layer), shape=15)+geom_line(aes(y=median_modules_per_layer), linetype='dashed')+
  facet_grid(PS~scenario, scales='free', labeller = my_labels)+
  scale_x_continuous(breaks = seq(0.05,0.95,0.1))+
  labs(x='Cut off %', y=' Mean or median modules per layer', title='Modules per layer')+
  scale_color_manual(values=scenario_cols)+mytheme
dev.off()

# Repertoires per module
pdf('Repertoires_per_module.pdf', width = 16, height = 10)
results_cutoff %>% 
  group_by(PS, scenario, cutoff_prob, module) %>% 
  summarise(repertoires_per_module=length(unique(strain_cluster))) %>% 
  group_by(PS, scenario, cutoff_prob) %>% 
  summarise(mean_repertoires_per_module=mean(repertoires_per_module),
            median_repertoires_per_layer=median(repertoires_per_module)) %>% 
  ggplot(aes(x=cutoff_prob, color=scenario))+
  geom_point(aes(y=mean_repertoires_per_module))+geom_line(aes(y=mean_repertoires_per_module))+
  geom_point(aes(y=median_repertoires_per_layer), shape=15)+geom_line(aes(y=median_repertoires_per_layer), linetype='dashed')+
  facet_grid(PS~scenario, scales='free', labeller = my_labels)+
  scale_x_continuous(breaks = seq(0.05,0.95,0.1))+
  labs(x='Cut off %', y=' Mean or median modules per layer', title='Repertoires per module')+
  scale_color_manual(values=scenario_cols)+mytheme
dev.off()
