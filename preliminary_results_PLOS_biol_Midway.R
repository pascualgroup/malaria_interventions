# Functions ---------------------------------------------------------------

source('functions.R')
prep.packages(c('tidyverse','magrittr','sqldf','igraph','data.table','utils'))
scenario_cols <- c('red','orange','blue')
ps_cols <- c('#0A97B7','#B70A97','#97B70A')



get_edge_disributions <- function(PS, scenario, exp, run, cutoff_prob, get_inter=T){
  x <- readLines(paste('Results/',PS,'_',scenario,'/PS',PS,'_',scenario,'_E',exp,'_R',run,'_',cutoff_prob,'_network_info.csv',sep=''))
  cutoff_value <- as.numeric(x[6])
  intra <- fread(paste('Results/',PS,'_',scenario,'/PS',PS,'_',scenario,'_E',exp,'_R',run,'_',cutoff_prob,'_intralayer_no_cutoff.csv',sep=''))  
  if (get_inter){
    inter <- fread(paste('Results/',PS,'_',scenario,'/PS',PS,'_',scenario,'_E',exp,'_R',run,'_',cutoff_prob,'_interlayer_no_cutoff.csv',sep=''))  
    x <- rbind(intra,inter)
  } else {
    x <- intra
  }
  x$PS=PS
  x$scenario=scenario
  x$exp=exp
  x$run=run
  x$cutoff_value=cutoff_value
  x$cutoff_prob=cutoff_prob
  # p <- as.tibble(x) %>% ggplot(aes(value))+geom_density()+geom_vline(xintercept = cutoff_value, color='red')+
  #   labs(title = paste('Cut-off quantile:',cutoff_prob,'Cut-off value:',cutoff_value))
  return(x)
}


get_results_for_cutoff <- function(cutoff_prob_seq=seq(0.05,0.95,0.05), scenario='S', run_range=1:10){
  design_cutoff <- expand.grid(PS=sprintf('%0.2d', 4:6),
                               scenario=scenario, 
                               exp='001',
                               run_range=run_range,
                               cutoff_prob=cutoff_prob_seq,
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
    write_csv(results_cutoff, paste('results_cutoff_',scenario,'_',run,'.csv',sep=''))
  }
  
  results_cutoff <- c()
  for (run in unique(design_cutoff$run_range)){
    file <- paste('results_cutoff_',scenario,'_R',run,'.csv',sep='')
    if(file.exists(file)){
      x <- read_csv(file, col_types = 'iicccicccid')
      results_cutoff <- rbind(results_cutoff, x)  
    } else {
      print(paste('File does not exist: ',file,sep=''))
    }
  }
  results_cutoff <- as.tibble(results_cutoff)
  return(results_cutoff)
}



calculate_module_diversity <- function(PS, scenario, exp, run, cutoff_prob){
  file_modules <- paste('Results/',PS,'_',scenario,'/PS',PS,'_',scenario,'_E',exp,'_R',run,'_',cutoff_prob,'_modules.csv',sep='')
  file_strains <- paste('Results/',PS,'_',scenario,'/PS',PS,'_',scenario,'_E',exp,'_R',run,'_',cutoff_prob,'_sampled_strains.csv',sep='')
  if(file.exists(file_modules) & file.exists(file_strains)){
    print(paste(PS,scenario,exp,run,cutoff_prob,sep=' | '))
    modules <- read_csv(file_modules, col_types = 'iiccciccccd')

    module_persistence <- modules %>% 
      select(scenario, PS, run, cutoff_prob, layer, module) %>% 
      group_by(scenario, PS,run,cutoff_prob,module) %>% 
      summarise(birth_layer=min(layer), death_layer=max(layer), persistence=death_layer-birth_layer+1) %>% 
      mutate(relative_persistence=persistence/(300-birth_layer+1))
  
    sampled_strains <- read_csv(file_strains, col_types = 'ccc')
    sampled_strains <-  sampled_strains[,-3]
    suppressMessages(modules %<>% select(scenario, PS, scenario, exp, run, cutoff_prob, module, strain_id) %>% left_join(sampled_strains))
    allele_freq <- xtabs(~module+allele_locus, modules)
    module_diversity <- vegan::diversity(allele_freq)/log(ncol(allele_freq))
    
   module_persistence$D <- module_diversity
   module_persistence$statistic <- module_diversity*module_persistence$relative_persistence
   return(module_persistence)
  } else {
    print(paste('One file does not exist:',file_modules,file_strains))
  }
}


# Edge weights and cutoffs -----------------------------------------------
# print('Getting edge weight distributions...')
# ## @knitr Edge_weights_distributions_S
# edges_S_04 <- get_edge_disributions(PS = '04',scenario = 'S',exp = '001',1, 0.25)
# edges_S_05 <- get_edge_disributions(PS = '05',scenario = 'S',exp = '001',1, 0.7)
# edges_S_06 <- get_edge_disributions(PS = '06',scenario = 'S',exp = '001',1, 0.9)
# x <- rbind(edges_S_04,edges_S_05,edges_S_06)
# x <- as.tibble(x)
# cutoff_df <- x %>% distinct(PS, cutoff_value, cutoff_prob)
# cutoff_df$label <- paste('PS: ',cutoff_df$PS,'--Cutoff quantile: ',cutoff_df$cutoff_prob,sep='')
# my_labels <- as_labeller(c(`04` = cutoff_df$label[1],
#                            `05` = cutoff_df$label[2],
#                            `06` = cutoff_df$label[3]))
# png('Results/Edge_weights_distributions_S_R1.png', width = 1920, height = 1080)
# x %>% ggplot(aes(value,fill=PS))+
#   geom_density()+
#   geom_vline(data=cutoff_df, aes(xintercept = cutoff_value), color='red', size=2)+
#   scale_fill_manual(values=ps_cols)+
#   facet_wrap(~PS,scales = 'free',labeller = my_labels)+
#   theme_bw(base_size=26)
# dev.off()
# 
# edges_G_04 <- get_edge_disributions(PS = '04',scenario = 'G',exp = '001',1, 0.25)
# edges_G_05 <- get_edge_disributions(PS = '05',scenario = 'G',exp = '001',1, 0.7)
# edges_G_06 <- get_edge_disributions(PS = '06',scenario = 'G',exp = '001',1, 0.9)
# x <- rbind(edges_G_04,edges_G_05,edges_G_06)
# x <- as.tibble(x)
# cutoff_df <- x %>% distinct(PS, cutoff_value, cutoff_prob)
# cutoff_df$label <- paste('PS: ',cutoff_df$PS,'--Cutoff quantile: ',cutoff_df$cutoff_prob,sep='')
# my_labels <- as_labeller(c(`04` = cutoff_df$label[1],
#                            `05` = cutoff_df$label[2],
#                            `06` = cutoff_df$label[3]))
# png('Results/Edge_weights_distributions_G_R1.png', width = 1920, height = 1080)
# x %>% ggplot(aes(value,fill=PS))+
#   geom_density()+
#   geom_vline(data=cutoff_df, aes(xintercept = cutoff_value), color='orange', size=2)+
#   scale_fill_manual(values=ps_cols)+
#   facet_wrap(~PS,scales = 'free',labeller = my_labels)+
#   theme_bw(base_size=26)
# dev.off()
# 
# edges_N_04 <- get_edge_disributions(PS = '04',scenario = 'N',exp = '001',1, 0.25, get_inter = F)
# edges_N_05 <- get_edge_disributions(PS = '05',scenario = 'N',exp = '001',1, 0.7, get_inter = F)
# edges_N_06 <- get_edge_disributions(PS = '06',scenario = 'N',exp = '001',1, 0.9, get_inter = F)
# x <- rbind(edges_N_04,edges_N_05,edges_N_06)
# x <- as.tibble(x)
# cutoff_df <- x %>% distinct(PS, cutoff_value, cutoff_prob)
# cutoff_df$label <- paste('PS: ',cutoff_df$PS,'--Cutoff quantile: ',cutoff_df$cutoff_prob,sep='')
# my_labels <- as_labeller(c(`04` = cutoff_df$label[1],
#                            `05` = cutoff_df$label[2],
#                            `06` = cutoff_df$label[3]))
# png('Results/Edge_weights_distributions_N_R1.png', width = 1920, height = 1080)
# x %>% ggplot(aes(value,fill=PS))+
#   geom_density()+
#   geom_vline(data=cutoff_df, aes(xintercept = cutoff_value), color='blue', size=2)+
#   scale_fill_manual(values=ps_cols)+
#   facet_wrap(~PS,scales = 'free',labeller = my_labels)+
#   theme_bw(base_size=26)
# dev.off()
# edges_N <- x
# 

# Edge cutoff (SELECTION) --------------------------------------------------------
print('Getting edge cutoffs for selection only (Boxplots)')
# Get results
results_cutoff <- get_results_for_cutoff(cutoff_prob_seq = seq(0.05,0.95,0.05), scenario = 'S', run_range = 1:10)
write_csv(results_cutoff, 'Results/results_cutoff_S.csv')

# Examples for modules
png('Results/module_examples.png', width = 1920, height = 1080)
results_cutoff %>% 
  filter(run==2) %>%
  distinct(module,layer, PS, cutoff_prob) %>% 
  # filter(scenario=='S') %>% 
  ggplot(aes(x=layer, y=module, group=PS, color=PS))+
  geom_point(size=1)+
  scale_color_manual(values = ps_cols)+
  scale_x_continuous(breaks = seq(0,300,50))+
  labs(y= 'module ID', x='Time (months)', title='Structure example')+
  facet_grid(cutoff_prob~PS, scales='free')+
  mytheme
dev.off()

# Module and repertoire persistence
module_persistence <- results_cutoff %>% 
  select(PS, run, cutoff_prob, layer, module) %>% 
  group_by(PS,run,cutoff_prob,module) %>% 
  summarise(birth_layer=min(layer), death_layer=max(layer), persistence=death_layer-birth_layer+1) %>% 
  mutate(relative_persistence=persistence/(300-birth_layer+1)) %>% 
  mutate(type='Module') %>% 
  rename(id=module) %>% mutate(id=as.character(id))

strain_persistence <- results_cutoff %>% 
  select(PS, run, cutoff_prob, layer, strain_cluster) %>% 
  group_by(PS,run,cutoff_prob,strain_cluster) %>% 
  summarise(birth_layer=min(layer), death_layer=max(layer), persistence=death_layer-birth_layer+1) %>% 
  mutate(relative_persistence=persistence/(300-birth_layer+1)) %>% 
  mutate(type='Repertoire') %>% 
  rename(id=strain_cluster)

persistence_df <- module_persistence %>% bind_rows(strain_persistence) 

ps_labels=as_labeller(c(`04`='Low',
                      `05`='Medium',
                      `06`='High',
                      `Repertoire`='Repertoire',
                      `Module`='Module'))


png('Results/persistence_boxplots.png', width = 1920, height = 1080)
persistence_df %>% 
  ggplot(aes(x=cutoff_prob, y=persistence, group=cutoff_prob, fill=PS))+
  geom_boxplot(outlier.shape = NA)+
  stat_summary(fun.y=mean, color="red", geom="point", 
               shape=18, size=3,show.legend = FALSE)+
  facet_grid(PS~type, scales='free', labeller = ps_labels)+
  scale_x_continuous(breaks = seq(0.05,0.95,0.1))+
  labs(x='Cut off', y=' Persistence')+
  scale_fill_manual(values=ps_cols)+mytheme
dev.off()

png('Results/relative_persistence_boxplots.png', width = 1920, height = 1080)
persistence_df %>% 
  ggplot(aes(x=cutoff_prob, y=relative_persistence, group=cutoff_prob, fill=PS))+
  geom_boxplot(outlier.shape = NA)+
  stat_summary(fun.y=mean, color="red", geom="point", 
               shape=18, size=3,show.legend = FALSE)+  
  facet_grid(PS~type, scales='free', labeller = ps_labels)+
  scale_x_continuous(breaks = seq(0.05,0.95,0.1))+
  labs(x='Cut off', y=' Relative persistence')+
  scale_fill_manual(values=ps_cols)+mytheme
dev.off()

# Modules per layer
png('Results/modules_per_layer_boxplots.png', width = 1920, height = 1080)
results_cutoff %>% 
  group_by(PS, scenario, run, cutoff_prob, layer) %>% 
  summarise(modules_per_layer=length(unique(module))) %>% 
  ggplot(aes(x=cutoff_prob, y=modules_per_layer, group=cutoff_prob, fill=PS))+
  geom_boxplot(outlier.shape = NA)+
  stat_summary(fun.y=mean, color="red", geom="point", 
               shape=18, size=3,show.legend = FALSE)+
  facet_grid(~PS, scales='free', labeller = ps_labels)+
  scale_x_continuous(breaks = seq(0.05,0.95,0.1))+
  labs(x='Cut off', y=' Modules per layer')+
  scale_fill_manual(values=ps_cols)+mytheme
dev.off()

# Reps per module
png('Results/Repertoires_per_module.png', width = 1920, height = 1080)
results_cutoff %>% 
  group_by(PS, scenario, run, cutoff_prob, module) %>% 
  summarise(repertoires_per_module=length(unique(strain_cluster))) %>% 
  # group_by(PS, scenario, run, cutoff_prob) %>% summarise(repertoires_per_module=mean(repertoires_per_module)) %>%
  ggplot(aes(x=cutoff_prob, y=repertoires_per_module, group=cutoff_prob, fill=PS))+
  geom_boxplot(outlier.size=0)+
  stat_summary(fun.y=mean, color="red", geom="point", 
               shape=18, size=3,show.legend = FALSE)+
  facet_grid(~PS, scales='free', labeller = ps_labels)+
  scale_x_continuous(breaks = seq(0.05,0.95,0.1))+
  labs(x='Cut off', y='Repertoires per module')+
  scale_fill_manual(values=ps_cols)+mytheme
dev.off()






# quantiles_95 <- function(x) {
#   r <- quantile(x, probs=c(0.05, 0.25, 0.5, 0.75, 0.95))
#   names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
#   r
# }
# results_cutoff %>% 
#   group_by(PS, scenario, run, cutoff_prob, module) %>% 
#   summarise(repertoires_per_module=length(unique(strain_cluster))) %>% 
#   # group_by(PS, scenario, run, cutoff_prob) %>% summarise(repertoires_per_module=mean(repertoires_per_module)) %>%
#   ggplot(aes(x=cutoff_prob, y=repertoires_per_module, group=cutoff_prob, fill=PS))+
#   stat_summary(fun.data = quantiles_95, geom="boxplot")+
#   stat_summary(fun.y=mean, color="red", geom="point", 
#                shape=18, size=3,show.legend = FALSE)+
#   facet_grid(~PS, scales='free', labeller = ps_labels)+
#   scale_x_continuous(breaks = seq(0.05,0.95,0.1))+
#   labs(x='Cut off', y='Repertoires per module')+
#   scale_fill_manual(values=ps_cols)+mytheme


# Edge cutoffs (all scenarios) --------------------------------------------
print('Getting edge distributions for all scenarios...')
my_labels <- as_labeller(c(`04` = 'Low',
                           `05` = 'Medium',
                           `06` = 'High',
                           `S` = 'Selection',
                           `G` = 'Generalized immunity',
                           `N` = 'Complete neutrality'))

results_cutoff <- get_results_for_cutoff(cutoff_prob_seq = seq(0.05,0.95,0.05), scenario = 'N', run_range = 1:10)
write_csv(results_cutoff, 'Results/results_cutoff_N.csv')

results_cutoff <- get_results_for_cutoff(cutoff_prob_seq = seq(0.05,0.95,0.05), scenario = 'G', run_range = 1:10)
write_csv(results_cutoff, 'Results/results_cutoff_G.csv')


s <- read_csv('Results/results_cutoff_S.csv', col_types = 'iicccicccid')
n <- read_csv('Results/results_cutoff_N.csv', col_types = 'iicccicccid')
g <- read_csv('Results/results_cutoff_G.csv', col_types = 'iicccicccid')

results_cutoff <- rbind(s,n,g)
results_cutoff %<>% mutate(scenario=factor(scenario, levels=c('S','N','G')))

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
  rename(id=strain_cluster)

persistence_df <- module_persistence %>% bind_rows(strain_persistence) %>% 
  group_by(scenario, PS, cutoff_prob, type) %>% 
  summarise(mean_persistence=mean(persistence),
            mean_relative_persistence=mean(relative_persistence),
            median_persistence=median(persistence),
            median_relative_persistence=median(relative_persistence)
            )
  
png('Results/module_persistence_scenarios.png', width = 1920, height = 1080)
persistence_df %>% 
  filter(type=='Module') %>% 
  ggplot(aes(x=cutoff_prob, color=scenario))+
  geom_point(aes(y=mean_persistence))+geom_line(aes(y=mean_persistence))+
  geom_point(aes(y=median_persistence), shape=15)+geom_line(aes(y=median_persistence), linetype='dashed')+
  facet_grid(PS~scenario, scales='free', labeller = my_labels)+
  scale_x_continuous(breaks = seq(0.05,0.95,0.1))+
  labs(x='Cut off', y=' Mean or median MODULE persistence')+
  scale_color_manual(values=scenario_cols)+mytheme
dev.off()

png('Results/repertoire_persistence_scenarios.png', width = 1920, height = 1080)
persistence_df %>% 
  filter(type=='Repertoire') %>% 
  ggplot(aes(x=cutoff_prob, color=scenario))+
  geom_point(aes(y=mean_persistence))+geom_line(aes(y=mean_persistence))+
  geom_point(aes(y=median_persistence), shape=15)+geom_line(aes(y=median_persistence), linetype='dashed')+
  facet_grid(PS~scenario, scales='free', labeller = my_labels)+
  scale_x_continuous(breaks = seq(0.05,0.95,0.1))+
  labs(x='Cut off', y=' Mean or median REPERTOIRE persistence')+
  scale_color_manual(values=scenario_cols)+mytheme
dev.off()


png('Results/module_relative_persistence_scenarios.png', width = 1920, height = 1080)
persistence_df %>% 
  filter(type=='Module') %>% 
  ggplot(aes(x=cutoff_prob, color=scenario))+
  geom_point(aes(y=mean_relative_persistence))+geom_line(aes(y=mean_relative_persistence))+
  geom_point(aes(y=median_relative_persistence), shape=15)+geom_line(aes(y=median_relative_persistence), linetype='dashed')+
  facet_grid(PS~scenario, scales='free', labeller = my_labels)+
  scale_x_continuous(breaks = seq(0.05,0.95,0.1))+
  labs(x='Cut off', y=' Mean or median MODULE relative persistence')+
  scale_color_manual(values=scenario_cols)+mytheme
dev.off()

png('Results/repertoire_relative_persistence_scenarios.png', width = 1920, height = 1080)
persistence_df %>% 
  filter(type=='Repertoire') %>% 
  ggplot(aes(x=cutoff_prob, color=scenario))+
  geom_point(aes(y=mean_relative_persistence))+geom_line(aes(y=mean_relative_persistence))+
  geom_point(aes(y=median_relative_persistence), shape=15)+geom_line(aes(y=median_relative_persistence), linetype='dashed')+
  facet_grid(PS~scenario, scales='free', labeller = my_labels)+
  scale_x_continuous(breaks = seq(0.05,0.95,0.1))+
  labs(x='Cut off', y=' Mean or median REPERTOIRE relative persistence')+
  scale_color_manual(values=scenario_cols)+mytheme
dev.off()

# Modules per layer
png('Results/modules_per_layer_scenarios.png', width = 1920, height = 1080)
results_cutoff %>% 
  group_by(PS, scenario, run, cutoff_prob, layer) %>% 
  summarise(modules_per_layer=length(unique(module))) %>% 
    group_by(PS, scenario, cutoff_prob) %>%
    summarise(mean_modules_per_layer=mean(modules_per_layer),
              median_modules_per_layer=median(modules_per_layer)) %>%
    ggplot(aes(x=cutoff_prob, color=scenario))+
    geom_point(aes(y=mean_modules_per_layer))+geom_line(aes(y=mean_modules_per_layer))+
    geom_point(aes(y=median_modules_per_layer), shape=15)+geom_line(aes(y=median_modules_per_layer), linetype='dashed')+
    facet_grid(PS~scenario, scales='free', labeller = my_labels)+
    scale_x_continuous(breaks = seq(0.05,0.95,0.1))+
    labs(x='Cut off', y=' Mean or median modules per layer')+
    scale_color_manual(values=scenario_cols)+mytheme
dev.off()

# Reps per module
png('Results/repertoires_per_module_scenarios.png', width = 1920, height = 1080)
results_cutoff %>%
  group_by(PS, scenario, run, cutoff_prob, module) %>%
  summarise(repertoires_per_module=length(unique(strain_cluster))) %>%
  group_by(PS, scenario, cutoff_prob) %>%
  summarise(mean_repertoires_per_module=mean(repertoires_per_module),
            median_repertoires_per_layer=median(repertoires_per_module)) %>%
  ggplot(aes(x=cutoff_prob, color=scenario))+
  geom_point(aes(y=mean_repertoires_per_module))+geom_line(aes(y=mean_repertoires_per_module))+
  geom_point(aes(y=median_repertoires_per_layer), shape=15)+geom_line(aes(y=median_repertoires_per_layer), linetype='dashed')+
  facet_grid(PS~scenario, scales='free', labeller = my_labels)+
  scale_x_continuous(breaks = seq(0.05,0.95,0.1))+
  labs(x='Cut off', y=' Mean or median Repertoires per module')+
  scale_color_manual(values=scenario_cols)+mytheme
dev.off()


# Calculate diversity -----------------------------------------------------
print ('Calculating diversity and the statistic')
design_cutoff <- expand.grid(PS=sprintf('%0.2d', 4:6),
                             scenario=c('S','N','G'), 
                             exp='001',
                             run_range=1:10,
                             cutoff_prob=seq(0.05,0.95,0.05),
                             stringsAsFactors = F)
statistic_results <- c()
for (i in 1:nrow(design_cutoff)){
  PS <- design_cutoff[i,1]
  scenario <- design_cutoff[i,2]
  exp <- design_cutoff[i,3]
  run <- design_cutoff[i,4]
  cutoff_prob <- design_cutoff[i,5]
  x <- calculate_module_diversity(PS, scenario, exp, run, cutoff_prob)
  statistic_results <- rbind(statistic_results, x)
  # print(paste(object.size(statistic_results), nrow(statistic_results)))
}

write_csv(statistic_results, 'Results/statistic_results_S.csv')


png('Results/statistic_results_S.png', width = 1920, height = 1080)
statistic_results %>% 
  ggplot(aes(x=cutoff_prob, y=statistic, group=cutoff_prob, fill=PS))+
  geom_boxplot(outlier.size=0)+
  stat_summary(fun.y=mean, color="red", geom="point", 
               shape=18, size=3,show.legend = FALSE)+
  facet_grid(~PS, scales='free', labeller = ps_labels)+
  scale_x_continuous(breaks = seq(0.05,0.95,0.1))+
  labs(x='Cut off', y='Repertoires per module')+
  scale_fill_manual(values=ps_cols)+mytheme
dev.off()


design_cutoff <- expand.grid(PS=sprintf('%0.2d', 4:6),
                             scenario=c('N'), 
                             exp='001',
                             run_range=1:10,
                             cutoff_prob=seq(0.05,0.95,0.05),
                             stringsAsFactors = F)
statistic_results <- c()
for (i in 1:nrow(design_cutoff)){
  PS <- design_cutoff[i,1]
  scenario <- design_cutoff[i,2]
  exp <- design_cutoff[i,3]
  run <- design_cutoff[i,4]
  cutoff_prob <- design_cutoff[i,5]
  x <- calculate_module_diversity(PS, scenario, exp, run, cutoff_prob)
  statistic_results <- rbind(statistic_results, x)
  # print(paste(object.size(statistic_results), nrow(statistic_results)))
}

write_csv(statistic_results, 'Results/statistic_results_N.csv')

png('statistic_results_N.png', width = 1920, height = 1080)
statistic_results %>% 
  ggplot(aes(x=cutoff_prob, y=statistic, group=cutoff_prob, fill=PS))+
  geom_boxplot(outlier.size=0)+
  stat_summary(fun.y=mean, color="orange", geom="point", 
               shape=18, size=3,show.legend = FALSE)+
  facet_grid(~PS, scales='free', labeller = ps_labels)+
  scale_x_continuous(breaks = seq(0.05,0.95,0.1))+
  labs(x='Cut off', y='Repertoires per module')+
  scale_fill_manual(values=ps_cols)+mytheme
dev.off()



design_cutoff <- expand.grid(PS=sprintf('%0.2d', 4:6),
                             scenario=c('G'), 
                             exp='001',
                             run_range=1:10,
                             cutoff_prob=seq(0.05,0.95,0.05),
                             stringsAsFactors = F)
statistic_results <- c()
for (i in 1:nrow(design_cutoff)){
  PS <- design_cutoff[i,1]
  scenario <- design_cutoff[i,2]
  exp <- design_cutoff[i,3]
  run <- design_cutoff[i,4]
  cutoff_prob <- design_cutoff[i,5]
  x <- calculate_module_diversity(PS, scenario, exp, run, cutoff_prob)
  statistic_results <- rbind(statistic_results, x)
  # print(paste(object.size(statistic_results), nrow(statistic_results)))
}

write_csv(statistic_results, 'Results/statistic_results_G.csv')

png('Results/statistic_results_G.png', width = 1920, height = 1080)
statistic_results %>% 
  ggplot(aes(x=cutoff_prob, y=statistic, group=cutoff_prob, fill=PS))+
  geom_boxplot(outlier.size=0)+
  stat_summary(fun.y=mean, color="blue", geom="point", 
               shape=18, size=3,show.legend = FALSE)+
  facet_grid(~PS, scales='free', labeller = ps_labels)+
  scale_x_continuous(breaks = seq(0.05,0.95,0.1))+
  labs(x='Cut off', y='Repertoires per module')+
  scale_fill_manual(values=ps_cols)+mytheme
dev.off()



s <- read_csv('Results/statistic_results_S.csv', col_types = 'ccidiiiiddd')
n <- read_csv('Results/statistic_results_N.csv', col_types = 'ccidiiiiddd')
g <- read_csv('Results/statistic_results_G.csv', col_types = 'ccidiiiiddd')

statistic_results <- rbind(s,n,g)

x <- statistic_results %>% 
  mutate(scenario=factor(scenario, levels=c('S','N','G'))) %>% 
  group_by(scenario, PS, cutoff_prob) %>% 
  summarise(mean_D=mean(D),
            median_D=median(D),
            sd_D=sd(D),
            mean_statistic=mean(statistic),
            median_statistic=median(statistic),
            sd_statistic=sd(statistic)) 

png('Results/D_scenarios.png', width = 1920, height = 1080)
x %>% 
  ggplot(aes(x=cutoff_prob, color=scenario))+
    geom_point(aes(y=mean_D))+geom_line(aes(y=mean_D))+
    geom_point(aes(y=median_D), shape=15)+geom_line(aes(y=median_D), linetype='dashed')+
    facet_grid(PS~scenario, scales='free', labeller = my_labels)+
    scale_x_continuous(breaks = seq(0.05,0.95,0.1))+
    labs(x='Cut off', y=' Mean or median D')+
    scale_color_manual(values=scenario_cols)+mytheme
dev.off()

png('Results/statistic_scenarios.png', width = 1920, height = 1080)
x %>% 
  ggplot(aes(x=cutoff_prob, color=scenario))+
  geom_point(aes(y=mean_statistic))+geom_line(aes(y=mean_statistic))+
  geom_point(aes(y=median_statistic), shape=15)+geom_line(aes(y=median_statistic), linetype='dashed')+
  facet_grid(PS~scenario, scales='free', labeller = my_labels)+
  scale_x_continuous(breaks = seq(0.05,0.95,0.1))+
  labs(x='Cut off', y=' Mean or median Statistic')+
  scale_color_manual(values=scenario_cols)+mytheme
dev.off()

# line plots --------------------------------------------------------------



# 
# x <- module_persistence %>%
#   bind_rows(strain_persistence) %>%
#   group_by(PS,scenario,cutoff_prob,type) %>%
#   summarise(mean_persistence=mean(persistence),
#             median_persistence=median(persistence),
#             mean_relative_persistence=mean(relative_persistence),
#             median_relative_persistence=median(relative_persistence),
#             sd_persistence=sd(persistence),
#             sd_relative_persistence=sd(relative_persistence))
# # 
# 
# pdf('persistence.pdf', width = 16, height = 10)
# x %>% 
#   filter(type=='module') %>% 
#   filter(scenario=='S') %>% 
#   ggplot(aes(x=cutoff_prob, color=PS))+
#   geom_point(aes(y=mean_persistence))+geom_line(aes(y=mean_persistence))+
#   geom_point(aes(y=median_persistence), shape=15)+geom_line(aes(y=median_persistence), linetype='dashed')+
#   facet_grid(PS~scenario, scales='free', labeller = my_labels)+
#   scale_x_continuous(breaks = seq(0.05,0.95,0.1))+
#   labs(x='Cut off %', y=' Mean or median persistence', title='MODULE persistence')+
#   scale_color_manual(values=scenario_cols)+mytheme
# 
# x %>% 
#   filter(type=='module') %>% 
#   ggplot(aes(x=cutoff_prob, color=scenario))+
#   geom_point(aes(y=mean_relative_persistence))+geom_line(aes(y=mean_relative_persistence))+
#   geom_point(aes(y=median_relative_persistence), shape=15)+geom_line(aes(y=median_relative_persistence), linetype='dashed')+
#   facet_grid(PS~scenario, scales='free', labeller = my_labels)+
#   scale_x_continuous(breaks = seq(0.05,0.95,0.1))+
#   labs(x='Cut off %', y=' Mean or median relative persistence', title='MODULE relative persistence')+
#   scale_color_manual(values=scenario_cols)+mytheme
# 
# x %>% 
#   filter(type=='repertoire') %>% 
#   ggplot(aes(x=cutoff_prob, color=scenario))+
#   geom_point(aes(y=mean_persistence))+geom_line(aes(y=mean_persistence))+
#   geom_point(aes(y=median_persistence), shape=15)+geom_line(aes(y=median_persistence), linetype='dashed')+
#   facet_grid(PS~scenario, scales='free', labeller = my_labels)+
#   scale_x_continuous(breaks = seq(0.05,0.95,0.1))+
#   labs(x='Cut off %', y=' Mean or median persistence', title='REPERTOIRE persistence')+
#   scale_color_manual(values=scenario_cols)+mytheme
# 
# x %>% 
#   filter(type=='repertoire') %>% 
#   ggplot(aes(x=cutoff_prob, color=scenario))+
#   geom_point(aes(y=mean_relative_persistence))+geom_line(aes(y=mean_relative_persistence))+
#   geom_point(aes(y=median_relative_persistence), shape=15)+geom_line(aes(y=median_relative_persistence), linetype='dashed')+
#   facet_grid(PS~scenario, scales='free', labeller = my_labels)+
#   labs(x='Cut off %', y=' Mean or median relative persistence', title='REPERTOIRE relative persistence')+
#   scale_color_manual(values=scenario_cols)+mytheme
# 
# dev.off()

# # Modules per layer
# pdf('modules_per_layer.pdf', width = 16, height = 10)
# results_cutoff %>% 
#   group_by(PS, scenario, cutoff_prob, layer) %>% 
#   summarise(modules_per_layer=length(unique(module))) %>% 
#   group_by(PS, scenario, cutoff_prob) %>% 
#   summarise(mean_modules_per_layer=mean(modules_per_layer),
#             median_modules_per_layer=median(modules_per_layer)) %>% 
#   ggplot(aes(x=cutoff_prob, color=scenario))+
#   geom_point(aes(y=mean_modules_per_layer))+geom_line(aes(y=mean_modules_per_layer))+
#   geom_point(aes(y=median_modules_per_layer), shape=15)+geom_line(aes(y=median_modules_per_layer), linetype='dashed')+
#   facet_grid(PS~scenario, scales='free', labeller = my_labels)+
#   scale_x_continuous(breaks = seq(0.05,0.95,0.1))+
#   labs(x='Cut off %', y=' Mean or median modules per layer', title='Modules per layer')+
#   scale_color_manual(values=scenario_cols)+mytheme
# dev.off()
# 
# # Repertoires per module
# pdf('Repertoires_per_module.pdf', width = 16, height = 10)
# results_cutoff %>% 
#   group_by(PS, scenario, cutoff_prob, module) %>% 
#   summarise(repertoires_per_module=length(unique(strain_cluster))) %>% 
#   group_by(PS, scenario, cutoff_prob) %>% 
#   summarise(mean_repertoires_per_module=mean(repertoires_per_module),
#             median_repertoires_per_layer=median(repertoires_per_module)) %>% 
#   ggplot(aes(x=cutoff_prob, color=scenario))+
#   geom_point(aes(y=mean_repertoires_per_module))+geom_line(aes(y=mean_repertoires_per_module))+
#   geom_point(aes(y=median_repertoires_per_layer), shape=15)+geom_line(aes(y=median_repertoires_per_layer), linetype='dashed')+
#   facet_grid(PS~scenario, scales='free', labeller = my_labels)+
#   scale_x_continuous(breaks = seq(0.05,0.95,0.1))+
#   labs(x='Cut off %', y=' Mean or median modules per layer', title='Repertoires per module')+
#   scale_color_manual(values=scenario_cols)+mytheme
# dev.off()
