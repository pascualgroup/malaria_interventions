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
design <- loadExperiments_GoogleSheets(local = F, workBookName = 'PLOS_Biol_design', sheetID = 2) 
ps_range <- sprintf('%0.2d', 7:9)
exp_range <- sprintf('%0.3d', 1)
run_range <- 1
monitored_variables <- c('prevalence', 'meanMOI','n_circulating_strains', 'n_circulating_genes', 'n_alleles', 'n_total_bites')
ps_cols <- c('#0A97B7','#B70A97','#97B70A')
scenario_cols <- c('red','orange','blue')


compare_ps <- function(ps_range=c('04','05','06'), scenario, exp, run_range, cutoff_prob=c(0.25,0.7,0.9)){
  
  cases <- expand.grid(ps=ps_range, scenario=scenario, exp=exp, run=run_range)
  cases$cutoff_prob <- rep(cutoff_prob, length(run_range))
  ps_comparison <- c()
  for (i in 1:nrow(cases)){
    print(paste('PS: ',cases$ps[i],' | Scenario: ',cases$scenario[i],' | exp: ',cases$exp[i], ' | run: ',cases$run[i],sep=''))
    tmp <- get_data(parameter_space = cases$ps[i], scenario = scenario, experiment = cases$exp[i], run = cases$run[i], cutoff_prob = cases$cutoff_prob[i], use_sqlite = F, tables_to_get = 'summary_general')[[1]]
    ps_comparison <- rbind(ps_comparison, tmp)
  }
    
  time_range <- c(28800,max(ps_comparison$time))
  
  p1 <- ps_comparison %>%
    left_join(subset(design, select=c(PS,BITING_RATE_MEAN,N_GENES_INITIAL)), by='PS') %>% 
    mutate(N_GENES_INITIAL=as.factor(N_GENES_INITIAL)) %>% 
    select(-year, -month, -n_infected) %>% 
    filter(time>time_range[1]&time<time_range[2]) %>%
    gather(variable, value, -pop_id, -time, -exp, -PS, -scenario, -run,-N_GENES_INITIAL, -BITING_RATE_MEAN) %>% 
    group_by(pop_id, time, exp, PS, scenario, N_GENES_INITIAL, BITING_RATE_MEAN, variable) %>%
    summarise(value_mean=mean(value), value_sd=sd(value)) %>% # Need to average across runs
    filter(variable %in% monitored_variables) %>%
    ggplot(aes(x=time, y=value_mean, color=N_GENES_INITIAL))+
    geom_line()+
    scale_color_manual(values=ps_cols)+
    scale_x_continuous(breaks=pretty(x=subset(ps_comparison, time>time_range[1]&time<time_range[2])$time,n=5))+
    facet_wrap(~variable, scales='free')+
    mytheme
  
  # EIR
  p2 <- ps_comparison %>% 
    left_join(design, by='PS') %>% 
    mutate(N_GENES_INITIAL=as.factor(N_GENES_INITIAL)) %>% 
    ggplot(aes(x=month, y=EIR, color=N_GENES_INITIAL))+
    geom_boxplot()+
    stat_summary(aes(group=N_GENES_INITIAL), fun.y=mean, geom="line", size=1)+
    facet_wrap(~N_GENES_INITIAL)+
    scale_color_manual(values=ps_cols)+
    scale_y_continuous(breaks=seq(0,20,2))+
    mytheme
  
  return(list(p1,p2))
}

compare_scenarios <- function(PS, scenarios=c('S','N'), run_range, cutoff_prob){
  cases <- expand.grid(scenario=scenarios, exp=sprintf('%0.3d',1), run=run_range)
  scenario_comparison <- c()
  for (i in 1:nrow(cases)){
    print(paste('Scenario: ',cases$scenario[i],' | exp: ',cases$exp[i], ' | run: ',cases$run[i],sep=''))
    tmp <- get_data(parameter_space = PS, scenario = cases$scenario[i], experiment = cases$exp[i], run = cases$run[i], cutoff_prob = cutoff_prob, use_sqlite = F, tables_to_get = 'summary_general')[[1]]
    scenario_comparison <- rbind(scenario_comparison, tmp)
  }
  
  time_range <- c(28800,max(scenario_comparison$time))
  
  p <- scenario_comparison %>%
    mutate(scenario=factor(scenario, levels=scenarios)) %>% 
    select(-year, -month, -n_infected) %>% 
    filter(time>time_range[1]&time<time_range[2]) %>%
    gather(variable, value, -time, -exp, -PS, -scenario, -run) %>% 
    group_by(time, exp, scenario, variable) %>%
    summarise(value_mean=mean(value), value_sd=sd(value)) %>% # Need to average across runs
    filter(variable %in% c('prevalence','n_circulating_strains', 'n_circulating_genes', 'n_alleles')) %>%
    ggplot(aes(x=time, y=value_mean, color=scenario))+
    geom_line()+
    scale_color_manual(values=scenario_cols)+
    facet_wrap(~variable, scales='free')+
    mytheme
  
  return(p)
}

get_edge_disributions <- function(PS, scenario, exp, run, cutoff_prob, get_inter=T){
  x <- readLines(paste('/media/Data/PLOS_Biol/Results/',PS,'_',scenario,'/PS',PS,'_',scenario,'_E',exp,'_R',run,'_',cutoff_prob,'_network_info.csv',sep=''))
  cutoff_value <- as.numeric(x[6])
  intra <- fread(paste('/media/Data/PLOS_Biol/Results/',PS,'_',scenario,'/PS',PS,'_',scenario,'_E',exp,'_R',run,'_',cutoff_prob,'_intralayer_no_cutoff.csv',sep=''))  
  if (get_inter){
    inter <- fread(paste('/media/Data/PLOS_Biol/Results/',PS,'_',scenario,'/PS',PS,'_',scenario,'_E',exp,'_R',run,'_',cutoff_prob,'_interlayer_no_cutoff.csv',sep=''))  
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
## @knitr END

# Diversity and epi comparisons between parameter spaces within an experiment -------------------

## @knitr basic_variables_S
basic_variable_S <- compare_ps(ps_range=c('04','05','06'), 'S', exp = '001', 1:10, cutoff_prob=c(0.9,0.9,0.9))
basic_variable_S[[1]]

## @knitr basic_variables_N
basic_variable_N <- compare_ps(ps_range=c('04','05','06'), 'N', exp = '001', 1:10, cutoff_prob=c(0.9,0.9,0.9))
basic_variable_N[[1]]

## @knitr basic_variables_G
basic_variable_G <- compare_ps(ps_range=c('04','05','06'), 'G', exp = '001', 1:10, cutoff_prob=c(0.9,0.9,0.9))
basic_variable_G[[1]]

## @knitr EIR_S
basic_variable_S[[2]]

## @knitr EIR_N
basic_variable_N[[2]]

## @knitr EIR_G
basic_variable_G[[2]]
## @knitr END


# Diversity and epi comparisons between scenarios within a parameter space and experiment -------

## @knitr compare_scenarios_01
x <- compare_scenarios(PS = '04', scenarios = c('S','N','G'), run_range = 1:10, cutoff_prob = 0.9)
x
## @knitr compare_scenarios_02
x <- compare_scenarios(PS = '05', scenarios = c('S','N','G'), run_range = 1:10, cutoff_prob = 0.9)
x
## @knitr compare_scenarios_03
x <- compare_scenarios(PS = '06', scenarios = c('S','N','G'), run_range = 1:10, cutoff_prob = 0.9)
x
## @knitr END




# Edge weights and cutoffs -----------------------------------------------

## @knitr Edge_weights_distributions_S
edges_S_01 <- get_edge_disributions(PS = '01',scenario = 'S',exp = '001',1, 0.25)
edges_S_02 <- get_edge_disributions(PS = '02',scenario = 'S',exp = '001',1, 0.7)
edges_S_03 <- get_edge_disributions(PS = '03',scenario = 'S',exp = '001',1, 0.9)
x <- rbind(edges_S_01,edges_S_02,edges_S_03)
x <- as.tibble(x)
cutoff_df <- x %>% distinct(PS, cutoff_value, cutoff_prob)
cutoff_df$label <- paste('PS: ',cutoff_df$PS,'--Cutoff quantile: ',cutoff_df$cutoff_prob,sep='')
my_labels <- as_labeller(c(`01` = cutoff_df$label[1],
                           `02` = cutoff_df$label[2],
                           `03` = cutoff_df$label[3]))
png('Edge_weights_distributions_S.png', width = 1920, height = 1080)
x %>% ggplot(aes(value,fill=PS))+
  geom_density()+
  geom_vline(data=cutoff_df, aes(xintercept = cutoff_value), color='red', size=2)+
  scale_fill_manual(values=ps_cols)+
  facet_wrap(~PS,scales = 'free',labeller = my_labels)+
  theme_bw(base_size=26)
dev.off()
edges_S <- x
## @knitr END



## @knitr Edge_weights_distributions_G
edges_G_01 <- get_edge_disributions(PS = '01',scenario = 'G',exp = '001',1, 0.25)
edges_G_02 <- get_edge_disributions(PS = '02',scenario = 'G',exp = '001',1, 0.7)
edges_G_03 <- get_edge_disributions(PS = '03',scenario = 'G',exp = '001',1, 0.9)
x <- rbind(edges_G_01,edges_G_02,edges_G_03)
x <- as.tibble(x)
cutoff_df <- x %>% distinct(PS, cutoff_value, cutoff_prob)
cutoff_df$label <- paste('PS: ',cutoff_df$PS,'--Cutoff quantile: ',cutoff_df$cutoff_prob,sep='')
my_labels <- as_labeller(c(`01` = cutoff_df$label[1],
                           `02` = cutoff_df$label[2],
                           `03` = cutoff_df$label[3]))
png('Edge_weights_distributions_G.png', width = 1920, height = 1080)
x %>% ggplot(aes(value,fill=PS))+
  geom_density()+
  geom_vline(data=cutoff_df, aes(xintercept = cutoff_value), color='orange', size=2)+
  scale_fill_manual(values=ps_cols)+
  facet_wrap(~PS,scales = 'free',labeller = my_labels)+
  theme_bw(base_size=26)
dev.off()
edges_G <- x
## @knitr END

# NOTE THAT IT IS IMPOSSIBLE TO PLOT EDGE WEIGHTS FOR NEUTRALITY: THERE IS NOT ENOUGH MEMORY TO LOAD ALL THE INTERACTIONS
## @knitr Edge_weights_distributions_N
edges_N_01 <- get_edge_disributions(PS = '01',scenario = 'N',exp = '001',1, 0.25, get_inter = F)
edges_N_02 <- get_edge_disributions(PS = '02',scenario = 'N',exp = '001',1, 0.7, get_inter = F)
# edges_N_03 <- get_edge_disributions(PS = '03',scenario = 'N',exp = '001',1, 0.9, get_inter = F)
x <- rbind(edges_N_01,edges_N_02)
x <- as.tibble(x)
cutoff_df <- x %>% distinct(PS, cutoff_value, cutoff_prob)
cutoff_df$label <- paste('PS: ',cutoff_df$PS,'--Cutoff quantile: ',cutoff_df$cutoff_prob,sep='')
my_labels <- as_labeller(c(`01` = cutoff_df$label[1],
                           `02` = cutoff_df$label[2],
                           `03` = cutoff_df$label[3]))
png('Edge_weights_distributions_N.png', width = 1920, height = 1080)
x %>% ggplot(aes(value,fill=PS))+
  geom_density()+
  geom_vline(data=cutoff_df, aes(xintercept = cutoff_value), color='blue', size=2)+
  scale_fill_manual(values=ps_cols)+
  facet_wrap(~PS,scales = 'free',labeller = my_labels)+
  theme_bw(base_size=26)
dev.off()
edges_N <- x
## @knitr END



# Infomap -----------------------------------------------------------------

## @knitr Infomap_load
design_basic <- expand.grid(PS=sprintf('%0.2d', 4:6),
                            scenario=c('S','N','G'), 
                            exp='001',
                            run_range=1,
                            # cutoff=seq(0.2,0.95,0.05),
                            stringsAsFactors = F)
design_basic$cutoff_prob <- rep(c(0.25,0.7,0.9),length(unique(design_basic$scenario))*length(unique(run_range)))

# PS <- design_basic[i,1]
# scenario <- design_basic[i,2]
# exp <- design_basic[i,3]
# run <- design_basic[i,4]
# cutoff_prob <- design_basic[i,5]

module_results <- c()
for (i in 1:nrow(design_basic)){
  PS <- design_basic[i,1]
  scenario <- design_basic[i,2]
  exp <- design_basic[i,3]
  run <- design_basic[i,4]
  cutoff_prob <- design_basic[i,5]
  # x <- infomap_readTreeFile(PS,scenario,exp,run,cutoff_prob)
  # write_csv(x$modules, paste('/media/Data/PLOS_Biol/Results/',PS,'_',scenario,'/PS',PS,'_',scenario,'_E',exp,'_R',run,'_',cutoff_prob,'_modules.csv',sep=''))
  # write_csv(x$sampled_strains, paste('/media/Data/PLOS_Biol/Results/',PS,'_',scenario,'/PS',PS,'_',scenario,'_E',exp,'_R',run,'_',cutoff_prob,'_sampled_strains.csv',sep=''))
  # write_csv(x$sampled_alleles, paste('/media/Data/PLOS_Biol/Results/',PS,'_',scenario,'/PS',PS,'_',scenario,'_E',exp,'_R',run,'_',cutoff_prob,'_sampled_alleles.csv',sep=''))
  file <- paste('/media/Data/PLOS_Biol/Results/',PS,'_',scenario,'/PS',PS,'_',scenario,'_E',exp,'_R',run,'_',cutoff_prob,'_modules.csv',sep='')
  if(file.exists(file)){
    print(paste(PS,scenario,exp,run,cutoff_prob,sep=' | '))
    x <- read_csv(file, col_types = 'iiccciccccd')  
    module_results <- rbind(module_results, x)
  } else {
    print(paste('File does not exist: ',file,sep=''))
  }
}

# write_csv(module_results,'~/Dropbox/Qixin_Shai_Malaria/PLOS_Biol/modularity_results.csv')

my_labels <- as_labeller(c(`04` = 'Low',
                           `05` = 'Medium',
                           `06` = 'High',
                           `S` = 'Selection',
                           `G` = 'Generalized immunity',
                           `N` = 'Complete neutrality'))

## @knitr Infomap_module_example
module_results %>% 
  distinct(module,layer,PS, scenario,cutoff_prob) %>% 
  mutate(scenario=factor(scenario, levels=c('S','N','G'))) %>%
  ggplot(aes(x=layer, y=module, color=scenario))+
  geom_point(size=2)+
  scale_color_manual(values = scenario_cols)+
  labs(y= 'module ID', x='Time (months)')+
  facet_grid(PS~scenario, scales='free', labeller = my_labels)+mytheme


## @knitr Infomap_relative_persistence

# Module and repertoire persistence
module_persistence <- module_results %>% 
  mutate(scenario=factor(scenario, levels=c('S','N','G'))) %>% 
  mutate(type='module') %>% 
  group_by(PS,scenario,type,module) %>% 
  summarise(birth_layer=min(layer), death_layer=max(layer), persistence=death_layer-birth_layer+1) %>% 
  rename(id=module) %>% mutate(id=as.character(id))
strain_persistence <- module_results %>% 
  mutate(scenario=factor(scenario, levels=c('S','N','G'))) %>% 
  mutate(type='repertoire') %>% 
  group_by(PS,scenario,type,strain_cluster) %>% 
  summarise(birth_layer=min(layer), death_layer=max(layer), persistence=death_layer-birth_layer+1) %>% 
  rename(id=strain_cluster)

x <- module_persistence %>% 
  bind_rows(strain_persistence) %>%
  # filter(PS %in% c('01','02','03')) %>% 
  mutate(relative_persistence=persistence/(300-birth_layer+1)) 
x %>% 
  ggplot()+
  geom_density(data=subset(x, type=='module'), aes(relative_persistence, y=..scaled.., fill=scenario))+
  geom_density(data=subset(x, type=='repertoire'), aes(relative_persistence, y=..scaled..),fill='gray',alpha=0.4)+
  # geom_rug(aes(x=relative_persistence, y=0), position = position_jitter(height = 0))+
  scale_fill_manual(values = scenario_cols)+
  labs(x='Relative persistence', y='Density (scaled)')+
  facet_grid(PS~scenario, scales='free_y', labeller = my_labels)+mytheme

## @knitr Infomap_persistence
x %>% 
  ggplot()+
  geom_density(data=subset(x, type=='module'), aes(persistence, y=..scaled.., fill=scenario))+
  geom_density(data=subset(x, type=='repertoire'), aes(persistence, y=..scaled..),fill='gray',alpha=0.4)+
  # geom_rug(aes(x=persistence, y=0), position = position_jitter(height = 0))+
  scale_fill_manual(values = scenario_cols)+
  labs(x='Absolute persistence (months)', y='Density (scaled)')+
  facet_grid(PS~scenario, scales='free_y', labeller = my_labels)+mytheme

## @knitr END

# Modules per layer
module_results %>% 
  mutate(scenario=factor(scenario, levels=c('S','N','G'))) %>% 
  group_by(PS, scenario, layer) %>% 
  summarise(modules_per_layer=length(unique(module))) %>% 
  ggplot(aes(x=layer, y=modules_per_layer, color=scenario))+
  geom_point()+geom_line()+
  facet_wrap(~PS, scales='free', labeller = my_labels)+
  scale_color_manual(values=scenario_cols)+mytheme

# Repertoires per module
module_results %>% 
  mutate(scenario=factor(scenario, levels=c('S','N','G'))) %>% 
  group_by(PS, scenario, module) %>% 
  summarise(repertoires_per_module=length(strain_cluster)) %>% 
  # filter(module==1) %>% 
  ggplot(aes(repertoires_per_module, fill=scenario))+
  geom_density(alpha=0.5)+
  facet_wrap(~PS, scales='free', labeller = my_labels)+
  scale_fill_manual(values=scenario_cols)+mytheme
  


# Edge cutoff -------------------------------------------------------------
#  !!! See updated code in preliminary_results_PLOS_biol_Midway

# Within-module diversity --------------------------------------------------------
#  !!! See updated code in preliminary_results_PLOS_biol_Midway

# Between-module diversity --------------------------------------------------------
#  !!! See updated code in preliminary_results_PLOS_biol_Midway


# Seasonality -------------------------------------------------------------

cases <- expand.grid(ps=c('11','12','13'), scenario='S', exp='001', run=1)
cases$cutoff_prob <- 0.85
ps_comparison <- c()
for (i in 1:nrow(cases)){
  print(paste('PS: ',cases$ps[i],' | Scenario: ',cases$scenario[i],' | exp: ',cases$exp[i], ' | run: ',cases$run[i],sep=''))
  tmp <- get_data(parameter_space = cases$ps[i], scenario = cases$scenario[i], experiment = cases$exp[i], run = cases$run[i], cutoff_prob = cases$cutoff_prob[i], use_sqlite = T, tables_to_get = 'summary_general')[[1]]
  ps_comparison <- rbind(ps_comparison, tmp)
}

time_range <- c(28800,max(ps_comparison$time))

ps_comparison %>%
  select(-year, -month, -n_infected) %>% 
  filter(time>time_range[1]&time<time_range[2]) %>%
  gather(variable, value, -pop_id, -time, -exp, -PS, -scenario, -run) %>% 
  group_by(pop_id, time, exp, PS, scenario, variable) %>%
  summarise(value_mean=mean(value), value_sd=sd(value)) %>% # Need to average across runs
  filter(variable %in% monitored_variables) %>%
  ggplot(aes(x=time, y=value_mean, color=PS))+
  geom_line()+
  scale_x_continuous(breaks=pretty(x=subset(ps_comparison, time>time_range[1]&time<time_range[2])$time,n=5))+
  facet_wrap(~variable, scales='free')+
  mytheme

# EIR
ps_comparison %>% 
  ggplot(aes(x=month, y=EIR, color=PS))+
  geom_boxplot()+
  stat_summary(aes(group=PS), fun.y=mean, geom="line", size=1)+
  facet_wrap(~PS)+
  scale_y_continuous(breaks=seq(0,20,2))+
  mytheme




# This repeats the diversity and epi comparisons, but with seasonality

## @knitr basic_variables_S_seasonality
basic_variable_S <- compare_ps(ps_range=c('09','10'), 'S', exp = '001', 1:10, cutoff_prob=c(0,0.85))
basic_variable_S[[1]]

## @knitr basic_variables_G_seasonality
basic_variable_G <- compare_ps(ps_range=c('07','08','09'), 'G', exp = '001', 1:10, cutoff_prob=c(0.3,0.6,0.85))
basic_variable_G[[1]]

## @knitr basic_variables_N_seasonality
basic_variable_N <- compare_ps(ps_range=c('07','08','09'), 'N', exp = '001', 1:10, cutoff_prob=c(0.3,0.6,0.85))
basic_variable_N[[1]]

## @knitr EIR_S_seasonality
basic_variable_S[[2]]

## @knitr EIR_G_seasonality
basic_variable_G[[2]]
## @knitr END

## @knitr EIR_N_seasonality
basic_variable_N[[2]]


## @knitr compare_scenarios_01_seasonality
x <- compare_scenarios(PS = '04', scenarios = c('S','N','G'), run_range = 1:10, cutoff_prob = 0.9)
x
## @knitr compare_scenarios_02_seasonality
x <- compare_scenarios(PS = '05', scenarios = c('S','N','G'), run_range = 1:10, cutoff_prob = 0.9)
x
## @knitr compare_scenarios_03_seasonality
x <- compare_scenarios(PS = '06', scenarios = c('S','N','G'), run_range = 1:10, cutoff_prob = 0.9)
x
## @knitr END