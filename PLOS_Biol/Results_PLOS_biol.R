library(tidyverse, quietly = T)
library(magrittr, quietly = T)
library(sqldf, quietly = T)
library(igraph, quietly = T)
library(data.table, quietly = T)
library(googlesheets, quietly = T)
library(utils, quietly = T)
library(cowplot)
library(grid)
library(gridExtra)

# Functions ---------------------------------------------------------------

source('~/Documents/malaria_interventions/functions.R')

# Initialize important variables ------------------------------------------
setwd('/media/Data/PLOS_Biol/')

exp <- '001'

scenario_cols <- c('red','blue','orange') # Order is: S, N, G
ps_cols <- c('#0A97B7','#B70A97','#97B70A')
gg_labels <- as_labeller(c(`04` = 'Low',
                           `05` = 'Medium',
                           `06` = 'High',
                           `S` = 'Selection',
                           `G` = 'Generalized immunity',
                           `N` = 'Complete neutrality',
                           `Module` = 'Module',
                           `Repertoire` = 'Repertoire'))

all_experiments <- expand.grid(PS=sprintf('%0.2d', c(4:6,18)),
                           scenario=c('S','N','G'), 
                           exp='001',
                           run=1:50,
                           # cutoff_prob=seq(0.3,0.95,0.05),
                           stringsAsFactors = F)
cutoff_df <- tibble(PS=sprintf('%0.2d', c(4:6,18)),cutoff_prob=c(0.3,0.6,0.85,0.85))
all_experiments <- left_join(all_experiments,cutoff_df)


files_list <- list.files('/media/Data/PLOS_Biol/Results/cutoff_to_use')
files <- tibble(file=files_list,
                PS = sapply(str_split(files_list,'_'),function (x) parse_number(x[1])),
                scenario=sapply(str_split(files_list,'_'),function (x) x[2]),
                experiment=sapply(str_split(files_list,'_'),function (x) str_sub(x[3],2,4)),
                run= sapply(str_split(files_list,'_'),function (x) parse_number(x[4])),
                cutoff_prob= sapply(str_split(files_list,'_'),function (x) parse_number(x[5])),
                type=sapply(str_split(files_list,'_'),function (x) paste(x[6:8],collapse='_'))
                )
files$PS <- sprintf('%0.2d', files$PS)
files$type <- str_remove_all(files$type,'_NA')

files %>% filter(PS=='06') %>% group_by(PS,scenario,cutoff_prob) %>% count(type) %>% print(n=20)

# Fig. 2 ------------------------------------------------------------------
# Figure has 3 scenarios in high diversity, showing the following:
# Module examples, Relative persistence, Temporal Diversity, mFst.
# Number of repertoires per module is not so useful.

PS_for_figure <- '18'

experiments <- subset(all_experiments, PS==PS_for_figure & run%in%1:10)

module_results <- c()
for (i in 1:nrow(experiments)){
  PS <- experiments$PS[i]
  scenario <- experiments$scenario[i]
  run <- experiments$run[i]
  cutoff_prob <- experiments$cutoff_prob[i]
  x <- get_modularity_results(PS,scenario,run,cutoff_prob)
  module_results <- rbind(module_results,x)
}
write_csv(module_results, paste('/media/Data/PLOS_Biol/Results/module_results_',PS_for_figure,'.csv',sep=''))
module_results$scenario <- factor(module_results$scenario, levels=c('S','G','N'))
module_results$strain_cluster <- as.integer(module_results$strain_cluster)

# Module examples Can only be plotted for one specific run There may be gaps in
# the layers for a given module because modules are based on sampled
# layer-repertoire tuples, and repertoires may not be sampled in every layer in the ABM.
# So there are two options to plot: one that shows the gaps, and one that does not.
module_results %>% 
  filter(run==1) %>%
  distinct(module, layer, PS, scenario) %>% 
  ggplot(aes(x=layer, y=module, color=scenario))+
  geom_point(size=2)+
  scale_color_manual(values = scenario_cols)+
  labs(y= 'module ID', x='Time (months)')+
  facet_wrap(~scenario, scales='free', labeller = gg_labels)+mytheme

module_results %>% 
  filter(run==1) %>%
  distinct(module, layer, PS, scenario) %>% 
  group_by(scenario,module) %>% summarise(birth_layer=min(layer),death_layer=max(layer)+1) %>% 
  ggplot(aes(xmin=birth_layer, xmax=death_layer, ymin=module, ymax=module, color=scenario))+
  geom_rect(size=2)+
  scale_color_manual(values = scenario_cols)+
  labs(y= 'module ID', x='Time (months)')+
  facet_wrap(~scenario, scales='free', labeller = gg_labels)+manuscript_theme

module_results %>% 
  group_by(PS,scenario,run,layer) %>% 
  summarise(repertoires=length(unique(strain_cluster))) %>% 
  group_by(PS,scenario,layer) %>% 
  summarise(repertoires_mean=mean(repertoires),
            repertoires_sd=sd(repertoires)) %>% 
  ggplot()+
  geom_errorbar(aes(x=layer, ymax=repertoires_mean+repertoires_sd, ymin=repertoires_mean-repertoires_sd, group=scenario),color='gray50')+
  geom_line(aes(x=layer, y=repertoires_mean, color=scenario))+
  scale_color_manual(values = scenario_cols)+
  labs(y= '# repertoires', x='Time (months)')+mytheme

module_results %>% 
  group_by(PS,scenario,run,layer) %>% 
  summarise(modules=length(unique(module))) %>% 
  group_by(PS,scenario,layer) %>% 
  summarise(modules_mean=mean(modules),
            modules_sd=sd(modules)) %>% 
  ggplot()+
  geom_errorbar(aes(x=layer, ymax=modules_mean+modules_sd, ymin=modules_mean-modules_sd, group=scenario),color='gray50')+
  geom_line(aes(x=layer, y=modules_mean, color=scenario))+
  scale_color_manual(values = scenario_cols)+
  labs(y= '# modules', x='Time (months)')+mytheme
    
  

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
persistence_df$scenario <- factor(persistence_df$scenario, levels=c('S','G','N'))
write_csv(persistence_df, paste('/media/Data/PLOS_Biol/Results/persistence_df_',PS_for_figure,'.csv',sep=''))

persistence_df %>% 
  ggplot()+
  geom_density(data=subset(persistence_df, type=='Module'), aes(relative_persistence, fill=scenario))+
  geom_density(data=subset(persistence_df, type=='Repertoire'), aes(relative_persistence),fill='gray',alpha=0.4)+
  # geom_rug(aes(x=relative_persistence, y=0, color=scenario), position = position_jitter(height = 0))+
  scale_fill_manual(values = scenario_cols)+
  labs(x='Relative persistence', y='Density')+
  facet_wrap(~scenario, scales='free_y', labeller = gg_labels)+mytheme

persistence_df %>% 
  ggplot()+
  geom_boxplot(aes(x=type, y=relative_persistence, fill=scenario), notch=F)+
  scale_fill_manual(values = scenario_cols)+
  labs(y='Relative persistence')+
  facet_grid(~scenario, scales='free_y', labeller = gg_labels)+mytheme

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

# Create the figure

module_results <- fread(paste('/media/Data/PLOS_Biol/Results/module_results_',PS_for_figure,'.csv',sep=''), colClasses = list(integer=c(1,2,3,4,6,10),character=c(5,7,8,9),double=11))
module_results <- as.tibble(module_results)
module_results$scenario <- factor(module_results$scenario, levels=c('S','G','N'))
module_results$strain_cluster <- as.integer(module_results$strain_cluster)

persistence_df <- read_csv(paste('/media/Data/PLOS_Biol/Results/persistence_df_',PS_for_figure,'.csv',sep=''))
persistence_df$scenario <- factor(persistence_df$scenario, levels=c('S','G','N'))
persistence_df$id <- as.character(persistence_df$id)


# Temporal diversity
temporal_diversity <- c()
for (i in 1:nrow(experiments)){
  PS <- experiments$PS[i]
  scenario <- experiments$scenario[i]
  run <- experiments$run[i]
  cutoff_prob <- experiments$cutoff_prob[i]
  x <- get_temporal_diversity(PS,scenario,run,cutoff_prob)
  temporal_diversity <- rbind(temporal_diversity,x)
}
temporal_diversity <- as.tibble(temporal_diversity) %>% mutate(scenario=factor(scenario, levels=c('S','G','N')))
# temporal_diversity %>% 
#   ggplot()+
#   geom_density(aes(statistic, fill=scenario))+
#   # geom_rug(aes(x=statistic, y=0, color=scenario), position = position_jitter(height = 0))+
#   scale_fill_manual(values = scenario_cols)+
#   scale_color_manual(values = scenario_cols)+
#   facet_wrap(~scenario, scales='free_y', labeller = gg_labels)+
#   labs(x='Temporal diversity', y='Density')+mytheme
# temporal_diversity %>% 
#   ggplot()+
#   geom_boxplot(aes(x=scenario, y=statistic, fill=scenario))+
#   scale_fill_manual(values = scenario_cols)+
#   labs(y='Temporal diversity')+mytheme

# mFst
mFst <- c()
for (i in 1:nrow(experiments)){
  PS <- experiments$PS[i]
  scenario <- experiments$scenario[i]
  run <- experiments$run[i]
  cutoff_prob <- experiments$cutoff_prob[i]
  x <- get_mFst(PS,scenario,run,cutoff_prob)
  mFst <- rbind(mFst,x)
}
mFst <- as.tibble(mFst) %>% mutate(scenario=factor(scenario, levels=c('S','G','N')))
# mFst %>% 
#   filter(scenario!='N') %>%
#   ggplot()+
#   geom_density(aes(mFst, fill=scenario), alpha=0.3)+
#   scale_fill_manual(values = scenario_cols)+
#   scale_color_manual(values = scenario_cols)+
#   labs(x='mFst', y='Density')+
#   mytheme
# mFst %>% 
#   ggplot()+
#   geom_boxplot(aes(x=scenario, y=mFst, fill=scenario))+
#   scale_fill_manual(values = scenario_cols)+
#   labs(y='mFst')+mytheme


panel_A <- module_results %>% 
  filter(run==1) %>%
  filter(scenario=='S') %>% 
  distinct(module, layer) %>% 
  group_by(module) %>% summarise(birth_layer=min(layer),death_layer=max(layer)+1) %>% 
  ggplot(aes(xmin=birth_layer, xmax=death_layer, ymin=module, ymax=module))+
  geom_rect(size=2, color=scenario_cols[1])+
  manuscript_theme+theme(axis.title = element_blank())
panel_B <- module_results %>% 
  filter(run==1) %>%
  filter(scenario=='G') %>% 
  distinct(module, layer) %>% 
  group_by(module) %>% summarise(birth_layer=min(layer),death_layer=max(layer)+1) %>% 
  ggplot(aes(xmin=birth_layer, xmax=death_layer, ymin=module, ymax=module))+
  geom_rect(size=2, color=scenario_cols[2])+
  manuscript_theme+theme(axis.title = element_blank())
panel_C <- module_results %>% 
  filter(run==1) %>%
  filter(scenario=='N') %>% 
  distinct(module, layer) %>% 
  group_by(module) %>% summarise(birth_layer=min(layer),death_layer=max(layer)+1) %>% 
  ggplot(aes(xmin=birth_layer, xmax=death_layer, ymin=module, ymax=module))+
  geom_rect(size=2, color=scenario_cols[3])+
  manuscript_theme+theme(axis.title = element_blank())
Fig <- plot_grid(panel_A,panel_B,panel_C, labels=c('A','B','C'), ncol=3, align='vh', label_size = 18, scale=0.95)
y.grob <- textGrob("Module ID", gp=gpar(fontface="bold", col="black", fontsize=16), rot=90)
x.grob <- textGrob("Time (months)", gp=gpar(fontface="bold", col="black", fontsize=16), vjust = -0.8)
Fig_2ABC <- grid.arrange(arrangeGrob(Fig, left = y.grob, bottom = x.grob))


panel_D <- persistence_df %>%
  ggplot()+
  geom_density(data=subset(persistence_df, type=='Repertoire' & scenario=='S'), aes(relative_persistence),fill='gray')+
  geom_density(data=subset(persistence_df, type=='Module' & scenario=='S'), aes(relative_persistence), fill=scenario_cols[1])+
  manuscript_theme+theme(axis.title = element_blank())
panel_E <- persistence_df %>%
  ggplot()+
  geom_density(data=subset(persistence_df, type=='Repertoire' & scenario=='G'), aes(relative_persistence),fill='gray')+
  geom_density(data=subset(persistence_df, type=='Module' & scenario=='G'), aes(relative_persistence), fill=scenario_cols[2])+
  manuscript_theme+theme(axis.title = element_blank())
panel_F <- persistence_df %>%
  ggplot()+
  geom_density(data=subset(persistence_df, type=='Repertoire' & scenario=='N'), aes(relative_persistence),fill='gray')+
  geom_density(data=subset(persistence_df, type=='Module' & scenario=='N'), aes(relative_persistence), fill=scenario_cols[3])+
  manuscript_theme+theme(axis.title = element_blank())
Fig <- plot_grid(panel_D,panel_E,panel_F, labels=c('D','E','F'), ncol=3, align='vh', label_size = 18, scale=0.95)
y.grob <- textGrob("Density", gp=gpar(fontface="bold", col="black", fontsize=16), rot=90)
x.grob <- textGrob("Relative persistence", gp=gpar(fontface="bold", col="black", fontsize=16), vjust = -0.8)
Fig_2DEF <- grid.arrange(arrangeGrob(Fig, left = y.grob, bottom = x.grob))


panel_G <- temporal_diversity %>% 
  filter(scenario=='S') %>% 
  ggplot()+
  geom_density(aes(statistic), fill=scenario_cols[1])+
  manuscript_theme+theme(axis.title = element_blank())
panel_H <- temporal_diversity %>% 
  filter(scenario=='G') %>% 
  ggplot()+
  geom_density(aes(statistic), fill=scenario_cols[2])+
  manuscript_theme+theme(axis.title = element_blank())
panel_I <- temporal_diversity %>% 
  filter(scenario=='N') %>% 
  ggplot()+
  geom_density(aes(statistic), fill=scenario_cols[3])+
  manuscript_theme+theme(axis.title = element_blank())
Fig <- plot_grid(panel_G,panel_H,panel_I, labels=c('G','H','I'), ncol=3, align='vh', label_size = 18, scale=0.95)
y.grob <- textGrob("Density", gp=gpar(fontface="bold", col="black", fontsize=16), rot=90)
x.grob <- textGrob("Temporal diversity", gp=gpar(fontface="bold", col="black", fontsize=16), vjust = -0.8)
Fig_2GHI <- grid.arrange(arrangeGrob(Fig, left = y.grob, bottom = x.grob))



panel_J <- mFst  %>% 
  filter(scenario=='S') %>% 
  ggplot()+
  geom_density(aes(mFst), fill=scenario_cols[1])+
  manuscript_theme+theme(axis.title = element_blank())
panel_K <- mFst %>% 
  filter(scenario=='G') %>% 
  ggplot()+
  geom_density(aes(mFst), fill=scenario_cols[2])+
  manuscript_theme+theme(axis.title = element_blank())
panel_L <- mFst %>% 
  filter(scenario=='N') %>% 
  ggplot()+
  geom_density(aes(mFst), fill=scenario_cols[3])+
  manuscript_theme+theme(axis.title = element_blank())
Fig <- plot_grid(panel_J,panel_K,panel_L, labels=c('J','K','L'), ncol=3, align='vh', label_size = 18, scale=0.95)
y.grob <- textGrob("Density", gp=gpar(fontface="bold", col="black", fontsize=16), rot=90)
x.grob <- textGrob("mFst", gp=gpar(fontface="bold", col="black", fontsize=16), vjust = -0.8)
Fig_2JKL <- grid.arrange(arrangeGrob(Fig, left = y.grob, bottom = x.grob))

plot_grid(Fig_2ABC,Fig_2DEF,Fig_2GHI,Fig_2JKL, nrow=4, align='vh')



# Cutoff plots (Selection) ------------------------------------------------

results_cutoff <- fread('/media/Data/PLOS_Biol/Results/results_cutoff_S.csv', 
                        colClasses = list(integer=c(1,2,3,4,6,10),
                                          character=c(5,7,8,9),
                                          double=11))
results_cutoff <- as.tibble(results_cutoff)

# ## Examples for modules
# png('Results/module_examples_S.png', width = 1920, height = 1080)
# results_cutoff %>% 
#   filter(run==2) %>%
#   distinct(module,layer, PS, cutoff_prob) %>% 
#   # filter(scenario=='S') %>% 
#   ggplot(aes(x=layer, y=module, group=PS, color=PS))+
#   geom_point(size=1)+
#   scale_color_manual(values = ps_cols)+
#   scale_x_continuous(breaks = seq(0,300,50))+
#   labs(y= 'module ID', x='Time (months)', title='Structure example')+
#   facet_grid(cutoff_prob~PS, scales='free')+
#   mytheme
# dev.off()

## Modules per layer
png('Results/modules_per_layer_boxplots.png', width = 1920, height = 1080)
results_cutoff %>% 
  group_by(PS, scenario, run, cutoff_prob, layer) %>% 
  summarise(modules_per_layer=length(unique(module))) %>% 
  ggplot(aes(x=cutoff_prob, y=modules_per_layer, group=cutoff_prob, fill=PS))+
  geom_boxplot(outlier.shape = NA)+
  stat_summary(fun.y=mean, color="red", geom="point", 
               shape=18, size=3,show.legend = FALSE)+
  facet_grid(~PS, scales='free', labeller = gg_labels)+
  scale_x_continuous(breaks = seq(0.05,0.95,0.1))+
  labs(x='Cut off', y=' Modules per layer')+
  scale_fill_manual(values=ps_cols)+mytheme
dev.off()

## Reps per module
png('Results/Repertoires_per_module_boxplots.png', width = 1920, height = 1080)
results_cutoff %>% 
  group_by(PS, scenario, run, cutoff_prob, module) %>% 
  summarise(repertoires_per_module=length(unique(strain_cluster))) %>% 
  # group_by(PS, scenario, run, cutoff_prob) %>% summarise(repertoires_per_module=mean(repertoires_per_module)) %>%
  ggplot(aes(x=cutoff_prob, y=repertoires_per_module, group=cutoff_prob, fill=PS))+
  geom_boxplot(outlier.size=0)+
  stat_summary(fun.y=mean, color="red", geom="point", 
               shape=18, size=3,show.legend = FALSE)+
  facet_grid(~PS, scales='free', labeller = gg_labels)+
  scale_x_continuous(breaks = seq(0.05,0.95,0.1))+
  labs(x='Cut off', y='Repertoires per module')+
  scale_fill_manual(values=ps_cols)+mytheme
dev.off()

## Module and repertoire persistence
module_persistence_cutoff <- results_cutoff %>% 
  select(PS, run, cutoff_prob, layer, module) %>% 
  group_by(PS,run,cutoff_prob,module) %>% 
  summarise(birth_layer=min(layer), death_layer=max(layer), persistence=death_layer-birth_layer+1) %>% 
  mutate(relative_persistence=persistence/(300-birth_layer+1)) %>% 
  mutate(type='Module') %>% 
  rename(id=module) %>% mutate(id=as.character(id))

strain_persistence_cutoff <- results_cutoff %>% 
  select(PS, run, cutoff_prob, layer, strain_cluster) %>% 
  group_by(PS,run,cutoff_prob,strain_cluster) %>% 
  summarise(birth_layer=min(layer), death_layer=max(layer), persistence=death_layer-birth_layer+1) %>% 
  mutate(relative_persistence=persistence/(300-birth_layer+1)) %>% 
  mutate(type='Repertoire') %>% 
  rename(id=strain_cluster) %>% mutate(id=as.character(id))

persistence_df_cutoff <- module_persistence_cutoff %>% bind_rows(strain_persistence_cutoff) 
write_csv(persistence_df_cutoff, '/media/Data/PLOS_Biol/Results/persistence_df_cutoff_S.csv')

persistence_df_cutoff <- read_csv('/media/Data/PLOS_Biol/Results/persistence_df_cutoff_S.csv')

# png('Results/persistence_boxplots_S.png', width = 1920, height = 1080)
# persistence_df_cutoff %>% 
#   ggplot(aes(x=cutoff_prob, y=persistence, group=cutoff_prob, fill=PS))+
#   geom_boxplot(outlier.shape = NA)+
#   stat_summary(fun.y=mean, color="red", geom="point", 
#                shape=18, size=3,show.legend = FALSE)+
#   facet_grid(PS~type, scales='free', labeller = gg_labels)+
#   scale_x_continuous(breaks = seq(0.05,0.95,0.1))+
#   labs(x='Cut off', y=' Persistence')+
#   scale_fill_manual(values=ps_cols)+mytheme
# dev.off()

png('Results/relative_persistence_boxplots_S.png', width = 1920, height = 1080)
persistence_df_cutoff %>% 
  ggplot(aes(x=cutoff_prob, y=relative_persistence, group=cutoff_prob, fill=PS))+
  geom_boxplot(outlier.shape = NA)+
  stat_summary(fun.y=mean, color="red", geom="point", 
               shape=18, size=3,show.legend = FALSE)+  
  facet_grid(PS~type, scales='free', labeller = gg_labels)+
  scale_x_continuous(breaks = seq(0.05,0.95,0.1))+
  labs(x='Cut off', y=' Relative persistence')+
  scale_fill_manual(values=ps_cols)+mytheme
dev.off()


# Temporal diveristy
statistic_results <- read_csv('/media/Data/PLOS_Biol/Results/statistic_results_S.csv', col_types = 'ccidiiiiddd')
statistic_results %>% 
  ggplot(aes(x=cutoff_prob, y=statistic, group=cutoff_prob, fill=PS))+
  geom_boxplot()+
  stat_summary(fun.y=mean, color="red", geom="point", 
               shape=18, size=3,show.legend = FALSE)+  
  facet_grid(~PS, scales='free', labeller = gg_labels)+
  scale_x_continuous(breaks = seq(0.05,0.95,0.1))+
  labs(x='Cut off', y='Temporal diversity')+
  scale_fill_manual(values=ps_cols)+mytheme


# Cutoff plots (all scenarios) --------------------------------------------

results_cutoff_N <- read_csv(results_cutoff_N, 'Results/results_cutoff_N.csv')
png('Results/module_examples_N.png', width = 1920, height = 1080)
results_cutoff_N %>% 
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


results_cutoff_G <- read_csv(results_cutoff_G, 'Results/results_cutoff_G.csv')
png('Results/module_examples_G.png', width = 1920, height = 1080)
results_cutoff_G %>% 
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


## Join the scenarios to one data frame
s <- read_csv('Results/results_cutoff_S.csv', col_types = 'iicccicccid')
# s <- data.table::fread('Results/results_cutoff_S_R1.csv')
n <- read_csv('Results/results_cutoff_N.csv', col_types = 'iicccicccid')
# n <- data.table::fread('Results/results_cutoff_N_R1.csv')
g <- read_csv('Results/results_cutoff_G.csv', col_types = 'iicccicccid')
# g <- data.table::fread('Results/results_cutoff_G_R1.csv')

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
  rename(id=strain_cluster) %>% mutate(id=as.character(id))

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
  facet_grid(PS~scenario, scales='free', labeller = gg_labels)+
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
  facet_grid(PS~scenario, scales='free', labeller = gg_labels)+
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
  facet_grid(PS~scenario, scales='free', labeller = gg_labels)+
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
  facet_grid(PS~scenario, scales='free', labeller = gg_labels)+
  scale_x_continuous(breaks = seq(0.05,0.95,0.1))+
  labs(x='Cut off', y=' Mean or median REPERTOIRE relative persistence')+
  scale_color_manual(values=scenario_cols)+mytheme
dev.off()

## Modules per layer
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
  facet_grid(PS~scenario, scales='free', labeller = gg_labels)+
  scale_x_continuous(breaks = seq(0.05,0.95,0.1))+
  labs(x='Cut off', y=' Mean or median modules per layer')+
  scale_color_manual(values=scenario_cols)+mytheme
dev.off()

## Reps per module
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
  facet_grid(PS~scenario, scales='free', labeller = gg_labels)+
  scale_x_continuous(breaks = seq(0.05,0.95,0.1))+
  labs(x='Cut off', y=' Mean or median Repertoires per module')+
  scale_color_manual(values=scenario_cols)+mytheme
dev.off()


# Temporal diveristy
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
            sd_statistic=sd(statistic),
            n=length(statistic)) %>% 
  mutate(ymin_D_sd=mean_D-sd_D, 
         ymax_D_sd=mean_D+sd_D,
         ymin_D_ci=mean_D-qnorm(0.975)*sd_D/sqrt(n),
         ymax_D_ci=mean_D+qnorm(0.975)*sd_D/sqrt(n),
         ymin_statistic_sd=mean_statistic-sd_statistic, 
         ymax_statistic_sd=mean_statistic+sd_statistic,
         ymin_statistic_ci=mean_statistic-qnorm(0.975)*sd_statistic/sqrt(n),
         ymax_statistic_ci=mean_statistic+qnorm(0.975)*sd_statistic/sqrt(n))

png('Results/D_scenarios.png', width = 1920, height = 1080)
x %>% 
  ggplot(aes(x=cutoff_prob, color=scenario))+
  geom_point(aes(y=mean_D))+geom_line(aes(y=mean_D))+
  geom_point(aes(y=median_D), shape=15)+geom_line(aes(y=median_D), linetype='dashed')+
  geom_errorbar(aes(ymin=ymin_D_ci,ymax=ymax_D_ci),alpha=0.6, width=0.01, size=1)+
  facet_grid(PS~scenario, scales='free', labeller = gg_labels)+
  scale_x_continuous(breaks = seq(0.05,0.95,0.1))+
  labs(x='Cut off', y=' Mean or median D')+
  scale_color_manual(values=scenario_cols)+mytheme
dev.off()

png('Results/statistic_scenarios.png', width = 1920, height = 1080)
x %>% 
  ggplot(aes(x=cutoff_prob, color=scenario))+
  geom_point(aes(y=mean_statistic))+geom_line(aes(y=mean_statistic))+
  geom_point(aes(y=median_statistic), shape=15)+geom_line(aes(y=median_statistic), linetype='dashed')+
  geom_errorbar(aes(ymin=ymin_statistic_ci,ymax=ymax_statistic_ci),alpha=0.6, width=0.01, size=1)+
  facet_grid(PS~scenario, scales='free', labeller = gg_labels)+
  scale_x_continuous(breaks = seq(0.05,0.95,0.1))+
  labs(x='Cut off', y=' Mean or median Temporal Diversity')+
  scale_color_manual(values=scenario_cols)+mytheme
dev.off()



# Sensitivity analysis ----------------------------------------------------

sensitivity_experiments <- expand.grid(PS=100:183,
                               scenario='S', 
                               exp='001',
                               run=1,
                               cutoff_prob=0.85,
                               stringsAsFactors = F)

scenario <- 'S'
run <- 1
cutoff_prob <- 0.85
module_results_sensitivity <- c()
for (i in 1:nrow(sensitivity_experiments)){
  PS <- sensitivity_experiments$PS[i]
  x <- get_modularity_results(PS,scenario,run,cutoff_prob, folder = '/media/Data/PLOS_Biol/Results/sensitivity_analysis/')
  module_results_sensitivity <- rbind(module_results_sensitivity,x)
}
write_csv(module_results_sensitivity, '/media/Data/PLOS_Biol/Results/module_results_sensitivity.csv')
module_results_sensitivity$strain_cluster <- as.integer(module_results_sensitivity$strain_cluster)

# Module examples Can only be plotted for one specific run There may be gaps in
# the layers for a given module because modules are based on sampled
# layer-repertoire tuples, and repertoires may not be sampled in every layer in the ABM.
# So there are two options to plot: one that shows the gaps, and one that does not.
module_results_sensitivity %>% 
  filter(PS==150) %>%
  distinct(module, layer, PS, scenario) %>% 
  group_by(scenario,module) %>% summarise(birth_layer=min(layer),death_layer=max(layer)+1) %>% 
  ggplot(aes(xmin=birth_layer, xmax=death_layer, ymin=module, ymax=module, color=scenario))+
  geom_rect(size=2)+
  scale_color_manual(values = scenario_cols)+
  labs(y= 'module ID', x='Time (months)')+
  facet_wrap(~scenario, scales='free', labeller = gg_labels)+manuscript_theme

# Relative persistence
module_persistence <- module_results_sensitivity %>% 
  select(scenario, PS, run, cutoff_prob, layer, module) %>% 
  group_by(scenario, PS,run,cutoff_prob,module) %>% 
  summarise(birth_layer=min(layer), death_layer=max(layer), persistence=death_layer-birth_layer+1) %>% 
  mutate(relative_persistence=persistence/(300-birth_layer+1)) %>% 
  mutate(type='Module') %>% 
  rename(id=module) %>% mutate(id=as.character(id))

strain_persistence <- module_results_sensitivity %>% 
  select(scenario, PS, run, cutoff_prob, layer, strain_cluster) %>% 
  group_by(scenario, PS,run,cutoff_prob,strain_cluster) %>% 
  summarise(birth_layer=min(layer), death_layer=max(layer), persistence=death_layer-birth_layer+1) %>% 
  mutate(relative_persistence=persistence/(300-birth_layer+1)) %>% 
  mutate(type='Repertoire') %>% 
  rename(id=strain_cluster) %>% mutate(id=as.character(id))

persistence_df_sensitivity <- module_persistence %>% bind_rows(strain_persistence)
write_csv(persistence_df_sensitivity, '/media/Data/PLOS_Biol/Results/persistence_df_sensitivity.csv')


# Compare module persistence of sensitivity to observed
x <- persistence_df_sensitivity %>% 
  bind_rows(persistence_df) %>% 
  filter(scenario=='S') %>% 
  mutate(grp=ifelse(PS=='06','Obs','Sens'))
# x %<>% filter(PS %in% c('06', '158'))
x %>% 
  ggplot()+
  geom_density(data=subset(x, PS !='06'), aes(relative_persistence), fill='#1E9B95')+
  geom_density(data=subset(x, PS =='06' & run==1), aes(relative_persistence), fill=scenario_cols[1])+
  facet_grid(~type, scales='free_y')+
  labs(x='Relative persistence', y='Density')+
  mytheme
x %>% 
  ggplot()+
  geom_boxplot(aes(x=grp, y=relative_persistence, fill=grp), notch=F)+
  scale_fill_manual(values = c(scenario_cols[1],'#1E9B95'))+
  labs(y='Relative persistence')+
  facet_grid(~type, scales='free_y', labeller = gg_labels)+
  mytheme

# Temporal diversity
temporal_diversity_sensitivity <- c()
scenario <- 'S'
run <- 1
cutoff_prob <- 0.85
for (i in 1:nrow(sensitivity_experiments)){
  PS <- sensitivity_experiments$PS[i]
  x <- get_temporal_diversity(PS,scenario,run,cutoff_prob,folder = '/media/Data/PLOS_Biol/Results/sensitivity_analysis/')
  temporal_diversity_sensitivity <- rbind(temporal_diversity_sensitivity,x)
}
temporal_diversity_sensitivity <- as.tibble(temporal_diversity_sensitivity)

x <- temporal_diversity_sensitivity %>% 
  bind_rows(temporal_diversity) %>% 
  filter(scenario=='S') %>% 
  mutate(grp=ifelse(PS=='06','Obs','Sens'))
x %>% 
  ggplot()+
  geom_density(data=subset(x, PS !='06'), aes(statistic), fill='#1E9B95', alpha=0.4)+
  geom_density(data=subset(x, PS =='06'), aes(statistic), fill=scenario_cols[1], alpha=0.4)+
  labs(x='Temporal diversity', y='Density')+
  mytheme

# mFst
mFst_sensitivity <- c()
scenario <- 'S'
run <- 1
cutoff_prob <- 0.85
for (i in 1:nrow(sensitivity_experiments)){
  PS <- sensitivity_experiments$PS[i]
  x <- get_mFst(PS,scenario,run,cutoff_prob,folder = '/media/Data/PLOS_Biol/Results/sensitivity_analysis/')
  mFst_sensitivity <- rbind(mFst_sensitivity,x)
}
mFst_sensitivity <- as.tibble(mFst_sensitivity)
x <- mFst_sensitivity %>% 
  bind_rows(mFst) %>% 
  filter(scenario=='S') %>% 
  mutate(grp=ifelse(PS=='06','Obs','Sens'))
x %>% 
  ggplot()+
  geom_density(data=subset(x, PS !='06'), aes(mFst), fill='#1E9B95', alpha=0.4)+
  geom_density(data=subset(x, PS =='06'), aes(mFst), fill=scenario_cols[1], alpha=0.4)+
  labs(x='mFst', y='Density')+
  mytheme
