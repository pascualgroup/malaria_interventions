# Functions ---------------------------------------------------------------

source('~/Documents/malaria_interventions/functions.R')
prep.packages(c('tidyverse','magrittr','sqldf','igraph','data.table','googlesheets','utils','cowplot','grid','gridExtra'), verbose = F)

# Initialize important variables ------------------------------------------
setwd('/media/Data/PLOS_Biol/')

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

files %>% filter(type=='modules.csv') %>% group_by(PS,scenario,cutoff_prob) %>% count(type) %>% print(n=Inf)

# Fig. 2 ------------------------------------------------------------------
# Figure has 3 scenarios in high diversity, showing the following:
# Module examples, Relative persistence, Temporal Diversity, mFst.
# Number of repertoires per module is not so useful.

PS_for_figure <- '06'
module_results <- persistence_df <- temporal_diversity <- mFst <- c()
experiments <- subset(all_experiments, PS==PS_for_figure)
exp <- experiments$exp[1]
for (scenario in unique(experiments$scenario)){
  print(paste(PS_for_figure,scenario,exp, sep=' | '))
  
  f <- paste('/media/Data/PLOS_Biol/Results/PS',PS_for_figure,'_',scenario,'_E',exp,'_module_results.csv',sep='')
  x <- fread(f)
  module_results <- rbind(module_results,x)
  
  f <- paste('/media/Data/PLOS_Biol/Results/PS',PS_for_figure,'_',scenario,'_E',exp,'_persistence_df.csv',sep='')
  x <- fread(f)
  persistence_df <- rbind(persistence_df,x)
  
  x <- read_csv(paste('/media/Data/PLOS_Biol/Results/PS',PS_for_figure,'_',scenario,'_E',exp,'_temporal_diversity.csv',sep=''))
  temporal_diversity <- rbind(temporal_diversity,x)
  
  # x <- read_csv(paste('/media/Data/PLOS_Biol/Results/PS',PS_for_figure,'_',scenario,'_E',exp,'_mFst.csv',sep=''))
  # mFst <- rbind(mFst,x)
}

module_results$scenario <- factor(module_results$scenario, levels=c('S','G','N'))
module_results <- as.tibble(module_results) 
module_results %<>% mutate(PS=str_pad(PS, 2, 'left', '0'),exp=str_pad(exp, 3, 'left', '0'))
persistence_df$scenario <- factor(persistence_df$scenario, levels=c('S','G','N'))
persistence_df <- as.tibble(persistence_df) 
persistence_df %<>% mutate(PS=str_pad(PS, 2, 'left', '0'),exp=str_pad(exp, 3, 'left', '0'))
temporal_diversity$scenario <- factor(temporal_diversity$scenario, levels=c('S','G','N'))
mFst$scenario <- factor(mFst$scenario, levels=c('S','G','N'))


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


title <- ggdraw() + draw_label("High diversity", size=20)
dev.off()
dev.off()
pdf('Results/Fig2_PS06.pdf', 16,12)
plot_grid(title, Fig_2ABC,Fig_2DEF,Fig_2GHI, nrow=4, align='vh', rel_heights = c(0.096,0.32,0.32,0.32))
# plot_grid(Fig_2JKL)
dev.off()


# Box plots
png('Results/Fig2_PS06_relative_persistence.png', 800,400)
persistence_df %>% 
  ggplot(aes(x=type, group=type, y=relative_persistence, fill=scenario))+
  geom_boxplot(alpha=0.8)+
  facet_grid(~scenario)+
  scale_fill_manual(values = scenario_cols)+
  labs(x='', y='Relative persistence', title='High diversity')+
  manuscript_theme
dev.off()

pdf('Results/Fig2_PS06_temporal_diversity.pdf', 16,10)
temporal_diversity %>% 
  ggplot(aes(x=scenario, group=scenario, y=statistic, fill=scenario))+
  geom_boxplot(alpha=0.8)+
  scale_fill_manual(values = scenario_cols)+
  labs(x='', y='Temporal diversity', title='High diversity')+
  manuscript_theme
dev.off()

# Time series (for seasonality) --------------------------------------------

module_time_series <- module_results %>% 
  filter(scenario=='S') %>% 
  group_by(PS,scenario,run,layer) %>% 
  summarise(modules=length(unique(module)),
            repertoires=length(unique(strain_cluster))) %>% 
  group_by(PS,scenario,layer) %>% 
  summarise(modules_mean=mean(modules),
            repertoires_mean=mean(repertoires)) 
module_ts_modules <- ts(module_time_series$modules_mean, start=c(1,1), deltat = 1/12)
module_ts_repertoires <- ts(module_time_series$repertoires_mean, start=c(1,1), deltat = 1/12)
decomposedRes_modules <- stats::decompose(module_ts_modules, type="mult") # use type = "additive" for additive components
decomposedRes_repertoires <- stats::decompose(module_ts_repertoires, type="mult") # use type = "additive" for additive components

dev.off()
pdf('Results/Fig2_PS18_time_series.pdf', 16,12)
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
  labs(y= '# repertoires', x='Time (months)', title='Repertoires time series')+mytheme

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
  labs(y= '# modules', x='Time (months)', title='Modules time series')+mytheme

ccf(module_ts_modules,module_ts_repertoires, main='CCF between # repertoires and # modules')
plot (decomposedRes_repertoires) # see plot below
plot(decomposedRes_modules) # see plot below
dev.off()


# Cutoff plots (Selection) ------------------------------------------------

results_cutoff_S <- fread('/media/Data/PLOS_Biol/Results/results_cutoff_S.csv', 
                        colClasses = list(integer=c(1,2,3,4,6,10),
                                          character=c(5,7,8,9),
                                          double=11))
results_cutoff_S <- as.tibble(results_cutoff_S)

persistence_df_cutoff_S <- read_csv('/media/Data/PLOS_Biol/Results/persistence_df_cutoff_S.csv', col_types = 'cidciiidc')
temporal_diversity_cutoff_S <- read_csv('/media/Data/PLOS_Biol/Results/temporal_diversity_cutoff_S.csv', col_types = 'ccidiiiiddd')

## Modules per layer
panel_A <- results_cutoff_S %>% 
  filter(PS=='04') %>% 
  filter(cutoff_prob>=0.25) %>% 
  group_by(scenario, run, cutoff_prob, layer) %>% 
  summarise(modules_per_layer=length(unique(module))) %>% 
  ggplot(aes(x=cutoff_prob, y=modules_per_layer, group=cutoff_prob))+
  geom_boxplot(outlier.size = 0, fill=ps_cols[1])+
  stat_summary(fun.y=mean, color="red", geom="point",shape=18, size=3,show.legend = FALSE)+
  scale_x_continuous(breaks = seq(0.05,0.95,0.1))+
  manuscript_theme+theme(axis.title = element_blank())
panel_B <- results_cutoff_S %>% 
  filter(PS=='05') %>% 
  filter(cutoff_prob>=0.25) %>% 
  group_by(scenario, run, cutoff_prob, layer) %>% 
  summarise(modules_per_layer=length(unique(module))) %>% 
  ggplot(aes(x=cutoff_prob, y=modules_per_layer, group=cutoff_prob))+
  geom_boxplot(outlier.size = 0, fill=ps_cols[2])+
  stat_summary(fun.y=mean, color="red", geom="point",shape=18, size=3,show.legend = FALSE)+
  scale_x_continuous(breaks = seq(0.05,0.95,0.1))+
  manuscript_theme+theme(axis.title = element_blank())
panel_C <- results_cutoff_S %>% 
  filter(PS=='06') %>% 
  filter(cutoff_prob>=0.25) %>% 
  group_by(scenario, run, cutoff_prob, layer) %>% 
  summarise(modules_per_layer=length(unique(module))) %>% 
  ggplot(aes(x=cutoff_prob, y=modules_per_layer, group=cutoff_prob))+
  geom_boxplot(outlier.size = 0, fill=ps_cols[3])+
  stat_summary(fun.y=mean, color="red", geom="point", shape=18, size=3,show.legend = FALSE)+
  scale_x_continuous(breaks = seq(0.05,0.95,0.1))+
  manuscript_theme+theme(axis.title = element_blank())
Fig <- plot_grid(panel_A,panel_B,panel_C, labels=c('A','B','C'), ncol=3, align='vh', label_size = 18, scale=0.95)
y.grob <- textGrob("Modules per layer", gp=gpar(fontface="bold", col="black", fontsize=16), rot=90)
# x.grob <- textGrob("Cut off", gp=gpar(fontface="bold", col="black", fontsize=16), vjust = -0.8)
Fig_cutoff_ABC <- grid.arrange(arrangeGrob(Fig, left = y.grob))

## Relative persistence
panel_D <- persistence_df_cutoff_S %>% 
  filter(PS=='04', type=='Module') %>% 
  filter(cutoff_prob>=0.25) %>% 
  ggplot(aes(x=cutoff_prob, y=relative_persistence, group=cutoff_prob))+
  geom_boxplot(outlier.size = 0, fill=ps_cols[1])+
  stat_summary(fun.y=mean, color="red", geom="point",shape=18, size=3,show.legend = FALSE)+  
  scale_x_continuous(breaks = seq(0.05,0.95,0.1))+
  manuscript_theme+theme(axis.title = element_blank())
panel_E <- persistence_df_cutoff_S %>% 
  filter(PS=='05', type=='Module') %>% 
  filter(cutoff_prob>=0.25) %>% 
  ggplot(aes(x=cutoff_prob, y=relative_persistence, group=cutoff_prob))+
  geom_boxplot(outlier.size = 0, fill=ps_cols[2])+
  stat_summary(fun.y=mean, color="red", geom="point",shape=18, size=3,show.legend = FALSE)+  
  scale_x_continuous(breaks = seq(0.05,0.95,0.1))+
  manuscript_theme+theme(axis.title = element_blank())
panel_F <- persistence_df_cutoff_S %>% 
  filter(PS=='06', type=='Module') %>% 
  filter(cutoff_prob>=0.25) %>% 
  ggplot(aes(x=cutoff_prob, y=relative_persistence, group=cutoff_prob))+
  geom_boxplot(outlier.size = 0, fill=ps_cols[3])+
  stat_summary(fun.y=mean, color="red", geom="point",shape=18, size=3,show.legend = FALSE)+  
  scale_x_continuous(breaks = seq(0.05,0.95,0.1))+
  manuscript_theme+theme(axis.title = element_blank())
Fig <- plot_grid(panel_D,panel_E,panel_F, labels=c('D','E','F'), ncol=3, align='vh', label_size = 18, scale=0.95)
y.grob <- textGrob("Relative persistence", gp=gpar(fontface="bold", col="black", fontsize=16), rot=90)
# x.grob <- textGrob("Cut off", gp=gpar(fontface="bold", col="black", fontsize=16), vjust = -0.8)
Fig_cutoff_DEF <- grid.arrange(arrangeGrob(Fig, left = y.grob))

## Temporal diveristy
panel_G <- temporal_diversity_cutoff_S %>%
  filter(PS=='04') %>% 
  ggplot(aes(x=cutoff_prob, y=statistic, group=cutoff_prob))+
  geom_boxplot(outlier.size = 0, fill=ps_cols[1])+
  stat_summary(fun.y=mean, color="red", geom="point",shape=18, size=3,show.legend = FALSE)+  
  scale_x_continuous(breaks = seq(0.05,0.95,0.1))+
  manuscript_theme+theme(axis.title = element_blank())
panel_H <- temporal_diversity_cutoff_S %>%
  filter(PS=='05') %>% 
  ggplot(aes(x=cutoff_prob, y=statistic, group=cutoff_prob))+
  geom_boxplot(outlier.size = 0, fill=ps_cols[2])+
  stat_summary(fun.y=mean, color="red", geom="point",shape=18, size=3,show.legend = FALSE)+  
  scale_x_continuous(breaks = seq(0.05,0.95,0.1))+
  manuscript_theme+theme(axis.title = element_blank())
panel_I <- temporal_diversity_cutoff_S %>%
  filter(PS=='06') %>% 
  ggplot(aes(x=cutoff_prob, y=statistic, group=cutoff_prob))+
  geom_boxplot(outlier.size = 0, fill=ps_cols[3])+
  stat_summary(fun.y=mean, color="red", geom="point",shape=18, size=3,show.legend = FALSE)+  
  scale_x_continuous(breaks = seq(0.05,0.95,0.1))+
  manuscript_theme+theme(axis.title = element_blank())
Fig <- plot_grid(panel_G,panel_H,panel_I, labels=c('G','H','I'), ncol=3, align='vh', label_size = 18, scale=0.95)
y.grob <- textGrob("Temporal diversity", gp=gpar(fontface="bold", col="black", fontsize=16), rot=90)
x.grob <- textGrob("Cut off", gp=gpar(fontface="bold", col="black", fontsize=16), vjust = -0.8)
Fig_cutoff_GHI <- grid.arrange(arrangeGrob(Fig, left = y.grob, bottom = x.grob))

title <- ggdraw() + draw_label("Cutoff results", size=20)
pdf('Results/Fig_SI_cutoff_S.pdf', 16,12)
plot_grid(title, Fig_cutoff_ABC,Fig_cutoff_DEF,Fig_cutoff_GHI, nrow=4, align='vh', rel_heights = c(0.096,0.32,0.32,0.32))
# Also plot an example of modules
results_cutoff_S %>% 
  # filter(PS=='06') %>% 
  filter(run==1) %>%
  filter(cutoff_prob>=0.5) %>% 
  distinct(PS,cutoff_prob, module, layer) %>% 
  group_by(PS,cutoff_prob,module) %>% summarise(birth_layer=min(layer),death_layer=max(layer)+1) %>% 
  ggplot(aes(xmin=birth_layer, xmax=death_layer, ymin=module, ymax=module, color=PS, fill=PS))+
  geom_rect(size=1)+
  scale_color_manual(values=ps_cols)+
  scale_fill_manual(values=ps_cols)+
  facet_grid(cutoff_prob~PS, scales = 'free')+
  mytheme+theme(axis.title = element_blank(), legend.position = 'none')
dev.off()


# Cutoff plots (all scenarios) --------------------------------------------
results_cutoff_N <- fread('/media/Data/PLOS_Biol/Results/results_cutoff_N.csv', 
                        colClasses = list(integer=c(1,2,3,4,6,10),
                                          character=c(5,7,8,9),
                                          double=11))
results_cutoff_N <- as.tibble(results_cutoff_N)
persistence_df_cutoff_N <- read_csv('/media/Data/PLOS_Biol/Results/persistence_df_cutoff_N.csv', col_types = 'cidciiidc')
temporal_diversity_cutoff_N <- read_csv('/media/Data/PLOS_Biol/Results/temporal_diversity_cutoff_N.csv', col_types = 'ccidiiiiddd')

results_cutoff_G <- fread('/media/Data/PLOS_Biol/Results/results_cutoff_G.csv', 
                          colClasses = list(integer=c(1,2,3,4,6,10),
                                            character=c(5,7,8,9),
                                            double=11))
results_cutoff_G <- as.tibble(results_cutoff_G)
persistence_df_cutoff_G <- read_csv('/media/Data/PLOS_Biol/Results/persistence_df_cutoff_G.csv', col_types = 'cidciiidc')
temporal_diversity_cutoff_G <- read_csv('/media/Data/PLOS_Biol/Results/temporal_diversity_cutoff_G.csv', col_types = 'ccidiiiiddd')

## Join the scenarios to one data frame
results_cutoff_SGN <- rbind(results_cutoff_S,results_cutoff_G,results_cutoff_N)
results_cutoff_SGN$scenario <- factor(results_cutoff_SGN$scenario, levels=c('S','G','N'))

persistence_df_cutoff_S$scenario <- 'S'
persistence_df_cutoff_G$scenario <- 'G'
persistence_df_cutoff_N$scenario <- 'N'
persistence_df_cutoff_SGN <- rbind(persistence_df_cutoff_S,persistence_df_cutoff_G,persistence_df_cutoff_N)
persistence_df_cutoff_SGN$scenario <- factor(persistence_df_cutoff_SGN$scenario, levels=c('S','G','N'))

temporal_diversity_cutoff_SGN <- rbind(temporal_diversity_cutoff_S,temporal_diversity_cutoff_G,temporal_diversity_cutoff_N)
temporal_diversity_cutoff_SGN$scenario <- factor(temporal_diversity_cutoff_SGN$scenario, levels=c('S','G','N'))

## Plotting

pdf('Results/Fig_SI_cutoff_SNG.pdf', 16,12)

persistence_df_cutoff_SGN %>% 
  filter(type=='Module') %>% 
  group_by(scenario, PS, cutoff_prob) %>% 
  summarise(mean_relative_persistence=mean(relative_persistence),
            sd_relative_persistence=sd(relative_persistence),
            median_relative_persistence=median(relative_persistence)
  ) %>% ggplot(aes(x=cutoff_prob, color=scenario))+
  geom_point(aes(y=mean_relative_persistence))+geom_line(aes(y=mean_relative_persistence))+
  # geom_errorbar(aes(x=cutoff_prob, ymin=mean_relative_persistence-sd_relative_persistence, ymax=mean_relative_persistence+sd_relative_persistence),color='gray')+
  geom_point(aes(y=median_relative_persistence), shape=15)+geom_line(aes(y=median_relative_persistence), linetype='dashed')+
  facet_grid(PS~scenario, scales='free', labeller = gg_labels)+
  scale_x_continuous(breaks = seq(0.05,0.95,0.1))+
  labs(x='Cut off', y=' Mean or median module persistence')+
  scale_color_manual(values=scenario_cols)+mytheme


## Modules per layer
results_cutoff_SGN %>% 
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

# Temporal diversity
x <- temporal_diversity_cutoff_SGN %>% 
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

x %>% 
  ggplot(aes(x=cutoff_prob, color=scenario))+
  geom_point(aes(y=mean_statistic))+geom_line(aes(y=mean_statistic))+
  geom_point(aes(y=median_statistic), shape=15)+geom_line(aes(y=median_statistic), linetype='dashed')+
  # geom_errorbar(aes(ymin=ymin_statistic_ci,ymax=ymax_statistic_ci),alpha=0.6, width=0.01, size=1)+
  facet_grid(PS~scenario, scales='free', labeller = gg_labels)+
  scale_x_continuous(breaks = seq(0.05,0.95,0.1))+
  labs(x='Cut off', y=' Mean or median Temporal Diversity')+
  scale_color_manual(values=scenario_cols)+mytheme


dev.off()



# Sensitivity analysis ----------------------------------------------------

sensitivity_experiments <- expand.grid(PS=as.character(100:183),
                               scenario='S', 
                               run=1,
                               cutoff_prob=0.85,
                               stringsAsFactors = F)

scenario <- 'S'
run <- 1
cutoff_prob <- 0.85
module_results_sensitivity <- c()
exp='001'
for (i in 1:nrow(sensitivity_experiments)){
  PS <- sensitivity_experiments$PS[i]
  x <- get_modularity_results(PS,scenario,run,cutoff_prob, folder = '/media/Data/PLOS_Biol/Results/sensitivity_analysis/')
  module_results_sensitivity <- rbind(module_results_sensitivity,x)
}
module_results_sensitivity$strain_cluster <- as.integer(module_results_sensitivity$strain_cluster)
write_csv(module_results_sensitivity, '/media/Data/PLOS_Biol/Results/module_results_sensitivity.csv')

# # Module examples Can only be plotted for one specific run There may be gaps in
# # the layers for a given module because modules are based on sampled
# # layer-repertoire tuples, and repertoires may not be sampled in every layer in the ABM.
# # So there are two options to plot: one that shows the gaps, and one that does not.
# module_results_sensitivity %>% 
#   filter(PS==150) %>%
#   distinct(module, layer, PS, scenario) %>% 
#   group_by(scenario,module) %>% summarise(birth_layer=min(layer),death_layer=max(layer)+1) %>% 
#   ggplot(aes(xmin=birth_layer, xmax=death_layer, ymin=module, ymax=module, color=scenario))+
#   geom_rect(size=2)+
#   scale_color_manual(values = scenario_cols)+
#   labs(y= 'module ID', x='Time (months)')+
#   facet_wrap(~scenario, scales='free', labeller = gg_labels)+manuscript_theme

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


# Read sensitivity results and plot
module_results_sensitivity <- read_csv('/media/Data/PLOS_Biol/Results/module_results_sensitivity.csv')
module_results_sensitivity$PS <- as.character(module_results_sensitivity$PS)
persistence_df_sensitivity <- read_csv('/media/Data/PLOS_Biol/Results/persistence_df_sensitivity.csv')
persistence_df_sensitivity$PS <- as.character(persistence_df_sensitivity$PS )

persistence_df <- fread('/media/Data/PLOS_Biol/Results/PS06_S_E001_persistence_df.csv') # get the PS06 persistence results
persistence_df <- as.tibble(persistence_df)
persistence_df$PS <- str_pad(persistence_df$PS, width = 2, side = 'left', pad = '0')
persistence_df_sensitivity$id <- as.integer(persistence_df_sensitivity$id)
persistence_df$id <- as.integer(persistence_df$id)

persistence_df_sensitivity %<>% 
  bind_rows(persistence_df) %>%
  filter(scenario=='S') %>% 
  mutate(grp=ifelse(PS=='06','Main result','Sensitivity'))
# Compare module persistence of sensitivity to observed
panel_A <- persistence_df_sensitivity %>% 
  filter(type=='Module') %>% 
  ggplot()+
  geom_boxplot(aes(x=grp, y=relative_persistence, fill=grp), notch = T)+
  # geom_density(data=subset(persistence_df_sensitivity, PS !='06'), aes(relative_persistence), fill='#1E9B95', alpha=0.5)+
  # geom_density(data=subset(persistence_df_sensitivity, PS =='06'), aes(relative_persistence), fill=scenario_cols[1], alpha=0.5)+
  # facet_grid(~type, scales='free_y')+
  labs(x='', y='Relative persistence')+
  manuscript_theme


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

temporal_diversity <- fread('/media/Data/PLOS_Biol/Results/PS06_S_E001_temporal_diversity.csv') # get the PS06 persistence results
temporal_diversity <- as.tibble(temporal_diversity)
temporal_diversity$PS <- str_pad(temporal_diversity$PS, width = 2, side = 'left', pad = '0')

temporal_diversity_sensitivity %<>% 
  bind_rows(temporal_diversity) %>% 
  filter(scenario=='S') %>% 
  mutate(grp=ifelse(PS=='06','Main result','Sensitivity'))
panel_B <- temporal_diversity_sensitivity %>% 
  ggplot()+
  geom_boxplot(aes(x=grp, y=statistic, fill=grp), notch = T)+
  # geom_density(data=subset(temporal_diversity_sensitivity, PS !='06'), aes(statistic), fill='#1E9B95', alpha=0.5)+
  # geom_density(data=subset(temporal_diversity_sensitivity, PS =='06'), aes(statistic), fill=scenario_cols[1], alpha=0.5)+
  labs(x='', y='Temporal diversity')+
  manuscript_theme
# panel_B <- panel_B + draw_label('Main result in red; \nsensitivity simulations in teal', 0.2, 7, hjust = 0, vjust = 0, size = 20)

pdf('Results/Fig_SI_sesitivity_selection.pdf', 16,12)
plot_grid(panel_A,panel_B, labels = c('A','B'))
dev.off()



# 
# # mFst
# mFst_sensitivity <- c()
# scenario <- 'S'
# run <- 1
# cutoff_prob <- 0.85
# for (i in 1:nrow(sensitivity_experiments)){
#   PS <- sensitivity_experiments$PS[i]
#   x <- get_mFst(PS,scenario,run,cutoff_prob,folder = '/media/Data/PLOS_Biol/Results/sensitivity_analysis/')
#   mFst_sensitivity <- rbind(mFst_sensitivity,x)
# }
# mFst_sensitivity <- as.tibble(mFst_sensitivity)
# x <- mFst_sensitivity %>% 
#   bind_rows(mFst) %>% 
#   filter(scenario=='S') %>% 
#   mutate(grp=ifelse(PS=='06','Obs','Sens'))
# x %>% 
#   ggplot()+
#   geom_density(data=subset(x, PS !='06'), aes(mFst), fill='#1E9B95', alpha=0.4)+
#   geom_density(data=subset(x, PS =='06'), aes(mFst), fill=scenario_cols[1], alpha=0.4)+
#   labs(x='mFst', y='Density')+
#   mytheme


# Statistical analysis for differences in persistence ----------
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