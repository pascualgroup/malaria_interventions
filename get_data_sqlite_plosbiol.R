library(tidyverse)
library(magrittr)
library(sqldf)
library(igraph)
library(data.table)
library(googlesheets)
library(utils)

# Functions ---------------------------------------------------------------

## @knitr FUNCTIONS
source('~/Documents/malaria_interventions/functions.R')

## @knitr INITIALIZE

# Initialize important variables ------------------------------------------
setwd('/media/Data/PLOS_Biol/')
design <- loadExperiments_GoogleSheets(local = F, workBookName = 'PLOS_Biol_design', sheetID = 2) 
ps_range <- sprintf('%0.2d', 1:3)
exp_range <- sprintf('%0.3d', 1)
run_range <- 1
monitored_variables <- c('prevalence', 'meanMOI','n_circulating_strains', 'n_circulating_genes', 'n_alleles', 'n_total_bites')
ps_cols <- c('#0A97B7','#B70A97','#97B70A')
exp_labels=c('Control','2-yr','5-yr','10-yr')
scenario_cols <- c('red','blue','orange')


compare_variables_ps <- function(scenario, exp, run_range){
  
  ps_comparison <- map(run_range, function(r){
    map(ps_range, function(ps){
      # print(paste(ps,r,sep=' | '))
      tmp <- get_data(parameter_space = ps, scenario = scenario, experiment = exp, run = r, use_sqlite = F, tables_to_get = 'summary_general')
      return(tmp[[1]])
    }) %>% bind_rows()
  }) %>% bind_rows()
  
  time_range <- c(28800,max(ps_comparison$time))
  
  ## @knitr COMPARE_DIVERSITY_PLOT
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
  
  ps_comparison_eir <- map(run_range, function(r){
    map(sprintf('%0.2d', 1:3), function(ps){
      tmp <- get_data(parameter_space = ps, scenario = scenario, experiment = '001', run = r, use_sqlite = F, tables_to_get = 'summary_general')
      return(tmp[[1]])
    }) %>% bind_rows()
  }) %>% bind_rows()
  
  p2 <- ps_comparison_eir %>% 
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

compare_scenarios <- function(PS, scenarios=c('S','N'), run_range){
  cases <- expand.grid(scenario=scenarios, exp=sprintf('%0.3d',1), run=run_range)
  scenario_comparison <- c()
  for (i in 1:nrow(cases)){
    print(paste('Scenario: ',cases$scenario[i],' | exp: ',cases$exp[i], ' | run: ',cases$run[i],sep=''))
    tmp <- get_data(parameter_space = PS, scenario = cases$scenario[i], experiment = cases$exp[i], run = cases$run[i], use_sqlite = F, tables_to_get = 'summary_general')[[1]]
    scenario_comparison <- rbind(scenario_comparison, tmp)
  }
  
  time_range <- c(28800,max(scenario_comparison$time))
  
  p <- scenario_comparison %>%
    mutate(scenario=factor(scenario, levels=c('S','N'))) %>% 
    select(-year, -month, -n_infected) %>% 
    filter(time>time_range[1]&time<time_range[2]) %>%
    gather(variable, value, -time, -exp, -PS, -scenario, -run) %>% 
    group_by(time, exp, scenario, variable) %>%
    summarise(value_mean=mean(value), value_sd=sd(value)) %>% # Need to average across runs
    filter(variable %in% c('prevalence','n_circulating_strains', 'n_circulating_genes', 'n_alleles')) %>%
    ggplot(aes(x=time, y=value_mean, color=scenario))+
    geom_line()+
    scale_color_manual(values=scenario_cols)+
    facet_wrap(~variable+scenario, scales='free')+
    mytheme
  
  return(p)
}

get_edge_disributions <- function(PS, scenario, exp, run, cutoff_prob){
  x <- readLines(paste('/media/Data/PLOS_Biol/Results/',PS,'_',scenario,'/PS',PS,'_',scenario,'_E',exp,'_R',run,'_',cutoff_prob,'_network_info.csv',sep=''))
  cutoff_value <- as.numeric(x[6])
  intra <- fread(paste('/media/Data/PLOS_Biol/Results/',PS,'_',scenario,'/PS',PS,'_',scenario,'_E',exp,'_R',run,'_',cutoff_prob,'_intralayer_no_cutoff.csv',sep=''))  
  inter <- fread(paste('/media/Data/PLOS_Biol/Results/',PS,'_',scenario,'/PS',PS,'_',scenario,'_E',exp,'_R',run,'_',cutoff_prob,'_interlayer_no_cutoff.csv',sep=''))  
  x <- rbind(intra,inter)
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

# Compare between parameter spaces within an experiment -------------------

## @knitr basic_variables_selection
basic_variable_s <- compare_variables_ps('S', exp = '001', 1)
basic_variable_s[[1]]

## @knitr basic_variables_neutral
basic_variable_n <- compare_variables_ps('N', exp = '001', 1)
basic_variable_n[[1]]

## @knitr EIR_selection
basic_variable_s[[2]]

## @knitr EIR_neutral
basic_variable_n[[2]]

## @knitr END


# Compare between scenarios within a parameter space and experiment -------

## @knitr compare_scenarios_01
x <- compare_scenarios(PS = '01', scenarios = c('S','N'), run_range = 1)
x
## @knitr compare_scenarios_02
x <- compare_scenarios(PS = '02', scenarios = c('S','N'), run_range = 1)
x
## @knitr compare_scenarios_03
x <- compare_scenarios(PS = '03', scenarios = c('S','N'), run_range = 1)
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
x %>% ggplot(aes(value,fill=PS))+
  geom_density()+
  geom_vline(data=cutoff_df, aes(xintercept = cutoff_value), color='red')+
  scale_fill_manual(values=ps_cols)+
  facet_wrap(~PS,scales = 'free',labeller = my_labels)
## @knitr END

## @knitr Edge_weights_distributions_N
edges_N_01 <- get_edge_disributions(PS = '01',scenario = 'N',exp = '001',1, 0.25)
edges_N_02 <- get_edge_disributions(PS = '02',scenario = 'N',exp = '001',1, 0.7)
edges_N_03 <- get_edge_disributions(PS = '03',scenario = 'N',exp = '001',1, 0.9)
x <- rbind(edges_N_01,edges_N_02,edges_N_03)
x <- as.tibble(x)
cutoff_df <- x %>% distinct(PS, cutoff_value, cutoff_prob)
cutoff_df$label <- paste('PS: ',cutoff_df$PS,'--Cutoff quantile: ',cutoff_df$cutoff_prob,sep='')
my_labels <- as_labeller(c(`01` = cutoff_df$label[1],
                           `02` = cutoff_df$label[2],
                           `03` = cutoff_df$label[3]))
x %>% ggplot(aes(value,fill=PS))+
  geom_density()+
  geom_vline(data=cutoff_df, aes(xintercept = cutoff_value), color='blue')+
  scale_fill_manual(values=ps_cols)+
  facet_wrap(~PS,scales = 'free',labeller = my_labels)
## @knitr END


run <- 1
ps_comparison <- 
  map(c('S'), function(scenario){
    map(ps_range, function(PS){
      tmp <- get_data(parameter_space = PS, scenario = scenario, experiment = '001', run = run, use_sqlite = F, tables_to_get = 'summary_general')
      return(tmp[[1]])
    }) %>% bind_rows()
  }) %>% bind_rows()

x <- unique(ps_comparison$time)
layer_map <- tibble(layer=1:length(x), time=x)

layers_to_include <- 1:100

network_S_01 <- get_network_structure(ps = '01',scenario = 'S', exp='001', run=run, cutoff_prob = 0.9, layers_to_include=layers_to_include, parse_interlayer=F, plotit=F, folder='/media/Data/PLOS_Biol/')
network_S_02 <- get_network_structure(ps = '02',scenario = 'S', exp='001', run=run, cutoff_prob = 0.9, layers_to_include=layers_to_include, parse_interlayer=F, plotit=F, folder='/media/Data/PLOS_Biol')
network_S_03 <- get_network_structure(ps = '03',scenario = 'S', exp='001', run=run, cutoff_prob = 0.9, layers_to_include=layers_to_include, parse_interlayer=F, plotit=F, folder='/media/Data/PLOS_Biol')
edge_weights_df_S <- NULL
for (ps in ps_range){
  x <- get(paste('network_S_',ps,sep=''))
  x <- x$temporal_network
  for (i in layers_to_include){
    print(paste(ps,i,sep=' | '))
    if (class(x[[i]])!='igraph') {next} # If there is no layer because there was exticntion.
    tmp <- data.frame(ps=ps, layer=i, weight=E(x[[i]])$weight)
    edge_weights_df_S <- rbind(edge_weights_df_S, tmp)
  }
}

as.tibble(edge_weights_df_S) %>% 
  mutate(scenario=factor(scenario, levels=c('S','G'))) %>% 
  ggplot(aes(weight, fill=ps))+
  geom_density(alpha=0.6) +
  scale_fill_manual(values=ps_cols)+
  facet_grid(~ps)+
  mytheme



edge_weights_df_S$scenario <- 'S'
edge_weights_df_G$scenario <- 'G'
edge_weights_df <- rbind(edge_weights_df_S,edge_weights_df_G)

## @knitr Edge_weights_distributions_Plot
as.tibble(edge_weights_df) %>% 
  # filter(exp != '001') %>% 
  filter(layer %in% seq(1,120,by = 12)) %>% 
  # filter(scenario == 'G') %>%
  mutate(scenario=factor(scenario, levels=c('S','G'))) %>% 
  ggplot(aes(weight, fill=scenario, y=..scaled..))+
  geom_density(alpha=0.6) +
  scale_fill_manual(values=scenario_cols)+
  facet_grid(exp~layer)+
  mytheme

as.tibble(edge_weights_df) %>% 
  # filter(exp != '001') %>% 
  filter(layer %in% 26:38) %>% 
  # filter(scenario == 'G') %>%
  mutate(scenario=factor(scenario, levels=c('S','G'))) %>% 
  ggplot(aes(weight, fill=scenario, y=..scaled..))+
  geom_density(alpha=0.6) +
  scale_fill_manual(values=scenario_cols)+
  facet_grid(exp~layer)+
  mytheme
## @knitr END

exp_comparison %>%
  mutate(scenario=factor(scenario, levels=c('S','G'))) %>% 
  left_join(layer_map) %>% 
  select(-year, -month, -n_infected) %>% 
  filter(time>time_range[1]&time<time_range[2]) %>%
  gather(variable, value, -time, -layer, -exp, -PS, -scenario, -run, -pop_id) %>% 
  filter(variable %in% c('n_circulating_strains','n_circulating_genes')) %>%
  group_by(scenario, time, layer, PS, exp, variable) %>% # Need to average across runs
  summarise(value_mean=mean(value)) %>% 
  ggplot(aes(x=layer, y=value_mean, color=exp))+
  geom_line()+
  facet_wrap(scenario~variable, scales = 'free')+
  geom_vline(xintercept = c(12+c(0,24,60,120)), linetype='dashed')+
  scale_color_manual(values=exp_cols)+
  # scale_x_continuous(breaks=pretty(x=subset(d, time>time_range[1]&time<time_range[2])$time,n=5))+
  mytheme




# Allelic diversity within a genome --------------------------

## @knitr Allele_diversity_load

ps <- '36'
scenario <- 'G'
exp <- '002'
run <- 1

get_genetic_composition <- function(ps, scenario, exp, run, layers_to_include=1:300){
  # Define the sqlite file to use
  base_name <- paste('PS',ps,'_',scenario,'_E',exp,'_R',run,sep='')
  if (on_Midway()){
    sqlite_file <- paste('sqlite/',base_name,'.sqlite',sep='')
  } else {
    sqlite_file <- paste('/media/Data/malaria_interventions_data/sqlite_',scenario,'/',base_name,'.sqlite',sep='')
  }
  
  # Extract data from sqlite. variable names correspond to table names
  db <- dbConnect(SQLite(), dbname = sqlite_file)
  print('Getting genetic data from sqlite...')
  sampled_strains <- as.tibble(dbGetQuery(db, 'SELECT id, gene_id FROM sampled_strains'))
  names(sampled_strains)[1] <- c('strain_id')
  sampled_alleles <- as.tibble(dbGetQuery(db, 'SELECT * FROM sampled_alleles'))
  names(sampled_alleles)[3] <- c('allele_id')
  sampled_strains <- full_join(sampled_strains, sampled_alleles)
  sampled_strains$allele_locus <- paste(sampled_strains$allele_id,sampled_strains$locus,sep='_') # each allele in a locus is unique
  
  # Get sampled infections
  sampled_infections <- get_data(parameter_space = ps, experiment = exp, scenario = scenario, run = run)$sampled_infections
  x <- unique(sampled_infections$time)
  layer_map <- tibble(layer=1:length(x), time=x)
  sampled_infections %<>% left_join(layer_map)
  
  result <- sampled_infections %>% 
    filter(layer %in% layers_to_include) %>% 
    left_join(sampled_strains, by='strain_id')
  return(result)
}

genetic_composition_S <- get_genetic_composition(ps = '36', scenario = 'S', exp = '003', run = 50)
genetic_composition_G <- get_genetic_composition(ps = '36', scenario = 'G', exp = '003', run = 50)

## @knitr Allele_diversity_plot

genetic_composition_S %>% bind_rows(genetic_composition_G) %>% 
  mutate(scenario=factor(scenario, levels = c('S','G'))) %>% 
  filter(layer %in% 1:75) %>%
  select(scenario, time, layer, strain_id, allele_locus) %>% 
  group_by(scenario, time, layer,strain_id) %>% 
  summarise(num_unique_alleles=length(unique(allele_locus))) %>% 
  group_by(scenario, time, layer) %>% 
  summarise(mean_unique_alleles=mean(num_unique_alleles)) %>% 
  ggplot(aes(x=layer,y=mean_unique_alleles, color=scenario))+ 
  geom_vline(xintercept = c(12,24,72),color='purple')+
  scale_color_manual(values=scenario_cols) +
  geom_line()+
  mytheme

## @knitr END


# Network structure across runs -------------------------------------------
## @knitr Network_structure_load

layers_to_include <- 1:200
PS <- '36'
run <- 1

# network_structure_S <- analyze_networks_multiple(ps = '36',scenario = 'S',runs = 1, layers_to_include = layers_to_include, parse_interlayer = F)
# write.csv(network_structure_S,paste('/media/Data/malaria_interventions_data/Results/',PS,'_',scenario,'/',PS,'_',scenario,'_','R',run,'_network_properties.csv',sep=''),row.names = F)

# network_structure_G <- analyze_networks_multiple(ps = '36',scenario = 'G',runs = 1, layers_to_include = layers_to_include, parse_interlayer = F)
# scenario <- 'G'
# write.csv(network_structure_G,paste('/media/Data/malaria_interventions_data/Results/',PS,'_',scenario,'/',PS,'_',scenario,'_','R',run,'_network_properties.csv',sep=''),row.names = F)

scenario <- 'S'
network_structure_S <- read_csv(paste('/media/Data/malaria_interventions_data/Results/',PS,'_',scenario,'/',PS,'_',scenario,'_','R',run,'_network_properties.csv',sep=''))
scenario <- 'G'
network_structure_G <- read_csv(paste('/media/Data/malaria_interventions_data/Results/',PS,'_',scenario,'/',PS,'_',scenario,'_','R',run,'_network_properties.csv',sep=''))


network_properties <- c('Num_nodes','Num_edges','mean_edge_weight','GCC','density','mean_degree','diameter','M11--A<->B<->C','M16--A<->B<->C_A<->C')

## @knitr Network_structure_plot_001
network_structure_S %>% 
  bind_rows(network_structure_G) %>% 
  mutate(scenario=factor(scenario,levels = c('S','G'))) %>% 
  select(c('exp','layer','scenario',network_properties)) %>%
  gather(variable, value, -exp, -layer, -scenario) %>%
  filter(exp=='001') %>% 
  mutate(variable=factor(variable, levels=c('Num_nodes',
                                            'Num_edges',
                                            'density',
                                            'mean_degree',
                                            'mean_edge_weight',
                                            'diameter',
                                            # 'num_edges_il','density_il','mean_edge_weight_il','mean_strength_s,'
                                            'GCC','M11--A<->B<->C','M16--A<->B<->C_A<->C'))) %>% 
  ggplot(aes(layer, value, color=scenario))+
  geom_line()+
  scale_color_manual(values=scenario_cols)+
  scale_x_continuous(breaks = seq(1,max(layers_to_include),20))+
  facet_wrap(~variable, scales='free')+
  geom_vline(xintercept = c(13,13+24,13+60,13+120), color='gray', linetype='dashed')+
  mytheme
  # theme(panel.grid.minor = element_blank(), legend.position = 'none')

## @knitr Network_structure_plot_002
network_structure_S %>% 
  bind_rows(network_structure_G) %>% 
  mutate(scenario=factor(scenario,levels = c('S','G'))) %>% 
  select(c('exp','layer','scenario',network_properties)) %>%
  gather(variable, value, -exp, -layer, -scenario) %>%
  filter(exp=='002') %>% 
  mutate(variable=factor(variable, levels=c('Num_nodes',
                                            'Num_edges',
                                            'density',
                                            'mean_degree',
                                            'mean_edge_weight',
                                            'diameter',
                                            # 'num_edges_il','density_il','mean_edge_weight_il','mean_strength_s,'
                                            'GCC','M11--A<->B<->C','M16--A<->B<->C_A<->C'))) %>% 
  ggplot(aes(layer, value, color=scenario))+
  geom_line()+
  scale_color_manual(values=scenario_cols)+
  scale_x_continuous(breaks = seq(1,max(layers_to_include),20))+
  facet_wrap(~variable, scales='free')+
  geom_vline(xintercept = c(13,13+24,13+60,13+120), color='gray', linetype='dashed')+
  mytheme
# theme(panel.grid.minor = element_blank(), legend.position = 'none')

## @knitr END


plotLayer(network_S_003, l = 42, remove.loops = T, edge_weight_multiply = 1, coords = NULL)
plotLayer(network_G_003, l = 42, remove.loops = T, edge_weight_multiply = 1, coords = NULL)
g <- network_G_003$temporal_network[[30]]
cl <- cluster_infomap(as.undirected(g))
plot(cl, simplify(g), vertex.label=NA, vertex.size=4, edge.arrow.width=0.2,edge.arrow.size=0.2,edge.curved=0.5)



# DOI vs. infections curve ------------------------------------------------
PS <- '36'
doi_S <- get_duration_infection(parameter_space = PS, scenario = 'S', experiment = '001', run = 1)
doi_S <- subset(doi_S, time>=28815 & time <=39945)
doi_S$layer <- .bincode(round(doi_S$time), breaks = seq(28815,39945,by = 30))

doi_G <- get_duration_infection(parameter_space = PS, scenario = 'G', experiment = '001', run = 1)
doi_G <- subset(doi_G, time>=28815 & time <=39945)
doi_G$layer <- .bincode(round(doi_G$time), breaks = seq(28815,39945,by = 30))

doi_S$scenario <- 'S'
doi_G$scenario <- 'G'
doi_S %>% 
  bind_rows(doi_G) %>%
  mutate(scenario=factor(scenario, levels = c('S','G'))) %>%
  filter(layer %in% 1:12) %>% 
  filter(infection_id<=400) %>% 
  ggplot(aes(x=infection_id, y=duration, color=scenario))+
    geom_point(alpha=0.6)+
    scale_color_manual(values=scenario_cols)+
    facet_wrap(~layer)+
    mytheme


# Example for structure ---------------------------------------------------

# Plot ------------------------------------
# Create a color table for all the repertoires. All copies of each unique repertoire have the same color
node_names <- sort(unique(unlist(lapply(network_test$temporal_network, rownames))))
color_table <- data.frame(node_name=node_names, unique_name=splitText(node_names,splitchar = '_', after = F))
colors <- data.frame(unique_name=unique(color_table$unique_name), color=gg_color_hue(length(unique(color_table$unique_name))), stringsAsFactors = F)
color_table <- merge(color_table,colors)

# # Get fixed coordinates
# g <- graph.adjacency(network_test$similarityMatrix, mode = 'directed', weighted = T)
# l <- layout.fruchterman.reingold(g)
# coords <- data.frame(node_name=V(g)$name, x=l[,1], y=l[,2])
layers_to_plot <- c(1,7,13,21,36,50,65,70,75,78,85,91,93,100,114,116,120,191)
# pdf('networks_ctrl.pdf', width = 8,height = 8)
for (i in layers_to_plot){
  print(i)
  png(paste('png/layer_',sprintf('%0.3d', i),'.png',sep=''), width=1920, height=1080, res = 96)
  plotLayer(network_test, i, ver.col = color_table, coords = NULL, main=i)
  dev.off()
}
#system('convert -delay 0.5 *.png animation.mpg')
#----------------------------------------------




# Infomap -----------------------------------------------------------------

infomap_readTreeFile <- function(file, reorganize_modules=T, max_layers=372, remove_buggy_instances=T){
  require(splitstackshape)
  lines <- readLines(file)
  #lines <- readLines('Infomap_linux/output/S11_S0_w30_300_300_0.1_expanded.tree');length(lines)
  cat(lines[1]);cat('\n')
  x=read.table(file, skip = 2, stringsAsFactors = F)
  modules <- data.frame(module=rep(NA,nrow(x)),
                        strain=rep(NA,nrow(x)),
                        layer=rep(NA,nrow(x)),
                        flow=rep(NA,nrow(x)), stringsAsFactors = F)
  
  modules$path <- x[,1]
  x.module <- x[,1]
  x.module <- cSplit(as.data.table(x.module),'x.module',':')
  x.module <- as.data.frame(x.module)
  modules$module <- x.module[,1]
  cat(nrow(x),'state nodes','in',paste(max(modules$module),'modules, organized in',length(x.module),'levels...'));cat('\n')
  modules$flow <- x[,2]
  x.strain <- x[,3]
  x.strain <- read.table(text = x.strain, sep = "|", colClasses = "character", stringsAsFactors = F, strip.white = T)
  modules$strain <- x.strain$V1
  modules$layer <- as.numeric(x$V4)
  
  # There is a bug in Infomap that assigns nodes to layers which do not exist (with higher number than the existing layers).
  
  if(remove_buggy_instances){
    buggy <- modules[modules$layer>max_layers,]
    if (nrow(buggy)>=1){
      cat('\n')
      print('---------------------------------')
      print('Some buggy instances encountered!')
      print(paste('removed ', nrow(buggy),' instances which were assigned to layers which do not exist.',sep=''))
      print(paste('total flow of removed instance: ',sum(buggy$flow)))
      print('Buggy instances written to file')
      write.table(buggy, paste(str_sub(file, 1, str_locate(file, 'output/'))[2],'buggy_instances_infomap.txt',sep=''))
      modules <- modules[modules$layer<=max_layers,]
    }
  }
  
  if(reorganize_modules){  # organize the names of modules to be consecutive
    cat('Re-organizing modules by layer...');cat('\t')
    modules=modules[with(modules, order(layer,module)),]
    x=table(modules$module)[unique(modules$module)]
    modules$module_ordered <- NA
    for (i in 1:nrow(modules)){
      modules[i,'module_ordered'] <- which(names(x)==modules[i,'module'])
    }
    modules <- modules[,-1]
    modules <- modules[,c('module_ordered','strain','layer','flow','path')]
    names(modules)[1] <- 'module'
  }
  print('Done!')
  return(modules)
}


infomap_objects <- build_infomap_objects(network_test)

# Get infomap results and process them











# Compare curves firs GI and IS -------------------------------------------

# This part compares the fits of the duration curve of the selection and the generalized immunity.
parameter_space <- '03'
experiment <- '01'
run <- 1

sqlite_file <- paste('/home/shai/Documents/malaria_interventions_sqlite/','PS',parameter_space,'_S_E',experiment,'_R',run,'.sqlite',sep='')
db <- dbConnect(SQLite(), dbname = sqlite_file)
sampled_duration <- dbGetQuery(db, 'SELECT time, duration,infection_id FROM sampled_duration')

setwd('/home/shai/Documents/malaria_interventions')
fit <- set_generalized_immunity(parameter_space=parameter_space, run=run)[[1]]
generalImmunityParams <- python.get('generalImmunityParams')
a=generalImmunityParams[1]
b=generalImmunityParams[2]
c=generalImmunityParams[3]
d=generalImmunityParams[4]

p <- sampled_duration %>% ggplot(aes(infection_id, duration))+
  geom_point()
# Check fit
x.fit <- 0:max(sampled_duration$infection_id)
y.fit <- ((b*exp(-c*x.fit))/(d*x.fit+1)^d)+a
fit <- data.frame(infection_id=x.fit, duration=y.fit)
p <- p+geom_point(data=fit, color='red',size=3)


sqlite_file <- paste('/home/shai/Documents/malaria_interventions_sqlite/','PS',parameter_space,'_G_E',experiment,'_R',run,'.sqlite',sep='')
db <- dbConnect(SQLite(), dbname = sqlite_file)
sampled_duration <- dbGetQuery(db, 'SELECT time, duration,infection_id FROM sampled_duration')
generalized <- sampled_duration %>% group_by(infection_id) %>% summarise(meanDOI=mean(duration))
fit_generalized <-  data.frame(infection_id=generalized$infection_id, duration=generalized$meanDOI)

p <- p+geom_point(data=fit_generalized, color='blue',size=3)

png('scenario_comparison_2.png',1800,1000)
p+mytheme
dev.off()

# This compares the age distribution of infected hosts
d <- rbind(PS03_S_01[[2]],PS03_G_01[[2]],PS03_N_01[[2]])
png('scenario_comparison_3.png',1800,1000)
d %>% ggplot(aes(x=host_age, fill=scenario))+geom_histogram() + 
  labs(x='Infected host age (months)') + 
  geom_vline(xintercept = 60) +
  scale_fill_manual(values=c('blue','orange','red'))+
  mytheme
dev.off()










# Calendar ----------------------------------------------------------------

num_years <- 100
months_in_year <- rep(c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'), each=30)
calendar <- data.frame(running_day=seq(from = 1,to = 360*num_years,by=1),
                       year_sim=rep(1:num_years, each=360),
                       month_sim=rep(months_in_year,num_years),
                       day_sim=rep(1:30,num_years))
# calendar$layer <- ceiling((calendar$running_day-burnin)/30)
# calendar$burnin <- 'No'
# calendar$burnin[1:burnin] <- 'Yes'
calendar <- as_tibble(calendar)
calendar$layer <- .bincode(round(doi_S$time), breaks = seq(28815,39945,by = 30))
calendar %>% 
  filter(running_day>=28815 & running_day <=39945) %>% 
  mutate(layer=.bincode(running_day, breaks = seq(28815,39945,by = 30))) %>% 
  filter(layer==7)
