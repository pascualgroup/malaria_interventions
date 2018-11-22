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
scenario_cols <- c('red','orange','blue')


compare_ps <- function(ps_range=c('01','02','03'), scenario, exp, run_range, cutoff_prob=c(0.25,0.7,0.9)){
  
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

# Compare between parameter spaces within an experiment -------------------

## @knitr basic_variables_S
basic_variable_S <- compare_ps(ps_range=c('01','02','03'), 'S', exp = '001', 1, cutoff_prob=c(0.25,0.7,0.9))
basic_variable_S[[1]]

## @knitr basic_variables_N
basic_variable_N <- compare_ps(ps_range=c('01','02','03'), 'N', exp = '001', 1, cutoff_prob=c(0.25,0.7,0.9))
basic_variable_N[[1]]

## @knitr basic_variables_G
basic_variable_G <- compare_ps(ps_range=c('01','02','03'), 'G', exp = '001', 1, cutoff_prob=c(0.25,0.7,0.9))
basic_variable_G[[1]]

## @knitr EIR_S
basic_variable_S[[2]]

## @knitr EIR_N
basic_variable_N[[2]]

## @knitr EIR_G
basic_variable_G[[2]]
## @knitr END


# Compare between scenarios within a parameter space and experiment -------

## @knitr compare_scenarios_01
x <- compare_scenarios(PS = '01', scenarios = c('S','N','G'), run_range = 1, cutoff_prob = 0.25)
x
## @knitr compare_scenarios_02
x <- compare_scenarios(PS = '02', scenarios = c('S','N','G'), run_range = 1, cutoff_prob = 0.7)
x
## @knitr compare_scenarios_03
x <- compare_scenarios(PS = '03', scenarios = c('S','N','G'), run_range = 1, cutoff_prob = 0.9)
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
edges_N_03 <- get_edge_disributions(PS = '03',scenario = 'N',exp = '001',1, 0.9, get_inter = F)
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

PS <- design_basic[i,1]
scenario <- design_basic[i,2]
exp <- design_basic[i,3]
run <- design_basic[i,4]
cutoff_prob <- design_basic[i,5]
design_basic <- expand.grid(PS=sprintf('%0.2d', 1:6),
                            scenario=c('S','N','G'), 
                            exp='001',
                            run_range=1, 
                            stringsAsFactors = F)
design_basic$cutoff_prob <- rep(c(0.25,0.7,0.9),length(unique(design_basic$scenario))*length(unique(run_range)))

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
    x <- read_csv(file)  
    module_results <- rbind(module_results, x)
  }
}

module_results %>% 
  distinct(module,layer,PS, scenario) %>% 
  mutate(scenario=factor(scenario, levels=c('S','N','G'))) %>% 
  ggplot(aes(x=layer, y=module, color=scenario))+
  geom_point(size=2)+
  scale_color_manual(values = scenario_cols)+
  labs(y= 'module ID', x='Time (months)')+
  facet_grid(PS~scenario, scales='free')


# Module persistence
module_persistence <- module_results %>% 
  mutate(scenario=factor(scenario, levels=c('S','N','G'))) %>% 
  group_by(PS,scenario,module) %>% 
  summarise(birth_layer=min(layer), death_layer=max(layer), persistence=death_layer-birth_layer+1)
module_persistence %>% 
  filter(PS %in% c('01','02','03')) %>% 
  ggplot(aes(persistence, y=..count.., fill=scenario))+
  geom_density()+
  geom_rug(aes(x=persistence, y=0), position = position_jitter(height = 0))+
  scale_fill_manual(values = scenario_cols)+
  labs(x='Persistence (months)')+
  facet_grid(PS~scenario, scales='free_y')

# Strain persistence. For that need to cluster the strains
PS <- '03'
scenario <- 'S'
modules <- read_csv(paste('/media/Data/PLOS_Biol/Results/',PS,'_',scenario,'/PS',PS,'_',scenario,'_E',exp,'_R',run,'_',cutoff_prob,'_modules.csv',sep=''))
sampled_strains <- read_csv(paste('/media/Data/PLOS_Biol/Results/',PS,'_',scenario,'/PS',PS,'_',scenario,'_E',exp,'_R',run,'_',cutoff_prob,'_sampled_strains.csv',sep=''))
sampled_alleles <- read_csv(paste('/media/Data/PLOS_Biol/Results/',PS,'_',scenario,'/PS',PS,'_',scenario,'_E',exp,'_R',run,'_',cutoff_prob,'_sampled_alleles.csv',sep=''))

setequal(sampled_strains$strain_id, modules$strain_id)

sampled_alleles$allele_locus <- paste(sampled_alleles$allele, sampled_alleles$locus,sep='_')
sampled_alleles %<>% select(gene_id, allele_locus) %>% arrange(gene_id,allele_locus)
sampled_strains %<>% left_join(sampled_alleles) %>% distinct(strain_id,allele_locus)



y <- xtabs(~strain_id+allele_locus, sampled_strains)
y <- matrix(y, nrow = nrow(y), ncol=ncol(y), dimnames = list(rownames(y),colnames(y)))

# Make a file in format:
# >Otu1
# ATTTAAATTCCTTTTAGGATTAAT
# >Otu2
# TTCCGTGTAACCTAGAACTTTCAATTCTATAGTAGATTAT
# Where Otu is the strain and ATGC are the 0 and 1 in alleles (say A is 0, G is 1).

#Then classify it using: 
#system("./usearch10.0.240_i86linux32 -cluster_fast Yellowstone_unique_no_singletons.fa -id 0.8 -centroids Yellowstone_otus.fa -uc Yellowstone_uc.txt -relabel Otu") #find clusters of seqs

# Then cluster strains to OTUs using:
#system('./usearch10.0.240_i86linux32 -otutab Yellowstone_all_spacers_renamed.fasta -db Yellowstone_otus.fa -otutabout Yellowstone_otutab.txt -id 0.8 -dbmatched Yellowstone_otus_with_sizes.fa -sizeout -notmatched Mutnovsky_notmatched.fa -sample_delim ";"') #build otu table

otutab <- as.data.frame(y)
otutab <- cbind(rownames(otutab), otutab)
rownames(otutab) <- NULL
names(otutab)[1] <- '#OTU ID'
otutab[1:5,1:5]
fwrite(otutab, 'otutab.txt', sep = '\t', quote = F)

# Cluster Yellowstone using USEARCH





strain_persistence <- module_results %>% 
  mutate(scenario=factor(scenario, levels=c('S','N','G'))) %>% 
  group_by(PS,scenario,strain_id) %>% 
  summarise(birth_layer=min(layer), death_layer=max(layer), persistence=death_layer-birth_layer+1)
strain_persistence %>% ggplot(aes(persistence, fill=scenario))+
  geom_density()+
  scale_fill_manual(values = scenario_cols)+
  labs(x='Persistence (months)')+
  facet_grid(PS~scenario, scales='free_y')

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
