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
basic_variable_s <- compare_ps(ps_range=c('01','02','03'), 'S', exp = '001', 1, cutoff_prob=c(0.25,0.7,0.9))
basic_variable_s[[1]]

## @knitr basic_variables_neutral
basic_variable_n <- compare_ps(ps_range=c('01','02','03'), 'N', exp = '001', 1, cutoff_prob=c(0.25,0.7,0.9))
basic_variable_n[[1]]

## @knitr EIR_selection
basic_variable_s[[2]]

## @knitr EIR_neutral
basic_variable_n[[2]]

## @knitr END


# Compare between scenarios within a parameter space and experiment -------

## @knitr compare_scenarios_01
x <- compare_scenarios(PS = '01', scenarios = c('S','N'), run_range = 1, cutoff_prob = 0.25)
x
## @knitr compare_scenarios_02
x <- compare_scenarios(PS = '02', scenarios = c('S','N'), run_range = 1, cutoff_prob = 0.7)
x
## @knitr compare_scenarios_03
x <- compare_scenarios(PS = '03', scenarios = c('S','N'), run_range = 1, cutoff_prob = 0.9)
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
  geom_vline(data=cutoff_df, aes(xintercept = cutoff_value), color='red')+
  scale_fill_manual(values=ps_cols)+
  facet_wrap(~PS,scales = 'free',labeller = my_labels)+
  theme_bw(base_size=26)
dev.off()
## @knitr END

# NOTE THAT IT IS IMPOSSIBLE TO PLOT EDGE WEIGHTS FOR NEUTRALITY: THERE IS NOT ENOUGH MEMORY TO LOAD ALL THE INTERACTIONS
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



# Infomap -----------------------------------------------------------------
file <- '/media/Data/PLOS_Biol/Results/03_S/PS03_S_E001_R1_0.9_Infomap_multilayer_expanded.tree'
infomap_readTreeFile <- function(file){
  lines <- readLines(file)
  cat(lines[1]);cat('\n')
  x <- fread(file, skip = 2, stringsAsFactors = F) # Read results of Infomap
  
  # Create a data frame to store results
  modules <- tibble(module=rep(NA,nrow(x)),
                        strain=rep(NA,nrow(x)),
                        layer=rep(NA,nrow(x)),
                        path=x$V1)
  
  modules$module <- as.numeric(str_split(string = modules$path, pattern = ':', simplify = T)[,1])
  modules$strain <- str_trim(str_split(string = x$V3, pattern = '\\|', simplify = T)[,1])
  modules$layer <- as.numeric(str_trim(str_split(string = x$V3, pattern = '\\|', simplify = T)[,2]))
  
  cat(nrow(x),'state nodes','in',paste(max(modules$module),'modules'));cat('\n')
  
  # Rename modules because Infomap gives random names
  modules2 <- modules %>% 
    distinct(module,layer) %>% 
    arrange(module,layer)
  x <- c(1,table(modules2$module))
  module_birth_layers <- modules2 %>% slice(cumsum(x)) %>% arrange(layer,module)
  module_renaming <- data.frame(module=module_birth_layers$module, module_renamed = 1:max(module_birth_layers$module)) 
  
  modules2 %<>% left_join(module_renaming)
  modules2 %<>% full_join(modules) 
  modules2 %<>% select(-module, -path)
  names(modules2)[2] <- 'module'
  modules2 %<>% arrange(module, layer, strain)
 return(modules2)
}

modules <- infomap_readTreeFile(file)
modules %>% 
  distinct(module,layer) %>% 
  ggplot(aes(x=layer, y=module, color=module))+
  geom_point(size=2)+
  # scale_x_continuous(breaks = tickBreaks)+
  # scale_y_discrete(breaks = tickBreaks_y)+
  labs(y= 'module ID', x='Time (months)')

# Files to calculate statistic
write_csv(modules, '~/Dropbox/Qixin_Shai_Malaria/PS03_S_E001_R1_0.9_modules.csv')

if (use_sqlite){
  base_name <- paste('PS',parameter_space,'_',scenario,'_E',experiment,'_R',run,sep='')
  if (on_Midway()){
    sqlite_file <- paste('/scratch/midway2/pilosofs/PLOS_Biol/sqlite/',base_name,'.sqlite',sep='')
    print(sqlite_file)
    print(file.exists(sqlite_file))
  } else {
    sqlite_file <- paste('/media/Data/PLOS_Biol/sqlite_',scenario,'/',base_name,'.sqlite',sep='')
  }
  
  if (!file.exists(sqlite_file)) {
    print (paste(sqlite_file, ' does not exist, ignoring and returning NULL'))
    return(NULL)
  }
  # parameter_file <- paste(base_name,'.py',sep='') # This may be necessary so I leave it
  
  # Extract data from sqlite. variable names correspond to table names
  print('Connecting to sqlite file...')
  db <- dbConnect(SQLite(), dbname = sqlite_file)
  summary_general <- dbGetQuery(db, 'SELECT * FROM summary')
  summary_general$PS <- parameter_space
  summary_general$exp <- experiment
  summary_general$scenario <- scenario
  summary_general$run <- run
  summary_general$year <- gl(n = max(summary_general$time)/360, length = nrow(summary_general), k = 1)
  summary_general$month <- gl(n = 12, k = 1, length = nrow(summary_general),labels = c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'), ordered = F)
  
  summary_alleles <- dbGetQuery(db, 'SELECT * FROM summary_alleles')
  summary_alleles %<>% group_by(time) %>% summarise(n_alleles=sum(n_circulating_alleles))
  
  summary_general <- suppressMessages(left_join(summary_general, summary_alleles))
  
  if ('sampled_infections'%in%tables_to_get){
    sampled_infections <- dbGetQuery(db, 'SELECT * FROM sampled_infections')
    sampled_infections$PS <- parameter_space
    sampled_infections$exp <- experiment
    sampled_infections$scenario <- scenario
    sampled_infections$run <- run
  }
  




modulePersistence <- modules %>% group_by(module_renamed) %>% summarise(birth_layer=min(layer), death_layer=max(layer), persistence=death_layer-birth_layer+1)
modulePersistence %>% ggplot(aes(persistence))+geom_density()



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
