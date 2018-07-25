library(tidyverse)
library(magrittr)
library(sqldf)

# Functions ---------------------------------------------------------------
mytheme <- theme_bw() + theme(
  legend.title  = element_text(colour = "black", size=17),
   # legend.position = "none",
  #	legend.direction = "horizontal",
  legend.key = element_blank(),
  legend.text  = element_text(colour = "black", size=17),
  panel.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.text = element_text(color='black', family="Helvetica", size=10),
  strip.text.x = element_text(family = "Helvetica", size = 10),
  strip.text.y = element_text(family = "Helvetica", size = 10),
  panel.border = element_rect(colour = "black", size=1.3),
  axis.ticks = element_line(size = 1.3),
  strip.background = element_rect( fill = "transparent", size = 1.3, colour = "black"  ),
  strip.text = element_text(size = 19)
)

gg_color_hue <- function(n, hue_min = 10, hue_max = 280, l = 62, c = 100) {
  hues = seq(hue_min, hue_max, length=n+1)
  hcl(h=hues, l=l, c=c)[1:n]
}

chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 

# This function extracts the biting rates from the parameter file, which are in
# daily resolution and returns a vector of the averaged biting rates for a given
# period. (e.g. 30 would be biting rates averaged over a month and 1 will return
# the daily biting rates).
get_biting_rate <- function(parameter_file, sampling_period=30){
  x <- readLines(parameter_file)
  y <- x[grep('BITING_RATE_MEAN',x)[1]]
  BITING_RATE_MEAN <- parse_number(y)
  y <- x[grep('DAILY_BITING_RATE_DISTRIBUTION',x)[1]]
  DAILY_BITING_RATE_DISTRIBUTION <- eval(parse(text=paste('c(',(str_sub(y, str_locate(y, '\\[(.*?)\\]')[1]+1, str_locate(y, '\\[(.*?)\\]')[2]-1)),')',sep='')))
  BITING_RATE <- BITING_RATE_MEAN*DAILY_BITING_RATE_DISTRIBUTION
  BITING_RATE <- chunk2(BITING_RATE, 360/sampling_period)
  sapply(BITING_RATE, mean)
}

get_data <- function(parameter_space, scenario, experiment, run, sampling_period=30){
  require(sqldf)
  # Initialize
  base_name <- paste('PS',parameter_space,'_',scenario,'_E',experiment,'_R',run,sep='')
  sqlite_file <- paste(base_name,'.sqlite',sep='')
  # parameter_file <- paste(base_name,'.py',sep='') # This may be necessary so I leave it
  
  # Extract data from sqlite. variable names correspond to table names
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
  
  summary_general %<>% left_join(summary_alleles)
  
  sampled_infections <- dbGetQuery(db, 'SELECT * FROM sampled_infections')
  sampled_infections$PS <- parameter_space
  sampled_infections$exp <- experiment
  sampled_infections$scenario <- scenario
  sampled_infections$run <- run
  
  # Extract statistics
  if (nrow(summary_general)%%sampling_period==1){# The beginning of the data set in E00 has a timestap of 0 and this line is unnecesary. So remove it
    summary_general <- summary_general[-1,]
  }

  ## Prevalence
  summary_general$prevalence <- summary_general$n_infected/10^4
  
  ## Or get EIR from the table in the sqlite
  summary_general$EIR <- summary_general$n_infected_bites/10000 # 10000 is the size of the human population
  
  # ##EIR can also be calculated manually as EIR=biting_rate * prevalence
  # biting_rate <- get_biting_rate(parameter_file)
  # summary_general$b <- NA
  # summary_general$b[] <- biting_rate # The [] is for recycling the biting_rate
  # summary_general$EIR1 <- summary_general$prevalence*summary_general$b*sampling_period
  
  # # Annual EIR is given by dividing by the sampling period to get a daily EIR, then multiplying by 360
  # summary_general %>% mutate(eir_y=EIR2/sampling_period*360) %>% group_by(year) %>% summarise(EIR_Year=mean(eir_y))

  # MOI#
  meanMOI <- sampled_infections %>% group_by(time, host_id) %>% summarise(MOI=length(strain_id)) %>% group_by(time) %>% summarise(meanMOI=mean(MOI))
  summary_general <- inner_join(summary_general, meanMOI) # Add MOI to the summary
  
  # Host age structure
  hosts <- dbGetQuery(db, 'SELECT * FROM hosts')
  names(hosts)[1] <- 'host_id'
  hosts$lifespan <- round((hosts$death_time-hosts$birth_time))
  hosts <- subset(hosts, host_id%in%sampled_infections$host_id)
  sampled_infections <- left_join(sampled_infections, hosts, by='host_id')
  sampled_infections$host_age <- round((sampled_infections$time-sampled_infections$birth_time)/30)
  
  return(list(summary_general=as.tibble(summary_general), sampled_infections=as.tibble(sampled_infections)))
}



get_data2 <- function(parameter_space, scenario, experiment, run, sampling_period=30){
  require(sqldf)
  # Initialize
  base_name <- paste('PS',parameter_space,'_',scenario,'_E',experiment,'_R',run,sep='')
  sqlite_file <- paste(base_name,'.sqlite',sep='')
  # parameter_file <- paste(base_name,'.py',sep='') # This may be necessary so I leave it
  
  # Extract data from sqlite. variable names correspond to table names
  db <- dbConnect(SQLite(), dbname = sqlite_file)
  hosts <- as.tibble(dbGetQuery(db, 'SELECT * FROM hosts'))
  names(hosts)[1] <- 'host_id'
  sampled_alleles <- as.tibble(dbGetQuery(db, 'SELECT * FROM sampled_alleles'))
  sampled_duration <- as.tibble(dbGetQuery(db, 'SELECT * FROM sampled_duration'))
  sampled_infections <- as.tibble(dbGetQuery(db, 'SELECT * FROM sampled_infections'))
  sampled_genes <- as.tibble(dbGetQuery(db, 'SELECT * FROM sampled_genes'))
  names(sampled_genes)[1] <- 'gene_id'
  sampled_strains <- as.tibble(dbGetQuery(db, 'SELECT * FROM sampled_strains'))
  names(sampled_strains)[1] <- 'strain_id'
  
  summary_general <- left_join(sampled_infections,hosts)
  data %>% group_by(time) %>% summarise(n_infections=length(infection_id))
  
  
  
  summary_general <- dbGetQuery(db, 'SELECT * FROM summary')
  summary_general$PS <- parameter_space
  summary_general$exp <- experiment
  summary_general$scenario <- scenario
  summary_general$run <- run
  summary_general$year <- gl(n = max(summary_general$time)/360, length = nrow(summary_general), k = 1)
  summary_general$month <- gl(n = 12, k = 1, length = nrow(summary_general),labels = c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'), ordered = F)
  
  summary_alleles <- dbGetQuery(db, 'SELECT * FROM summary_alleles')
  summary_alleles %<>% group_by(time) %>% summarise(n_alleles=sum(n_circulating_alleles))
  
  summary_general %<>% left_join(summary_alleles)
  
  sampled_infections <- dbGetQuery(db, 'SELECT * FROM sampled_infections')
  sampled_infections$PS <- parameter_space
  sampled_infections$exp <- experiment
  sampled_infections$scenario <- scenario
  sampled_infections$run <- run
  
  # Extract statistics
  if (nrow(summary_general)%%sampling_period==1){# The beginning of the data set in E00 has a timestap of 0 and this line is unnecesary. So remove it
    summary_general <- summary_general[-1,]
  }
  
  ## Prevalence
  summary_general$prevalence <- summary_general$n_infected/10^4
  
  ## Or get EIR from the table in the sqlite
  summary_general$EIR <- summary_general$n_infected_bites/10000 # 10000 is the size of the human population
  
  # ##EIR can also be calculated manually as EIR=biting_rate * prevalence
  # biting_rate <- get_biting_rate(parameter_file)
  # summary_general$b <- NA
  # summary_general$b[] <- biting_rate # The [] is for recycling the biting_rate
  # summary_general$EIR1 <- summary_general$prevalence*summary_general$b*sampling_period
  
  # # Annual EIR is given by dividing by the sampling period to get a daily EIR, then multiplying by 360
  # summary_general %>% mutate(eir_y=EIR2/sampling_period*360) %>% group_by(year) %>% summarise(EIR_Year=mean(eir_y))
  
  # MOI#
  meanMOI <- sampled_infections %>% group_by(time, host_id) %>% summarise(MOI=length(strain_id)) %>% group_by(time) %>% summarise(meanMOI=mean(MOI))
  summary_general <- inner_join(summary_general, meanMOI) # Add MOI to the summary
  
  # Host age structure
  hosts <- dbGetQuery(db, 'SELECT * FROM hosts')
  names(hosts)[1] <- 'host_id'
  hosts$lifespan <- round((hosts$death_time-hosts$birth_time))
  hosts <- subset(hosts, host_id%in%sampled_infections$host_id)
  sampled_infections <- left_join(sampled_infections, hosts, by='host_id')
  sampled_infections$host_age <- round((sampled_infections$time-sampled_infections$birth_time)/30)
  
  return(list(summary_general=as.tibble(summary_general), sampled_infections=as.tibble(sampled_infections)))
}


# This function uses the list obtained by get_data() and generates relevant plots
generate_plots <- function(data, time_range=NULL){
  data1 <- data[[1]]
  data2 <- data[[2]]
  if (is.null(time_range)){
    time_range <- c(min(data1$time), max(data1$time))
  }
  plot_variables <- data1 %>% 
    select(-n_infected, -year, -month, -PS, -exp, -run, -scenario) %>% 
    filter(time>time_range[1]&time<time_range[2]) %>%    gather(variable, value, -time) %>% 
    ggplot(aes(time, value, color=variable))+
    geom_line()+
    facet_wrap(~variable, scales = 'free')+
    mytheme+theme(legend.position = 'none')
  
  plot_eir <- data1 %>% 
    ggplot(aes(x=month,y=EIR))+
    geom_boxplot()+
    geom_point(stat='summary', fun.y=mean, color='red')+
    stat_summary(fun.y=mean, geom="line")+mytheme
  
  # Plot the age structure of infected hosts
  plot_age_structure <- data2 %>% group_by()%>% ggplot(aes(x=host_age))+geom_histogram() + 
    labs(x='Infected host age (months)') + 
    geom_vline(xintercept = 60) +
    mytheme
  return(list(plot_variables, plot_eir, plot_age_structure))
}


# This function compares an experiment to control and calculates statistics for the intervention
post_intervention_stats <- function(PS, scenario='S', exp, run, post_intervention_lag=360, control_data=NULL, plot.it=F, design_irs){
  
  plots_out <- list()
  
  # If a data frame with control data is not included (NULL), then load the
  # control experiment
  if (is.null(control_data)){
    control_data <- get_data(parameter_space = PS, scenario = scenario, experiment = '001', run = run)[[1]]
  }
  
  # Calculate variable means (across time) for control
  control_means <- control_data %>% select(-year, -month, -n_infected) %>% 
    gather(variable, value, -time, -exp, -PS, -scenario, -run, -pop_id) %>% 
    group_by(PS, exp, variable) %>% summarise(mean_value=mean(value, na.rm = T))
  
  # Plot control and means
  if (plot.it){
    plots_out$control <- control_data %>% select(-year, -month, -n_infected) %>%
      gather(variable, value, -time, -exp, -PS, -scenario, -run, -pop_id) %>%
      ggplot(aes(x=time, y=value))+
      geom_line()+
      facet_wrap(~variable,scales='free')+
      geom_hline(aes(yintercept=mean_value), data=control_means, color='blue')+
      mytheme
  }
  
  # Load experiment data
  experiment_data <- get_data(parameter_space = PS, scenario = scenario, experiment = exp, run = run)[[1]]
  
  # Plot control vs experiment
  if (plot.it){
    plots_out$control_experiment <- bind_rows(control_data, experiment_data) %>%
      select(-year, -month, -n_infected) %>% 
      gather(variable, value, -time, -exp, -PS, -scenario, -run, -pop_id) %>% 
      ggplot(aes(x=time, y=value, color=exp))+
      geom_vline(xintercept = c(29160,29160+1800), linetype='dashed', color='tomato')+
      geom_line()+
      scale_color_manual(values=c('#504E4E','purple'))+
      scale_x_continuous(breaks=pretty(x=subset(d, time>time_range[1]&time<time_range[2])$time,n=5))+
      geom_hline(aes(yintercept=mean_value), data=control_means, color='tomato')+
      facet_wrap(~variable, scales = 'free')+mytheme
  }
  
  # Calculate absolute difference from control at every timepoint
  x <- control_data %>%
    select(-year, -month, -n_infected) %>% 
    gather(variable, value, -time, -exp, -PS, -scenario, -run, -pop_id)
  y <- experiment_data %>%
    select(-year, -month, -n_infected) %>% 
    gather(variable, value, -time, -exp, -PS, -scenario, -run, -pop_id)
  diff_control <- inner_join(x,y,by=c('pop_id','PS','scenario','run','time','variable')) %>% 
    mutate(diff=value.y-value.x, ratio=value.y/value.x) %>% 
    select(variable, time, diff, ratio)
  
  if (plot.it){
    plots_out$control_diff <- diff_control %>%
      ggplot(aes(time, diff))+
      geom_line()+
      facet_wrap(~variable,scales='free')+
      geom_hline(aes(yintercept=mean_value), data=control_means, color='blue')+
      geom_vline(xintercept = c(29160,29160+3600), color='red')+
      mytheme
  }
  
  # Maximum amplitudes after an intervention
  amplitude <- experiment_data %>% left_join(subset(design_irs, select=c(IRS_START_TIMES,IRS_length,exp)), by='exp') %>% 
    mutate(intervention_lift=as.numeric(IRS_START_TIMES)+as.numeric(IRS_length)) %>% 
    filter(time>intervention_lift) %>% 
    select(-year, -month, -n_infected) %>% 
    gather(variable, value, -time, -exp, -PS, -scenario, -run, -pop_id,-IRS_START_TIMES,-IRS_length) %>% 
    group_by(variable) %>% summarise(max_value=max(value))
  
  # New mean from some lag post intervention
  new_mean <- experiment_data %>% left_join(subset(design_irs, select=c(IRS_START_TIMES,IRS_length,exp)), by='exp') %>% 
    mutate(intervention_lift=as.numeric(IRS_START_TIMES)+as.numeric(IRS_length)) %>% 
    filter(time>intervention_lift+post_intervention_lag) %>% 
    select(-year, -month, -n_infected) %>% 
    gather(variable, value, -time, -exp, -PS, -scenario, -run, -pop_id,-IRS_START_TIMES,-IRS_length) %>% 
    group_by(variable) %>% summarise(mean_value=mean(value, na.rm = T))
  
  # Time when extinction occured (last time point with data)
  time_extinct <- experiment_data %>% left_join(subset(design_irs, select=c(IRS_START_TIMES,IRS_length,exp)), by='exp') %>% 
    select(-year, -month, -n_infected) %>% 
    gather(variable, value, -time, -exp, -PS, -scenario, -run, -pop_id,-IRS_START_TIMES,-IRS_length) %>% 
    group_by(variable) %>% summarise(time_extinct=max(time))
  
  summary_stats <- left_join(time_extinct, new_mean)
  summary_stats <- left_join(summary_stats, amplitude)
  summary_stats$PS <- PS
  summary_stats$scenario <- scenario
  summary_stats$exp <- exp
  summary_stats$run <- run  
  
  diff_control$PS <- PS
  diff_control$scenario <- scenario
  diff_control$exp <- exp
  diff_control$run <- run 
  
  if (plot.it){
    return(list(plots=plots_out,stats=list(diff_control=diff_control, summary_stats=summary_stats)))
  } else {
    return(list(diff_control=diff_control, summary_stats=summary_stats))  
  }
}


# Get data trials ---------------------------------------------------------
setwd('~/Documents/malaria_interventions_data/')
PS03_S_00 <- get_data(parameter_space = '03', scenario = 'S', experiment = '00', run = 1)
PS03_S_01 <- get_data(parameter_space = '03', scenario = 'S', experiment = '01', run = 1)
PS03_S_02 <- get_data(parameter_space = '03', scenario = 'S', experiment = '02', run = 1)
PS03_S_03 <- get_data(parameter_space = '03', scenario = 'S', experiment = '03', run = 1)

PS03_N_00 <- get_data(parameter_space = '03', scenario = 'N', experiment = '00', run = 1)
PS03_N_01 <- get_data(parameter_space = '03', scenario = 'N', experiment = '01', run = 1)
PS03_N_02 <- get_data(parameter_space = '03', scenario = 'N', experiment = '02', run = 1)
PS03_N_03 <- get_data(parameter_space = '03', scenario = 'N', experiment = '03', run = 1)

PS03_G_00 <- get_data(parameter_space = '03', scenario = 'G', experiment = '00', run = 1)
PS03_G_01 <- get_data(parameter_space = '03', scenario = 'G', experiment = '01', run = 1)
PS03_G_02 <- get_data(parameter_space = '03', scenario = 'G', experiment = '02', run = 1)
PS03_G_03 <- get_data(parameter_space = '03', scenario = 'G', experiment = '03', run = 1)

plots <- generate_plots(PS03_S_01)
splot <- plots[[3]]  
gplot <- plots[[3]]  
nplot <- plots[[3]]  
splot

rbind(PS03_S_01)

PS02_S_00 <- get_data(parameter_space = '02', scenario = 'S', experiment = '00', run = 1)
PS02_S_01 <- get_data(parameter_space = '02', scenario = 'S', experiment = '01', run = 1)
PS02_S_02 <- get_data(parameter_space = '02', scenario = 'S', experiment = '02', run = 1)

for (i in sprintf('%0.2d', 0:5)){
  print(i)
  assign(paste('PS04_S_',i,sep=''), get_data(parameter_space = '04', scenario = 'S', experiment = i, run = 1))
}


# Compare between PS ------------------------------------------------------
setwd('~/Documents/malaria_interventions_data/')
ps <- sprintf('%0.2d',27:39)[-11]
d <- map(1:5, function(r){
      map(ps, function(i){
        cat(i)
        tmp <- get_data(parameter_space = i, scenario = 'N', experiment = '002', run = r)
        return(tmp[[1]])
      }) %>% bind_rows()
    }) %>% bind_rows()


time_range <- c(28800,40000)

my_cols <- gg_color_hue(length(unique(d$PS)), hue_min = 10, hue_max = 280, l = 62, c = 200)
# my_cols <- c('red','green','blue','orange','black')
d %>%
  left_join(subset(design, select=c(PS,BITING_RATE_MEAN,N_GENES_INITIAL)), by='PS') %>% 
  mutate(N_GENES_INITIAL=as.factor(N_GENES_INITIAL)) %>% 
  select(-year, -month, -n_infected) %>% 
  filter(time>time_range[1]&time<time_range[2]) %>%
  gather(variable, value, -pop_id, -time, -exp, -PS, -scenario, -run,-N_GENES_INITIAL, -BITING_RATE_MEAN) %>% 
  group_by(pop_id, time, exp, PS, scenario, N_GENES_INITIAL, BITING_RATE_MEAN, variable) %>%
  summarise(value_mean=mean(value), value_sd=sd(value)) %>% 
  filter(variable %in% c('prevalence', 'meanMOI','n_circulating_strains', 'n_circulating_genes', 'n_alleles', 'n_infected_bites')) %>%
  # filter(variable %in% c('n_circulating_strains')) %>%
  ggplot(aes(x=time, y=value_mean, color=N_GENES_INITIAL))+
  geom_line()+
  geom_vline(xintercept = c(29160,29160+3600), linetype='dashed')+
  # geom_vline(xintercept = 29160, linetype='dashed')+
  scale_color_manual(values=my_cols)+
  scale_x_continuous(breaks=pretty(x=subset(d, time>time_range[1]&time<time_range[2])$time,n=5))+
  facet_wrap(~variable, scales='free')+
  mytheme


e <- sprintf('%0.2d', 6:10)
d <- map(e, function(i){
  cat(i)
  tmp <- get_data(parameter_space = i, scenario = 'S', experiment = '001', run = 1)
  return(tmp[[1]])
}) %>% bind_rows()
d %>% 
  ggplot(aes(x=month,y=EIR, color=PS, group=PS))+
  # geom_boxplot()+
  geom_point(stat='summary', fun.y=mean)+
  # scale_y_continuous(limits = c(0,10))+
  stat_summary(fun.y=mean, geom="line")+
  mytheme


# Compare between experiments ---------------------------------------------
control <- get_data(parameter_space = '39', scenario = 'S', experiment = '001', run = 1)[[1]]
control_means <- control %>% select(-year, -month, -n_infected) %>% 
  gather(variable, value, -time, -exp, -PS, -scenario, -run, -pop_id) %>% 
  group_by(PS, exp, variable) %>% summarise(mean_value=mean(value))
control %>% select(-year, -month, -n_infected) %>% 
  gather(variable, value, -time, -exp, -PS, -scenario, -run, -pop_id) %>% 
  ggplot(aes(x=time, y=value))+
  geom_line()+
  facet_wrap(~variable,scales='free')+
  geom_hline(aes(yintercept=mean_value), data=control_means, color='blue')+
  mytheme

e <- sprintf('%0.3d', 1:4)
d <- map(1:5, function(r){
  map(e, function(i){
    cat(i)
    tmp <- get_data(parameter_space = '39', scenario = 'S', experiment = i, run = r)
    return(tmp[[1]])
  }) %>% bind_rows()
}) %>% bind_rows()

# Add control to the design_irs data frame, which is created in build_parameter_files.R
design_irs <- create_intervention_scheme_IRS(PS_benchmark = '12', scenario_benchmark = 'S', IRS_START_TIMES = '29160', immigration_range=c(0,0.001), length_range=c(1800,3600), coverage_range=0.9, write_to_file = T, design_ref=design)
design_irs %<>% slice(rep(1, each = 1)) %>% bind_rows(design_irs)
design_irs[1, 'exp'] <- '001'
design_irs[1, grepl("IRS", names(design_irs))] <- 'control'
# design_irs[nrow(design_irs), 'IRS_coverage'] <- 'none'

# mintime=d %>% group_by(exp) %>% summarise(m=max(time)) %>% summarise(min(m))
# mintime=mintime[1,1]
# pdf('seasonal_comparison.pdf',16,10)
time_range <- c(28800,40000)

my_cols <- c('black',gg_color_hue(length(unique(design_irs$exp))-1, hue_min = 10, hue_max = 280, l = 62, c = 200))
d %>%
  select(-year, -month, -n_infected) %>% 
  filter(time>time_range[1]&time<time_range[2]) %>%
  gather(variable, value, -time, -exp, -PS, -scenario, -run) %>% 
  # filter(variable %in% c('prevalence', 'meanMOI','n_circulating_strains', 'n_circulating_genes', 'n_alleles', 'n_infected_bites')) %>%
  left_join(design_irs) %>%
  # filter(IRS_IMMIGRATION=='1' | IRS_IMMIGRATION=='control') %>%
  # filter(IRS_coverage=='0.9' | IRS_coverage=='control') %>%
  ggplot(aes(x=time, y=value, color=exp))+
  geom_vline(xintercept = c(29160,29160+3600), linetype='dashed')+
  # geom_vline(xintercept = seq(28800,39960,720),color='gray')+
  geom_line()+
  scale_color_manual(values=my_cols)+
  scale_x_continuous(breaks=pretty(x=subset(d, time>time_range[1]&time<time_range[2])$time,n=5))+
  geom_hline(aes(yintercept=mean_value), data=control_means, color='blue')+
  facet_wrap(~variable, scales = 'free')+mytheme

my_cols <- c(gg_color_hue(length(unique(design_irs$IRS_IMMIGRATION))-1, hue_min = 10, hue_max = 280, l = 62, c = 200),'black')
d %>%
  select(-year, -month, -n_infected) %>% 
  filter(time>time_range[1]&time<time_range[2]) %>%
  gather(variable, value, -time, -exp, -PS, -scenario, -run) %>% 
  filter(variable %in% c('prevalence', 'meanMOI','n_circulating_strains', 'n_circulating_genes')) %>%
  
  left_join(design_irs) %>% 
  filter(IRS_length=='7200' | IRS_length=='control') %>%
  filter(IRS_coverage=='0.9' | IRS_coverage=='control') %>%
  ggplot(aes(x=time, y=value, color=IRS_IMMIGRATION))+
  # geom_vline(xintercept = intervention_start+seq(0,7200,1800),color='gray')+
  # geom_vline(xintercept = seq(28800,39960,720),color='gray')+
  geom_line()+
  scale_color_manual(values=my_cols)+
  scale_x_continuous(breaks=pretty(x=subset(d, time>time_range[1]&time<time_range[2])$time,n=5))+
  facet_wrap(~variable, scales = 'free')+mytheme


d %>% 
  ggplot(aes(x=month,y=EIR, color=exp, group=exp))+
  # geom_boxplot()+
  geom_point(stat='summary', fun.y=mean)+
  # scale_y_continuous(limits = c(0,10))+
  stat_summary(fun.y=mean, geom="line")+
  mytheme
#dev.off()


# Calculate stats at the end of interventions -----------------------------

# These are some summary statistics that would quantify the effect of
# intervention on a particular metric (e.g., prevalence. numbner of circulating
# genes) after an intervention is lifted. The comparison is done to the control
# experiment. 


intervention_stats <- c()
intervention_stats_diff <- c()
PS <- sprintf('%0.2d', (27:39)[-11])
exp <- sprintf('%0.3d', 1:4)
design_irs <- create_intervention_scheme_IRS(PS_benchmark = '27', scenario_benchmark = 'S', IRS_START_TIMES = '29160', immigration_range=c(0), length_range=c(720,1800,3600), coverage_range=0.9, write_to_file = F, design_ref=design)
# Add control to the design_irs data frame, which is created in build_parameter_files.R
design_irs %<>% slice(rep(1, each = 1)) %>% bind_rows(design_irs)
design_irs[1, 'exp'] <- '001'
design_irs[1, 'IRS_length'] <- 0

for (run in 1:5){
  for (ps in PS){
    if (file.exists(paste('~/Documents/malaria_interventions_data/PS',ps,'_S_E001_R',run,'.sqlite',sep=''))){
      control_data <- get_data(parameter_space = ps, scenario = 'S', experiment = '001', run = run)[[1]]
    }
    for (e in exp){
      print(paste(Sys.time(),run,ps,e,sep=' | '))
      if (file.exists(paste('~/Documents/malaria_interventions_data/PS',ps,'_S_E',e,'_R',run,'.sqlite',sep=''))){
        x <- post_intervention_stats(PS = ps, scenario = 'S', exp=e, run = run, plot.it = F, control_data = control_data, design_irs = design_irs)
        intervention_stats <- rbind(intervention_stats, x$summary_stats)
        intervention_stats_diff <- rbind(intervention_stats_diff, x$diff_control)
      } else {
        cat('missing file: ');cat(paste('PS',ps,'_S_E',e,'_R',run,'.sqlite',sep=''));cat('\n')
      }
    }
  }
}

# Time to extinction as a function of the diversity
intervention_stats %>% left_join(subset(design, select=c(PS,BITING_RATE_MEAN,N_GENES_INITIAL)), by='PS') %>% 
  # filter(BITING_RATE_MEAN=='0.0005') %>% 
  group_by(PS, BITING_RATE_MEAN, N_GENES_INITIAL, scenario, exp) %>% 
  summarise(time_ext_max=max(time_extinct), time_ext_mean=mean(time_extinct)) %>% 
  ggplot()+
  geom_point(aes(x=PS, y=time_ext_max), size=3, color='red')+
  geom_point(aes(x=PS, y=time_ext_mean), size=3, color='blue')+
  labs(y='Mean and max time to extinction')+
  facet_grid(~exp)+
  mytheme

# Calculate probability of extinction as a proportion of runs which went extinct
intervention_stats %>% 
  left_join(subset(design, select=c(PS,BITING_RATE_MEAN,N_GENES_INITIAL)), by='PS') %>% 
  left_join(subset(design_irs, select=c(exp,IRS_START_TIMES,IRS_length))) %>%
  distinct(PS, exp, run, time_extinct,IRS_START_TIMES,IRS_length) %>% 
  mutate(extinct=ifelse(time_extinct<as.numeric(IRS_START_TIMES)+as.numeric(IRS_length),1,0)) %>%
  group_by(PS, exp) %>% 
  summarise(extinct_prob=sum(extinct)/max(intervention_stats$run)) %>% 
  ggplot(aes(x=PS, y=extinct_prob, group=exp, color=exp))+
  geom_point(size=4)+
  geom_line()+
  facet_wrap(~exp)+mytheme

         
# This calculates the ratio of the value of variables POST-intervention compared
# to control. For example, we can say that intervention has reduced the
# prevalence by 50% compared to control in a given PS and experiment.
intervention_stats_diff %>% 
  left_join(subset(design, select=c(PS,BITING_RATE_MEAN,N_GENES_INITIAL)), by='PS') %>% 
  left_join(subset(design_irs, select=c(exp,IRS_START_TIMES,IRS_length))) %>%
  mutate(time_threshold=as.numeric(IRS_START_TIMES)+as.numeric(IRS_length)+360) %>%
  filter(time>=time_threshold) %>% 
  
  filter(variable %in% c('n_alleles','n_circulating_genes','n_circulating_strains','prevalence','n_total_bites','meanMOI')) %>% 
  
  # select(variable, mean_value, PS, exp, run) %>% arrange(variable, PS, run, exp)
  group_by(PS, BITING_RATE_MEAN, N_GENES_INITIAL, scenario, exp, variable) %>% 
  summarise(mean_diff=mean(diff), 
            mean_ratio=mean(ratio),
            sd_ratio=sd(ratio)) %>%
  mutate(error_up=mean_ratio+sd_ratio,
         error_low=mean_ratio-sd_ratio) %>% 
  ggplot(aes(x=PS, y=mean_ratio, color=exp, group=exp))+
  geom_point(size=3)+
  geom_line()+
  # geom_errorbar(aes(ymin=error_low, ymax=error_up), width=0.1)+
  # geom_line(aes(x=PS, y=error_up, color=exp, group=exp), linetype='dashed')+
  # geom_line(aes(x=PS, y=error_low, color=exp, group=exp), linetype='dashed')+
  facet_wrap(BITING_RATE_MEAN~variable,scales='free_y')+
  mytheme

# Time series of difference from control (i.e., experiment-control)
intervention_stats_diff %>%
  left_join(subset(design, select=c(PS,BITING_RATE_MEAN,N_GENES_INITIAL)), by='PS') %>% 
  filter(variable %in% c('n_circulating_strains','prevalence')) %>% 
  # filter(exp=='004') %>% filter(PS=='27') %>% 
  group_by(time,variable,PS,scenario,exp,BITING_RATE_MEAN,N_GENES_INITIAL) %>% 
  summarise(diff=mean(diff)) %>% 
  ggplot(aes(x=time, y=diff, color=as.factor(N_GENES_INITIAL)))+
  scale_color_manual(values = gg_color_hue(length(PS),hue_min = 10, 400))+
  geom_line()+
  geom_hline(yintercept = 0, color='black')+
  facet_wrap(variable~BITING_RATE_MEAN+exp, scales='free_y')+
  mytheme




# Also can look at a bulk:
# First get all data across PS, experiments and runs
setwd('~/Documents/malaria_interventions_data/')
PS <- sprintf('%0.2d', 6:12)
exp <- sprintf('%0.3d', 1:4)
d <- c()
for (ps in PS){
  for (e in exp){
    cat(ps);cat('--');cat(e);cat('\n')
    tmp <- get_data(parameter_space = ps, scenario = 'S', experiment = e, run = 1)
    d <- rbind(d, tmp$summary_general)
  }
}
d <- as.tibble(d)
x <- d %>% select(PS, exp, time, prevalence) %>% 
  left_join(subset(design_irs, select=c(IRS_START_TIMES,IRS_length,exp)), by='exp') %>%
  select(PS, exp, time, prevalence, IRS_START_TIMES, IRS_length) %>% 
  mutate(intervention_lift=as.numeric(IRS_START_TIMES)+as.numeric(IRS_length)) %>% 
  filter(time>intervention_lift) %>% group_by(PS, exp) %>% summarise(max_prev=max(prevalence))

d %>% left_join(subset(design_irs, select=c(IRS_START_TIMES,IRS_length,exp)), by='exp') %>%
  mutate(intervention_lift=as.numeric(IRS_START_TIMES)+as.numeric(IRS_length)) %>% 
  filter(time>intervention_lift) %>% 
  select(-year, -month, -n_infected) %>% 
  gather(variable, value, -time, -exp, -PS, -scenario, -run,-IRS_START_TIMES,-IRS_length) %>% 
  group_by(PS, exp, variable) %>% summarise(max_value=max(value)) %>% 
  ggplot(aes(x=PS, y=max_value))+
  geom_point()+
  facet_grid(variable~exp, scales='free_y')+mytheme


sample_times <- c(intervention_start+seq(0,7200,1800))
length_range # This comes from build_parameter_files.R

x <- map(sample_times, function(s){
  d %>% left_join(design_irs) %>% 
    filter(IRS_length!='control') %>% 
    filter(time %in% s) %>% 
    filter(IRS_coverage==0.9)
}) %>% bind_rows()

x %>% ggplot(aes(x=IRS_IMMIGRATION, y=n_infections, color=as.numeric(IRS_IMMIGRATION)))+geom_point()+facet_wrap(~IRS_length)



# Compare between scenarios ---------------------------------------------
PS04_S_001 <- get_data(parameter_space = '04', scenario = 'S', experiment = '001', run = 1)
PS04_N_001 <- get_data(parameter_space = '04', scenario = 'N', experiment = '001', run = 1)

d <- rbind(PS04_S_001[[1]],
           PS04_N_001[[1]])

# mintime=d %>% group_by(exp) %>% summarise(m=max(time)) %>% summarise(min(m))
# mintime=mintime[1,1]
# pdf('seasonal_comparison.pdf',16,10)
time_range <- c(28000,36000)
# png('scenario_comparison_1.png',1800,1000)
d %>%
  select(-year, -month, -n_infected) %>% 
  filter(time>time_range[1]&time<time_range[2]) %>%
  gather(variable, value, -time, -exp, -PS, -scenario, -run) %>% 
  ggplot(aes(x=time, y=value, color=scenario, group=scenario))+
  geom_line()+
  # geom_vline(xintercept = c(21600,21960,22320,22680,23040,23400))+
  scale_x_continuous(breaks=pretty(x=subset(d, time>time_range[1]&time<time_range[2])$time,n=5))+
  scale_color_manual(values=c('blue','orange','red'))+
  facet_wrap(~variable, scales = 'free')+mytheme
# dev.off()


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




# Structure ---------------------------------------------------------------
# Initialize
parameter_space <- '03'
scenario <- 'S'
experiment <- '01'
run <- 1
base_name <- paste('PS',parameter_space,'_',scenario,'_E',experiment,'_R',run,sep='')
sqlite_file <- paste(base_name,'.sqlite',sep='')
# parameter_file <- paste(base_name,'.py',sep='') # This may be necessary so I leave it

# Extract data from sqlite. variable names correspond to table names
db <- dbConnect(SQLite(), dbname = sqlite_file)
sampled_strains <- as.tibble(dbGetQuery(db, 'SELECT id, gene_id FROM sampled_strains'))
names(sampled_strains)[1] <- c('strain_id')
sampled_alleles <- as.tibble(dbGetQuery(db, 'SELECT * FROM sampled_alleles'))
names(sampled_alleles)[3] <- c('allele_id')
sampled_strains <- full_join(sampled_strains, sampled_alleles)
sampled_strains$allele_locus <- paste(sampled_strains$allele_id,sampled_strains$locus,sep='_') # each allele in a locus is unique

sampled_infections <- PS03_S_01[[2]]
sampled_infections$layer <- group_indices(sampled_infections, time) # Add layer numbers. Each time point is a layer
sampled_infections %<>% arrange(layer,strain_id,host_id) %>% group_by(layer,strain_id) %>% mutate(strain_copy = row_number()) # add unique id for each strain copy within each layer. A strain copy is an instance of a strain in a particualr host
sampled_infections$strain_id_unique=paste(sampled_infections$strain_id,sampled_infections$strain_copy,sep='_') # Create the unique strains

# Integrate the strain composition into the infections table
identical(unique(sampled_infections$strain_id),unique(sampled_strains$strain_id)) # This pretty important
sampled_infections %<>% select(time, layer, host_id, strain_id, strain_copy, strain_id_unique) %>% left_join(sampled_strains)

#sampled_infections <- subset(sampled_infections, layer<=10)


overlapAlleleAdj<-function(mat){
  newmat<-tcrossprod(mat>0)
  newmat<-newmat/rowSums(mat>0)
  return(newmat)  
}

# A function to build the similarity matrix for a single layer and calculate some summary stats
build_layer <- function(l){
  sampled_infections_layer <- subset(sampled_infections, layer==l) # This is MUCH faster than sampled_infections_layer <- sampled_infections %>% filter(layer==l)
  sampled_infections_layer %<>% group_by(strain_id) %>% mutate(freq = n()/120) %>% arrange(strain_id_unique) # strain frequency (the number of strain copies should be equal to the frequency)
  # Calculate the edges
  similarity_matrix <- table(sampled_infections_layer$strain_id_unique, sampled_infections_layer$allele_locus)
  similarity_matrix <- overlapAlleleAdj(similarity_matrix)
  # Some summary
  layer_summary <- with(sampled_infections_layer, 
                        data.frame(hosts=length(unique(host_id)),
                                   strains=length(unique(strain_id)),
                                   total_infections=length(unique(strain_id_unique))
                                   ))
  return(list(similarity_matrix=similarity_matrix, infections=sampled_infections_layer, layer_summary=layer_summary))
}

Layers <- list()
for (l in 1:max(sampled_infections$layer)){
  cat(paste('[',Sys.time(), '] building layer ',l,'\n',sep=''))
  Layers[[l]] <- build_layer(l)  
}

# Get just the matrices
temporal_network <- unique(lapply(Layers, function(x) head(x, 1)$similarity_matrix))
sapply(temporal_network, dim)

# Apply a cutoff
x <- xtabs(~strain_id+allele_locus, subset(sampled_strains, strain_id%in%sampled_infections$strain_id))
similarityMatrix <- overlapAlleleAdj(x)
dim(similarityMatrix)
# # Write similarity matrix without cutoff. This used to plot the edge weight distributions.
# write.table(similarityMatrix, paste('../',filenameBase,'_similarityMatrix_nocutoff.csv',sep=''), row.names = T, col.names = T, sep=',')
cutoffValue <- quantile(as.vector(similarityMatrix), probs = 0.9)
print(cutoffValue)
# hist(as.vector(similarityMatrix));abline(v=cutoffValue,col='red')
similarityMatrix[similarityMatrix<cutoffValue] <- 0
# write.table(similarityMatrix, paste(filenameBase,'_similarityMatrix.csv',sep=''), row.names = T, col.names = T, sep=',')
# write.table(similarityMatrix, paste('../',filenameBase,'_similarityMatrix.csv',sep=''), row.names = T, col.names = T, sep=',')

for (i in 1:length(temporal_network)){
  print(i)
  x <- temporal_network[[i]]
  x[x<cutoffValue] <- 0
  temporal_network[[i]] <- x
}

# This function takes a list of temporal matrices and returns an intralayer and
# interlayer edge lists in a format: [layer_source, node_source, layer_target node_target, weight].
# It also returns the list of node names
build_infomap_objects <- function(temporal_network, write_to_infomap_file=T, return_objects=T){
  require(igraph)
  
  # Define functions
  
  ## A function to get the repertoire names from unique repertoire copy names
  splitText <- function(str,after=T,splitchar='\\.'){
    if (after){
      return(sapply(strsplit(str, split=splitchar), tail, 1))
    }
    if (after==F){
      return(sapply(strsplit(str, split=splitchar), head, 1))
    }
  }
  
  ##  A function that gets the layer as a matrix and writes it for infomap as an edge list
  matrix_to_infomap <- function(l){
    require(igraph)
    current_layer <- temporal_network[[l]]
    if(nrow(current_layer)<2){
      print(paste('No strains in layer',l,'!!! skipping intralayer edges'))
      next
    }
    g <- graph.adjacency(current_layer, mode = 'directed', weighted = T)
    current_layer_el <- as.tibble(as_data_frame(g, what = 'edges'))
    names(current_layer_el) <- c('node_s','node_t','w')
    current_layer_el$layer_s <- l
    current_layer_el$layer_t <- l
    current_layer_el$node_s <- nodeList$nodeID[match(current_layer_el$node_s,nodeList$nodeLabel)]
    current_layer_el$node_t <- nodeList$nodeID[match(current_layer_el$node_t,nodeList$nodeLabel)]
    current_layer_el %<>% select(layer_s, node_s, layer_t, node_t, w)
    print(paste('[',Sys.time(), '] Created edge list of layer ',l,' for Infomap | ', nrow(current_layer_el),' edges',sep=''))
    return(current_layer_el)
  }
  
  ## A function that calculates interlayer edges between layer t and t+1
  build_interlayer_edges_1step <- function(t){
    strain_copies_t <- rownames(temporal_network[[t]]) # repertoires at time t
    strain_copies_t1 <- rownames(temporal_network[[t+1]]) # repertoires at time t+1
    if (length(strain_copies_t)<2 | length(strain_copies_t1)<2){
      print(paste('No interlayer edges between layers',t, 'and',t+1,'because there are < 2 repertoires.'))
      next} # need minimum of 2 strains in t and t+1 to build a matrix
    # Pull the similarity values between the repertoires from the general similarity matrix, then write back the repertoire copy names
    inter_layer_edges_matrix <- similarityMatrix[splitText(strain_copies_t, after = F, '_'),splitText(strain_copies_t1, after = F, '_')] 
    rownames(inter_layer_edges_matrix) <- strain_copies_t
    colnames(inter_layer_edges_matrix) <- strain_copies_t1
    if (all(inter_layer_edges_matrix==0)){
      print(paste('All interlayer edges between layers',t, 'and',t+1,' are 0 (due to the cutoff). skipping.'))
      next
    }
    # transform to an edge list
    g <- graph.incidence(inter_layer_edges_matrix, directed = T, mode = 'out', weighted = T)
    inter_layer_edges <- as.tibble(as_data_frame(g, what = 'edges'))
    names(inter_layer_edges) <- c('node_s','node_t','w')
    inter_layer_edges$layer_s <- t
    inter_layer_edges$layer_t <- t+1
    inter_layer_edges$node_s <- nodeList$nodeID[match(inter_layer_edges$node_s,nodeList$nodeLabel)]
    inter_layer_edges$node_t <- nodeList$nodeID[match(inter_layer_edges$node_t,nodeList$nodeLabel)]
    inter_layer_edges %<>% select(layer_s, node_s, layer_t, node_t, w)
    print(paste('Built interlayer edges for layers: ',t,'-->',t+1,' | ', nrow(inter_layer_edges),' edges',sep=''))
    return(inter_layer_edges)
  }
  
  # Get the node list
  nodeLabel <- sort(unique(unlist(lapply(temporal_network,rownames))))
  nodeList <- data.frame(nodeID=1:length(nodeLabel), nodeLabel)
  
  # Create intralayer edge lists
  layers <- 1:length(temporal_network)
  infomap_intralayer <- map(layers, matrix_to_infomap) %>% bind_rows()
  
  # create interlayer edge lists
  layers <- layers[-length(layers)]
  infomap_interlayer <- map(layers, build_interlayer_edges_1step) %>% bind_rows()
  
  if (write_to_infomap_file){
    sink.reset <- function(){
      for(i in seq_len(sink.number())){
        sink(NULL)
      }
    }
    ## Write file for infomap
    print('Writing Infomap files')
    file <- paste(base_name,'_Infomap_multilayer','.txt',sep='')
    print(paste('Infomap file:',file))
    if (file.exists(file)){unlink(file)}
    sink(file, append = T)
    cat("# A network in a general multiplex format");cat('\n')
    cat(paste("*Vertices",nrow(nodeList)));cat('\n')
    write.table(nodeList, file, append = T,sep=' ', quote = T, row.names = F, col.names = F)
    cat("*Multiplex");cat('\n')
    cat("# layer node layer node [weight]");cat('\n')
    cat("# Intralayer edges");cat('\n')
    write.table(infomap_intralayer, file, sep = ' ', row.names = F, col.names = F, quote = F, append = T)
    cat("# Interlayer edges");cat('\n')
    write.table(infomap_interlayer, file, sep = ' ', row.names = F, col.names = F, quote = F, append = T)
    sink.reset()
  }
  
  if (return_objects){
    return(list(infomap_intralayer=infomap_intralayer, infomap_interlayer=infomap_interlayer, nodeList=nodeList))
  }
}

infomap_objects <- build_infomap_objects(temporal_network)



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
calendar %>% filter(running_day==10930)
