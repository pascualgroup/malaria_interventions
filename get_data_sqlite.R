library(tidyverse)
library(magrittr)
library(sqldf)

# Functions ---------------------------------------------------------------

## @knitr FUNCTIONS
loadExperiments_GoogleSheets <- function(workBookName='malaria_interventions_design',sheetID=4){
  require(googlesheets)
  GS <- gs_title(workBookName)
  col_types <- GS %>% gs_read(ws=1, col_names=T)
  col_t <- unname(as.list(col_types[1,]))
  experiments <- GS %>% gs_read(ws=sheetID, col_names=T, col_types=col_t)
  print(experiments)
  return(experiments)
}


create_intervention_scheme_IRS <- function(PS_benchmark, scenario_benchmark, IRS_START_TIMES=NULL, immigration_range, length_range, coverage_range, write_to_file=T, design_ref=NULL){
  # This fixes the 1. in the file name produced in Mathematica
  files <- list.files(path = '~/Documents/malaria_interventions_data/', full.names = T, pattern = '1\\._')
  if(length(files)>0){sapply(files, function(x){file.rename(from = x, to = str_replace(x, '\\.', ''))})}
  
  # PS_benchmark is the parameter space for which intervention is tested. it should already have the checkpoint (000) and control (001) experiments in the google sheet.
  # IRS_START_TIMES can also be something like: '29000,35000'
  if(is.null(design_ref)){
    design_ref <- loadExperiments_GoogleSheets() # Get data design
  }
  design_ref <- subset(design_ref, PS==PS_benchmark & scenario==scenario_benchmark)
  reference_row <- which(grepl(PS_benchmark, design_ref$PS))[2] # reference_row is the row in the design data set which is the reference for this set of IRS exepriments. Should be the line of the control experiment (not the checkpoint)
  
  scheme <- expand.grid(PS=design_ref$PS[reference_row],
                        scenario=design_ref$scenario[reference_row],
                        BITING_RATE_MEAN=design_ref$BITING_RATE_MEAN[reference_row],
                        DAILY_BITING_RATE_DISTRIBUTION=design_ref$DAILY_BITING_RATE_DISTRIBUTION[reference_row],
                        N_GENES_INITIAL=design_ref$N_GENES_INITIAL[reference_row],
                        N_LOCI=design_ref$N_LOCI[reference_row],
                        N_ALLELES_INITIAL=design_ref$N_ALLELES_INITIAL[reference_row],
                        N_POPULATIONS=design_ref$N_POPULATIONS[reference_row],
                        populations_dist=design_ref$populations_dist[reference_row],
                        T_BURNIN=design_ref$T_BURNIN[reference_row],
                        T_END=design_ref$T_END[reference_row],
                        IRS_START_TIMES=IRS_START_TIMES,
                        IRS_IMMIGRATION=immigration_range,
                        IRS_length=length_range,
                        IRS_coverage=coverage_range,
                        MDA_START=NA, # This so to not have an MDA (for function generate_files)
                        wall_time=design_ref$wall_time[reference_row],
                        mem_per_cpu=design_ref$mem_per_cpu[reference_row],
                        stringsAsFactors = F)
  scheme$exp <- sprintf('%0.3d', 2:(nrow(scheme)+1)) # Start with 002 because 000 and 001 are checkpoint and control
  # Important!!! The files for the IRS in the next line are generated separately in Mathematica, and should already exist.
  scheme$IRS_input <- paste('IRS_120_300_',scheme$IRS_coverage,'_',scheme$IRS_length,'.csv',sep='')
  if (write_to_file){
    write_csv(scheme, paste('PS',scheme$PS[1],'_',scheme$scenario[1],'_IRS_scheme.csv',sep=''))
  }
  return(scheme)
}


mytheme <- theme_bw() + theme(
  legend.title  = element_text(colour = "black", size=17),
  # legend.position = "none",
  #	legend.direction = "horizontal",
  legend.key = element_blank(),
  legend.text  = element_text(colour = "black", size=17),
  panel.background = element_blank(),
  panel.grid.major = element_blank(),
  # panel.grid.minor = element_blank(),
  axis.text = element_text(color='black', family="Helvetica", size=14),
  axis.title = element_text(color='black', family="Helvetica", size=14),
  strip.text.x = element_text(family = "Helvetica", size = 14),
  strip.text.y = element_text(family = "Helvetica", size = 14),
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

get_data <- function(parameter_space, scenario, experiment, run, sampling_period=30, host_age_structure=F){
  require(sqldf)
  # Initialize
  base_name <- paste('PS',parameter_space,'_',scenario,'_E',experiment,'_R',run,sep='')
  sqlite_file <- paste(base_name,'.sqlite',sep='')
  if (!file.exists(sqlite_file)) {
    print (paste(sqlite_file, ' does not exist, ignoring and terminating function'))
    return(NULL)
  }
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
  
  summary_general <- suppressMessages(left_join(summary_general, summary_alleles))
  
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
  
  ## Get EIR from the table in the sqlite
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
  summary_general <- suppressMessages(inner_join(summary_general, meanMOI)) # Add MOI to the summary
  
  
  # Host age structure
  if(host_age_structure) {
    hosts <- dbGetQuery(db, 'SELECT * FROM hosts')
    names(hosts)[1] <- 'host_id'
    hosts$lifespan <- round((hosts$death_time-hosts$birth_time))
    hosts <- subset(hosts, host_id%in%sampled_infections$host_id)
    sampled_infections <- suppressMessages(left_join(sampled_infections, hosts, by='host_id'))
    sampled_infections$host_age <- round((sampled_infections$time-sampled_infections$birth_time)/30)    
  }

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
  diff_control <- suppressMessages(inner_join(x,y,by=c('pop_id','PS','scenario','run','time','variable'))) %>% 
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
  
  summary_stats <- suppressMessages(left_join(time_extinct, new_mean))
  summary_stats <- suppressMessages(left_join(summary_stats, amplitude))
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

## @knitr END

# Some example to test ----------------------------------------------------

# Join pre-intervention (E000) and intervention (E003) time-series
ctrl <- get_data(parameter_space = '36', scenario = 'S', experiment = '000', run = 1)[[1]]
x <- get_data(parameter_space = '36', scenario = 'S', experiment = '003', run = 1)[[1]]
# Plot
svg('/home/shai/Google Drive/LabSync/FacultyJobs/UToronto 2018/time_series_intervention.svg', 2.4, 1.7)
x %>% bind_rows(subset(ctrl, time<min(x$time))) %>% 
 filter(time>28000 & time<35000) %>% 
  ggplot(aes(x=time, y=n_circulating_strains))+
  geom_line(color='#900C3F', size=1)+
  labs(x='Time',y='Diversity')+
  mytheme+
    theme(axis.text=element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank())
dev.off()


## @knitr INITIALIZE

# Initialize important variables ------------------------------------------
setwd('~/Documents/malaria_interventions_data/')
design <- loadExperiments_GoogleSheets() 
ps_range <- sprintf('%0.2d', 27:39)
exp_range <- sprintf('%0.3d', 1:4)
run_range <- 1:10
scenario <- 'S'
monitored_variables <- c('prevalence', 'meanMOI','n_circulating_strains', 'n_circulating_genes', 'n_alleles', 'n_total_bites')
exp_cols <- c('black','#0A97B7','#B70A97','#97B70A')
scenario_cols <- c('red','blue','orange')

## @knitr COMPARE_EXPERIMENTS_LOAD

# Compare between experiments within a parameter space --------------------
PS <- '36'
control <- get_data(parameter_space = PS, scenario = scenario, experiment = '001', run = 1)[[1]]
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

exp_comparison <- map(run_range, function(r){
  map(exp_range, function(e){
    print(paste(e,r,sep=' | '))
    tmp <- get_data(parameter_space = PS, scenario = scenario, experiment = e, run = r)
    return(tmp[[1]])
  }) %>% bind_rows()
}) %>% bind_rows()

# Add control to the design_irs data frame, which is created in build_parameter_files.R
design_irs <- create_intervention_scheme_IRS(PS_benchmark = PS, scenario_benchmark = scenario,IRS_START_TIMES = '29160', immigration_range=c(0), length_range=c(720,1800,3600), coverage_range=0.9, write_to_file = F, design_ref=design)
design_irs %<>% slice(rep(1, each = 1)) %>% bind_rows(design_irs)
design_irs[1, 'exp'] <- '001'
design_irs[1, grepl("IRS", names(design_irs))] <- 'control'

time_range <- c(28800,36000)

## @knitr COMPARE_EXPERIMENTS_PLOT

exp_comparison %>%
  select(-year, -month, -n_infected) %>% 
  filter(time>time_range[1]&time<time_range[2]) %>%
  gather(variable, value, -time, -exp, -PS, -scenario, -run, -pop_id) %>% 
  filter(variable %in% monitored_variables) %>%
  group_by(time, PS, exp, variable) %>% # Need to average across runs
  summarise(value_mean=mean(value)) %>% 
  ggplot(aes(x=time, y=value_mean, color=exp))+
  geom_line()+
  facet_wrap(~variable, scales = 'free')+
  geom_vline(xintercept = c(29160+c(0,720,2800,3600)), linetype='dashed')+
  scale_color_manual(values=exp_cols)+
  # scale_x_continuous(breaks=pretty(x=subset(d, time>time_range[1]&time<time_range[2])$time,n=5))+
  geom_hline(aes(yintercept=mean_value), data=subset(control_means, variable %in% monitored_variables) , color='blue')+
  mytheme


## @knitr COMPARE_DIVERSITY_LOAD

# Compare between parameter spaces within an experiment -------------------
exp <- '003'
ps_comparison <- map(run_range, function(r){
      map(ps_range, function(ps){
        print(paste(ps,r,sep=' | '))
        tmp <- get_data(parameter_space = ps, scenario = scenario, experiment = exp, run = r)
        return(tmp[[1]])
      }) %>% bind_rows()
    }) %>% bind_rows()

time_range <- c(28800,max(ps_comparison$time))
my_cols <- gg_color_hue(length(unique(ps_comparison$PS)), hue_min = 10, hue_max = 280, l = 62, c = 200)

## @knitr COMPARE_DIVERSITY_PLOT
ps_comparison %>%
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
  geom_vline(xintercept = c(29160,29160+5*360), linetype='dashed')+
  scale_color_manual(values=my_cols)+
  scale_x_continuous(breaks=pretty(x=subset(ps_comparison, time>time_range[1]&time<time_range[2])$time,n=5))+
  facet_wrap(~variable, scales='free')+
  mytheme
# +theme(legend.position = "none",
#                 axis.text = element_text(color='black', family="Helvetica", size=5),
#   strip.text.x = element_text(family = "Helvetica", size = 5),
#   strip.text.y = element_text(family = "Helvetica", size = 5)
#                 )

## @knitr END

# Compare EIR
ps_comparison_eir <- map(run_range, function(r){
  map(sprintf('%0.2d', c(27,30,36,39)), function(ps){
    tmp <- get_data(parameter_space = ps, scenario = scenario, experiment = '001', run = r)
    return(tmp[[1]])
  }) %>% bind_rows()
}) %>% bind_rows()

ps_comparison_eir %>% 
  left_join(design, by='PS') %>% 
  mutate(N_GENES_INITIAL=as.factor(N_GENES_INITIAL)) %>% 
  ggplot(aes(x=month, y=EIR, color=N_GENES_INITIAL))+
  geom_boxplot()+
  stat_summary(aes(group=N_GENES_INITIAL), fun.y=mean, geom="line", size=1)+
  facet_wrap(~N_GENES_INITIAL)+
  scale_y_continuous(breaks=seq(0,20,2))+
  mytheme
  

## @knitr INTERVENTION_STATS_LOAD

# Calculate stats at the end of interventions -----------------------------

# These are some summary statistics that would quantify the effect of
# intervention on a particular metric (e.g., prevalence. numbner of circulating
# genes) after an intervention is lifted. The comparison is done to the control
# experiment. 

intervention_stats <- c()
intervention_stats_diff <- c()
design_irs <- create_intervention_scheme_IRS(PS_benchmark = '27', scenario_benchmark = scenario, IRS_START_TIMES = '29160', immigration_range=c(0), length_range=c(720,1800,3600), coverage_range=0.9, write_to_file = F, design_ref=design)
# Add control to the design_irs data frame, which is created in build_parameter_files.R
design_irs %<>% slice(rep(1, each = 1)) %>% bind_rows(design_irs)
design_irs[1, 'exp'] <- '001'
design_irs[1, 'IRS_length'] <- 0

for (run in run_range){
  for (ps in ps_range){
    if (file.exists(paste('~/Documents/malaria_interventions_data/PS',ps,'_',scenario,'_E001_R',run,'.sqlite',sep=''))){
      control_data <- get_data(parameter_space = ps, scenario = scenario, experiment = '001', run = run)[[1]]
    }
    for (e in exp_range){
      print(paste(Sys.time(),run,ps,e,sep=' | '))
      if (file.exists(paste('~/Documents/malaria_interventions_data/PS',ps,'_',scenario,'_E',e,'_R',run,'.sqlite',sep=''))){
        x <- post_intervention_stats(PS = ps, scenario = scenario, exp=e, run = run, plot.it = F, control_data = control_data, design_irs = design_irs)
        intervention_stats <- rbind(intervention_stats, x$summary_stats)
        intervention_stats_diff <- rbind(intervention_stats_diff, x$diff_control)
      } else {
        cat('missing file: ');cat(paste('PS',ps,'_',scenario,'_E',e,'_R',run,'.sqlite',sep=''));cat('\n')
      }
    }
  }
}


## @knitr TIME_TO_EXTINCTION_PLOT

# Time to extinction as a function of the diversity
intervention_stats %>% 
  left_join(subset(design, select=c(PS,BITING_RATE_MEAN,N_GENES_INITIAL)), by='PS') %>% 
  group_by(PS, BITING_RATE_MEAN, N_GENES_INITIAL, scenario, exp) %>% 
  summarise(time_ext_max=max(time_extinct), time_ext_mean=mean(time_extinct)) %>% 
  ggplot()+
  geom_point(aes(x=as.numeric(N_GENES_INITIAL), y=time_ext_max), size=5, color='red')+
  geom_point(aes(x=as.numeric(N_GENES_INITIAL), y=time_ext_mean), size=5, color='blue')+
  geom_line(aes(x=as.numeric(N_GENES_INITIAL), y=time_ext_max), color='red')+
  geom_line(aes(x=as.numeric(N_GENES_INITIAL), y=time_ext_mean), color='blue')+
  scale_x_continuous("Gene pool size", breaks = seq(1200,15600,1200)) +
  labs(y='Mean and max time to extinction')+
  facet_grid(~exp)+
  mytheme+
  theme(axis.text.x = element_text(angle = 90, hjust=0))


## @knitr PROB_OF_EXTINCTION_PLOT

# Calculate probability of extinction as a proportion of runs which went extinct
intervention_stats %>% 
  left_join(subset(design, select=c(PS,BITING_RATE_MEAN,N_GENES_INITIAL)), by='PS') %>% 
  left_join(subset(design_irs, select=c(exp,IRS_START_TIMES,IRS_length))) %>%
  distinct(PS, N_GENES_INITIAL, exp, run, time_extinct,IRS_START_TIMES,IRS_length) %>% 
  mutate(extinct=ifelse(time_extinct<as.numeric(IRS_START_TIMES)+as.numeric(IRS_length),1,0)) %>%
  group_by(PS, N_GENES_INITIAL, exp) %>% 
  summarise(extinct_prob=sum(extinct)/max(intervention_stats$run)) %>% 
  ggplot(aes(x=N_GENES_INITIAL, y=extinct_prob, group=exp, color=exp))+
  geom_point(size=4)+
  geom_line()+
  scale_x_continuous("Gene pool size", breaks = seq(1200,15600,1200)) +
  scale_color_manual(values=exp_cols)+
  facet_wrap(~exp)+
  mytheme+
  theme(axis.text.x = element_text(angle = 90, hjust=0))

## @knitr IRS_EFFECT_TS_PLOT

# A time series of difference between control and experiment.
intervention_stats_diff %>%
  left_join(subset(design, select=c(PS,BITING_RATE_MEAN,N_GENES_INITIAL)), by='PS') %>% 
  filter(variable %in% c('n_circulating_strains','prevalence')) %>% 
  filter(exp=='003') %>% filter(PS=='36') %>%
  group_by(time,variable,PS,scenario,exp,BITING_RATE_MEAN,N_GENES_INITIAL) %>% 
  summarise(diff=mean(diff)) %>% 
  ggplot(aes(x=time, y=diff))+
  geom_line(color='#97B70A', size=1)+
  geom_hline(yintercept = 0, color='black')+
  facet_wrap(~variable, scales='free_y')+
  mytheme+theme(legend.position = 'none')


## @knitr IRS_EFFECT_RATIO_PLOT

# This calculates the ratio of the value of variables POST-intervention compared
# to control. For example, we can say that intervention has reduced the
# prevalence by 50% compared to control in a given PS and experiment.
intervention_stats_diff %>%
  left_join(subset(design, select=c(PS,BITING_RATE_MEAN,N_GENES_INITIAL)), by='PS') %>%
  left_join(subset(design_irs, select=c(exp,IRS_START_TIMES,IRS_length))) %>%
  mutate(time_threshold=as.numeric(IRS_START_TIMES)+as.numeric(IRS_length)+360) %>%
  filter(time>=time_threshold) %>%
  
  filter(variable %in% monitored_variables) %>%
  
  # select(variable, mean_value, PS, exp, run) %>% arrange(variable, PS, run, exp)
  group_by(PS, BITING_RATE_MEAN, N_GENES_INITIAL, scenario, exp, variable) %>%
  # filter(exp!='001') %>%
  summarise(mean_diff=mean(diff),
            mean_ratio=mean(ratio),
            sd_ratio=sd(ratio)) %>%
  mutate(error_up=mean_ratio+sd_ratio,
         error_low=mean_ratio-sd_ratio) %>%
  ggplot(aes(x=as.numeric(N_GENES_INITIAL), y=mean_ratio, color=exp, group=exp))+
  geom_point(size=3)+
  geom_line()+
  scale_x_continuous("Gene pool size", breaks = seq(1200,15600,1200)) +
  geom_errorbar(aes(ymin=error_low, ymax=error_up), width=0.1)+
  geom_line(aes(x=as.numeric(N_GENES_INITIAL), y=error_up, color=exp, group=exp), linetype='dashed')+
  geom_line(aes(x=as.numeric(N_GENES_INITIAL), y=error_low, color=exp, group=exp), linetype='dashed')+
  scale_color_manual(values=exp_cols)+
  facet_wrap(~variable,scales='free_y')+
  mytheme+
  theme(axis.text.x = element_text(angle = 90, hjust=0))


## @knitr END


# Compare between scenarios within a parameter space and experiment -------

## @knitr COMPARE_SCENARIOS_TS_LOAD

PS <- '36'
cases <- expand.grid(scenario=c('S','G','N'), exp=sprintf('%0.3d',1:4), run=run_range)
scenario_comparison <- c()
for (i in 1:nrow(cases)){
  print(paste('Scenario: ',cases$scenario[i],' | exp: ',cases$exp[i], ' | run: ',cases$run[i],sep=''))
  tmp <- get_data(parameter_space = PS, scenario = cases$scenario[i], experiment = cases$exp[i], run = cases$run[i])[[1]]
  scenario_comparison <- rbind(scenario_comparison, tmp)
}

time_range <- c(28800,max(scenario_comparison$time))

## @knitr COMPARE_SCENARIOS_TS_PLOT

scenario_comparison %>%
  mutate(scenario=factor(scenario, levels=c('S','N','G'))) %>% 
  select(-year, -month, -n_infected) %>% 
  filter(time>time_range[1]&time<time_range[2]) %>%
  gather(variable, value, -time, -exp, -PS, -scenario, -run) %>% 
  group_by(time, exp, scenario, variable) %>%
  summarise(value_mean=mean(value), value_sd=sd(value)) %>% # Need to average across runs
  filter(variable %in% c('prevalence','n_circulating_strains', 'n_circulating_genes', 'n_alleles')) %>%
  ggplot(aes(x=time, y=value_mean, color=scenario))+
  geom_line()+
  scale_color_manual(values=scenario_cols)+
  facet_grid(variable~exp, scales='free')+
  mytheme


## @knitr COMPARE_IRS_EFFECT_SCENARIOS_LOAD

# Compare the effects of intervention with summary stats

# These are some summary statistics that would quantify the effect of
# intervention on a particular metric (e.g., prevalence. numbner of circulating
# genes) after an intervention is lifted. The comparison is done to the control
# experiment. 

design_irs <- create_intervention_scheme_IRS(PS_benchmark = '27', scenario_benchmark = scenario, IRS_START_TIMES = '29160', immigration_range=c(0), length_range=c(720,1800,3600), coverage_range=0.9, write_to_file = F, design_ref=design)
# Add control to the design_irs data frame, which is created in build_parameter_files.R
design_irs %<>% slice(rep(1, each = 1)) %>% bind_rows(design_irs)
design_irs[1, 'exp'] <- '001'
design_irs[1, 'IRS_length'] <- 0

intervention_stats_scenarios <- c()
intervention_stats_diff_scenarios <- c()
for (scenario in c('S','G','N')){
  for (run in run_range){
    for (ps in ps_range){
      if (file.exists(paste('~/Documents/malaria_interventions_data/PS',ps,'_',scenario,'_E001_R',run,'.sqlite',sep=''))){
        control_data <- get_data(parameter_space = ps, scenario = scenario, experiment = '001', run = run)[[1]]
      }
      for (e in exp_range){
        print(paste(Sys.time(),scenario, run, ps, e,sep=' | '))
        if (file.exists(paste('~/Documents/malaria_interventions_data/PS',ps,'_',scenario,'_E',e,'_R',run,'.sqlite',sep=''))){
          x <- post_intervention_stats(PS = ps, scenario = scenario, exp=e, run = run, plot.it = F, control_data = control_data, design_irs = design_irs)
          intervention_stats_scenarios <- rbind(intervention_stats_scenarios, x$summary_stats)
          intervention_stats_diff_scenarios <- rbind(intervention_stats_diff_scenarios, x$diff_control)
        } else {
          cat('missing file: ');cat(paste('PS',ps,'_',scenario,'_E',e,'_R',run,'.sqlite',sep=''));cat('\n')
        }
      }
    }
  }
}


## @knitr TIME_TO_EXTINCTION_SCENARIOS_PLOT

# Time to extinction as a function of the diversity
intervention_stats_scenarios %>%
  mutate(scenario=factor(scenario, levels=c('S','N','G'))) %>%
  left_join(subset(design, select=c(PS,BITING_RATE_MEAN,N_GENES_INITIAL)), by='PS') %>%
  group_by(PS, BITING_RATE_MEAN, N_GENES_INITIAL, scenario, exp) %>%
  summarise(time_ext_max=max(time_extinct), time_ext_mean=mean(time_extinct)) %>%
  gather(variable, value, -exp, -PS, -scenario, -run, -BITING_RATE_MEAN, -N_GENES_INITIAL) %>%
  ggplot(aes(x=N_GENES_INITIAL, y=value, color=scenario))+
  geom_point(size=4)+
  geom_line()+
  scale_x_continuous("Gene pool size", breaks = seq(1200,15600,1200)) +
  labs(y='Mean and max time to extinction')+
  scale_color_manual(values=scenario_cols)+
  facet_grid(variable~exp)+
  mytheme+
  theme(axis.text.x = element_text(angle = 90, hjust=0))

## @knitr EXTINCTION_PROB_SCENARIOS_PLOT

# Calculate probability of extinction as a proportion of runs which went extinct
intervention_stats_scenarios %>%
  mutate(scenario=factor(scenario, levels=c('S','N','G'))) %>%
  left_join(subset(design, select=c(PS,BITING_RATE_MEAN,N_GENES_INITIAL)), by='PS') %>%
  left_join(subset(design_irs, select=c(exp,IRS_START_TIMES,IRS_length))) %>%
  distinct(scenario, PS, N_GENES_INITIAL, exp, run, time_extinct,IRS_START_TIMES,IRS_length) %>%
  mutate(extinct=ifelse(time_extinct<as.numeric(IRS_START_TIMES)+as.numeric(IRS_length),1,0)) %>%
  group_by(scenario,PS, N_GENES_INITIAL, exp) %>%
  summarise(extinct_prob=sum(extinct)/length(extinct)) %>% 
  ggplot(aes(x=N_GENES_INITIAL, y=extinct_prob, color=scenario))+
  geom_point(size=4)+
  geom_line()+
  scale_x_continuous("Gene pool size", breaks = seq(1200,15600,1200)) +
  scale_color_manual(values=scenario_cols)+
  labs(y='Extinction probability')+
  facet_wrap(~exp)+
  mytheme+
  theme(axis.text.x = element_text(angle = 90, hjust=0))

## @knitr IRS_EFFECT_RATIO_SCENARIOS_PLOT

intervention_stats_diff_scenarios %>%
  mutate(scenario=factor(scenario, levels=c('S','N','G'))) %>%
  left_join(subset(design, select=c(PS,BITING_RATE_MEAN,N_GENES_INITIAL)), by='PS') %>%
  left_join(subset(design_irs, select=c(exp,IRS_START_TIMES,IRS_length))) %>%
  mutate(time_threshold=as.numeric(IRS_START_TIMES)+as.numeric(IRS_length)+360) %>%
  filter(time>=time_threshold) %>%
  filter(variable %in% c('prevalence', 'meanMOI','n_circulating_strains', 'n_circulating_genes', 'n_alleles')) %>%
  group_by(scenario, PS, BITING_RATE_MEAN, N_GENES_INITIAL, scenario, exp, variable) %>%
  summarise(mean_diff=mean(diff),
            mean_ratio=mean(ratio),
            sd_ratio=sd(ratio)) %>%
  mutate(error_up=mean_ratio+sd_ratio,
         error_low=mean_ratio-sd_ratio) %>%
  ggplot(aes(x=as.numeric(N_GENES_INITIAL), y=mean_ratio, color=scenario, group=scenario))+
  geom_point(size=3)+
  geom_line()+
  scale_x_continuous("Gene pool size", breaks = seq(1200,15600,1200)) +
  geom_errorbar(aes(ymin=error_low, ymax=error_up), width=0.1)+
  # geom_line(aes(x=as.numeric(N_GENES_INITIAL), y=error_up, color=scenario, group=scenario), linetype='dashed')+
  # geom_line(aes(x=as.numeric(N_GENES_INITIAL), y=error_low, color=scenario, group=scenario), linetype='dashed')+
  scale_color_manual(values=scenario_cols)+
  facet_grid(exp~variable)+
  mytheme+
  theme(axis.text.x = element_text(angle = 90, hjust=0))




#----------------------------------------------------------------------------
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




# FUNCTIONS for network Structure ---------------------------------------------------------------


overlapAlleleAdj<-function(mat){
  newmat<-tcrossprod(mat>0)
  newmat<-newmat/rowSums(mat>0)
  return(newmat)  
}


# A function to build the similarity matrix for a single layer and calculate some summary stats
build_layer <- function(l){
  sampled_infections_layer <- subset(sampled_infections, layer==l) # This is MUCH faster than sampled_infections_layer <- sampled_infections %>% filter(layer==l)
  sampled_infections_layer %<>% group_by(strain_id) %>%
    mutate(freq = n()/120) %>% # strain frequency (the number of strain copies should be equal to the frequency)
    arrange(strain_id_unique) 
  # Calculate the edges
  similarity_matrix <- table(sampled_infections_layer$strain_id_unique, sampled_infections_layer$allele_locus)
  similarity_matrix <- overlapAlleleAdj(similarity_matrix)
  # Some summary
  layer_summary <- with(sampled_infections_layer, 
                        data.frame(hosts=length(unique(host_id)),
                                   repertoires=length(unique(strain_id)),
                                   total_infections=length(unique(strain_id_unique))
                        ))
  return(list(similarity_matrix=similarity_matrix, infections=sampled_infections_layer, layer_summary=layer_summary))
}


createTemporalNetwork <- function(ps, scenario, exp, run, cutoff_prob=0.9, write_files=F, layers_to_include=NULL){
  base_name <- paste('PS',ps,'_',scenario,'_E',exp,'_R',run,sep='')
  sqlite_file <- paste(base_name,'.sqlite',sep='')
  
  print('Getting data from sqlite...')
  # Extract data from sqlite. variable names correspond to table names
  db <- dbConnect(SQLite(), dbname = sqlite_file)
  sampled_strains <- as.tibble(dbGetQuery(db, 'SELECT id, gene_id FROM sampled_strains'))
  names(sampled_strains)[1] <- c('strain_id')
  sampled_alleles <- as.tibble(dbGetQuery(db, 'SELECT * FROM sampled_alleles'))
  names(sampled_alleles)[3] <- c('allele_id')
  sampled_strains <- full_join(sampled_strains, sampled_alleles)
  sampled_strains$allele_locus <- paste(sampled_strains$allele_id,sampled_strains$locus,sep='_') # each allele in a locus is unique
  dbDisconnect(db)
  # Get infection data
  sampled_infections <- get_data(parameter_space = ps, experiment = exp, scenario = scenario, run = run)$sampled_infections
  
  print('Building data set...')
  # Add layer numbers. Each time point is a layer
  sampled_infections$layer <- group_indices(sampled_infections, time) 
  sampled_infections %<>% arrange(layer,strain_id,host_id) %>% 
    group_by(layer,strain_id) %>% 
    mutate(strain_copy = row_number()) # add unique id for each strain copy within each layer. A strain copy is an instance of a strain in a particualr host
  sampled_infections$strain_id_unique <- paste(sampled_infections$strain_id,sampled_infections$strain_copy,sep='_') # Create the unique strains
  # Integrate the strain composition into the infections table
  if (all(unique(sampled_strains$strain_id)%in%unique(sampled_infections$strain_id))==F || all(unique(sampled_infections$strain_id)%in%unique(sampled_strains$strain_id))==F) {
    stop('There may be a mismatch in repertoires between the sampled_strains and sampled_infections data sets. Revise')
  }
  sampled_infections %<>% select(time, layer, host_id, strain_id, strain_copy, strain_id_unique) %>% 
    left_join(sampled_strains)
  
  print('Building layers...')
  Layers <- list()
  if (is.null(layers_to_include)){
    layers_to_include <- sort(unique(sampled_infections$layer))
  }
  for (l in layers_to_include){
    cat(paste('[',Sys.time(), '] building layer ',l,'\n',sep=''))
    Layers[[l]] <- build_layer(l)  
  }
  # Get just the matrices
  temporal_network <- unique(lapply(Layers, function(x) head(x, 1)$similarity_matrix))
  sapply(temporal_network, nrow)
  
  print('Applying cutoff...')
  # Apply a cutoff
  ## Get the similarity matrix for all the repertoires and alleles
  x <- xtabs(~strain_id+allele_locus, subset(sampled_strains, strain_id%in%sampled_infections$strain_id))
  similarityMatrix <- overlapAlleleAdj(x)
  dim(similarityMatrix)
  
  cutoffValue <- quantile(as.vector(similarityMatrix), probs = 0.9)
  print(cutoffValue)
  # hist(as.vector(similarityMatrix));abline(v=cutoffValue,col='red')
  similarityMatrix[similarityMatrix<cutoffValue] <- 0
  for (i in 1:length(temporal_network)){
    print(i)
    x <- temporal_network[[i]]
    x[x<cutoffValue] <- 0
    temporal_network[[i]] <- x
  }
  
  if(write_files){
    print('Writing similarit matrices to files...')
    # Write similarity matrix without cutoff. This is used to plot the edge weight distributions.
    write.table(similarityMatrix, paste('../',filenameBase,'_similarityMatrix_nocutoff.csv',sep=''), row.names = T, col.names = T, sep=',')
    # write.table(similarityMatrix, paste(filenameBase,'_similarityMatrix.csv',sep=''), row.names = T, col.names = T, sep=',')
    # write.table(similarityMatrix, paste('../',filenameBase,'_similarityMatrix.csv',sep=''), row.names = T, col.names = T, sep=',')
  }
  
  print('Done!')
  return(list(temporal_network=temporal_network, base_name = base_name, ps=ps, scenario=scenario, experiment=exp, run=run))
}


# This function takes a list of temporal matrices and returns an intralayer and
# interlayer edge lists in a format: [layer_source, node_source, layer_target node_target, weight].
# It also returns the list of node names
build_infomap_objects <- function(network_object, write_to_infomap_file=T, return_objects=T){
  temporal_network <- network_object$temporal_network
  base_name <- network_object$base_name
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
  # network_object is a list of matrices, each element in the list is a layer.
  matrix_to_infomap <- function(l, nodeList, network_object){
    require(igraph)
    current_layer <- network_object[[l]]
    if(nrow(current_layer)<2){
      print(paste('Less than 2 repertoires in layer',l,'!!! skipping intralayer edges'))
      next
    }
    g <- graph.adjacency(current_layer, mode = 'directed', weighted = T)
    current_layer_el <- as.tibble(as_data_frame(g, what = 'edges'))
    names(current_layer_el) <- c('node_s','node_t','w')
    current_layer_el$layer_s <- l
    current_layer_el$layer_t <- l
    current_layer_el$node_s <- nodeList$nodeID[match(current_layer_el$node_s,nodeList$nodeLabel)]
    current_layer_el$node_t <- nodeList$nodeID[match(current_layer_el$node_t,nodeList$nodeLabel)]
    current_layer_el %<>% select(layer_s, node_s, layer_t, node_t, w) # Re-arrange columns for the infomap input order
    print(paste('[',Sys.time(), '] Created edge list of layer ',l,' for Infomap | ', nrow(current_layer_el),' edges',sep=''))
    return(current_layer_el)
  }
  
  ## A function that calculates interlayer edges between layer t and t+1
  build_interlayer_edges_1step <- function(t, nodeList, network_object){
    strain_copies_t <- rownames(network_object[[t]]) # repertoires at time t
    strain_copies_t1 <- rownames(network_object[[t+1]]) # repertoires at time t+1
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
  infomap_intralayer <- lapply(layers, function (x) matrix_to_infomap(x, nodeList = nodeList, network_object = temporal_network))
  infomap_intralayer <- do.call("rbind", infomap_intralayer)
  
  # Create interlayer edge lists. Only connect consecutive layers (interlayer ONLY go from layer l to l+1).
  layers <- layers[-length(layers)]
  infomap_interlayer <- lapply(layers, function (x) build_interlayer_edges_1step(x, nodeList = nodeList, network_object = temporal_network))
  infomap_interlayer <- do.call("rbind", infomap_interlayer)
  
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


# Example for structure ---------------------------------------------------

ps <- '36'
scenario <- 'S'
exp <- '001'
run <- 1

# parameter_file <- paste(base_name,'.py',sep='') # This may be necessary so I leave it


network_test <- createTemporalNetwork(ps,scenario,exp,1, layers_to_include = c(1:20,30:40))
sapply(network_test, nrow)
infomap_objects <- build_infomap_objects(network_test[1:20])


# Run Infomap

# Get infomap results and process them






















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
