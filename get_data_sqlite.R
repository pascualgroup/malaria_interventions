library(tidyverse)
library(magrittr)
library(sqldf)
library(igraph)

# Functions ---------------------------------------------------------------

## @knitr FUNCTIONS
source('functions.R')
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

# Add control to the design_irs data frame
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





# Example for structure ---------------------------------------------------



setwd('~/Documents/malaria_interventions_data/')

ps <- '36'
scenario <- 'S'
exp <- '004'
run <- 1
x <- get_data(parameter_space = ps, scenario, experiment = exp, run)[[1]]
x$layer <- group_indices(x, time) 
layers_to_include = 1:max(x$layer)
x %>% 
  filter(time>28000 & time<35000) %>% 
  ggplot(aes(x=layer, y=n_circulating_strains))+
  geom_line(color='#900C3F', size=1)+
  labs(x='Time',y='Diversity')+
  scale_x_continuous(breaks = seq(0,max(x$layer),5))+
  mytheme

network_control <- createTemporalNetwork(ps, scenario, exp = '001', run, layers_to_include = layers_to_include, sampled_infections = NULL)
x <- as.tibble(network_control$layer_summary)
x %>% gather(variable, value, -layer) %>% 
  ggplot(aes(layer, value, color=variable))+
  geom_line()+
  scale_x_continuous(breaks = seq(1,max(x$layer)+5,5))+
  mytheme+
  theme(panel.grid.minor = element_blank())

network_test <- createTemporalNetwork(ps, scenario, exp, run, cutoff_value = network_control$cutoff_value, layers_to_include = layers_to_include, sampled_infections = NULL)
x <- as.tibble(network_test$layer_summary)
x %>% gather(variable, value, -layer) %>% 
  ggplot(aes(layer, value, color=variable))+
  geom_line()+
  scale_x_continuous(breaks = seq(1,max(x$layer)+5,5))+
  mytheme+
  theme(panel.grid.minor = element_blank())

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

network_test$base_name
network_control$base_name

netork_properties_control <- matrix(0,ncol=23, nrow=length(layers_to_include))
for (i in layers_to_include){
  print(i)
  tmp <- calculateFeatures(network_control, i)
  netork_properties_control[i,] <- tmp
}
colnames(netork_properties_control) <- names(tmp)
netork_properties_control <- as.tibble(netork_properties_control)
netork_properties_control$layer <- layers_to_include
netork_properties_control  %>%  gather(variable, value, -layer) %>% 
  ggplot(aes(layer, value, color=variable))+
  geom_line()+
  scale_x_continuous(breaks = seq(1,max(layers_to_include),5))+
  facet_wrap(~variable,scales='free')+
  geom_vline(xintercept = c(13,73))+
  mytheme+
  theme(panel.grid.minor = element_blank(), legend.position = 'none')

netork_properties <- matrix(0,ncol=23, nrow=length(layers_to_include))
for (i in layers_to_include){
  print(i)
  tmp <- calculateFeatures(network_test, i)
  netork_properties[i,] <- tmp
}
colnames(netork_properties) <- names(tmp)
netork_properties <- as.tibble(netork_properties)
netork_properties$layer <- layers_to_include
netork_properties  %>%  gather(variable, value, -layer) %>% 
  ggplot(aes(layer, value, color=variable))+
  geom_line()+
  scale_x_continuous(breaks = seq(1,max(layers_to_include),5))+
  facet_wrap(~variable,scales='free')+
  geom_vline(xintercept = c(13,73))+
  mytheme+
  theme(panel.grid.minor = element_blank(), legend.position = 'none')

netork_properties_control$exp <- '001'
netork_properties$exp <- '003'
comparison <- rbind(netork_properties,netork_properties_control)
comparison %>% 
  select('exp','layer','Num_nodes','GCC','density','mean_degree','diameter','M07--A<-B->C','M06--A<->B<-C','M10--A<->B->C','M11--A<->B<->C','M16--A<->B<->C_A<->C') %>%
  gather(variable, value, -exp, -layer) %>% 
  ggplot(aes(layer, value, color=exp))+
  geom_line()+
  scale_x_continuous(breaks = seq(1,max(layers_to_include),20))+
  facet_wrap(~variable,scales='free')+
  geom_vline(xintercept = c(13,13+120))+
  mytheme+
  theme(panel.grid.minor = element_blank(), legend.position = 'none')


unlist(lapply(network_test$temporal_network, nrow))



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
