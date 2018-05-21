library(sqldf)
library(tidyverse)


# Functions ---------------------------------------------------------------
mytheme <- theme_bw() + theme(
  legend.title  = element_text(colour = "black", size=17),
   # legend.position = "none",
  #	legend.direction = "horizontal",
  legend.key = element_blank(),
  legend.text  = element_text(colour = "black", size=17),
  panel.background = element_blank(),
  panel.grid.minor = element_blank(),
  axis.text = element_text(color='black', family="Helvetica", size=10),
  strip.text.x = element_text(family = "Helvetica", size = 10),
  strip.text.y = element_text(family = "Helvetica", size = 10),
  panel.border = element_rect(colour = "black", size=1.3),
  axis.ticks = element_line(size = 1.3),
  strip.background = element_rect( fill = "transparent", size = 1.3, colour = "black"  ),
  strip.text = element_text(size = 19)
)

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


# Arguments ---------------------------------------------------------------
setwd('~/Documents/malaria_interventions_sqlite')
# setwd('~/GitHub/')
parameter_space <- '03'
scenario <- 'N'
experiment <- '02' # 00 is for the checkpoint and control
run <- 1
base_name <- paste('PS',parameter_space,'_',scenario,'_E',experiment,'_R',run,sep='')
sqlite_file <- paste(base_name,'.sqlite',sep='')
parameter_file <- paste(base_name,'.py',sep='')

# Extract from sqlite -----------------------------------------------------
db <- dbConnect(SQLite(), dbname = sqlite_file)
# sampled_hosts <- dbGetQuery(db, 'SELECT * FROM sampled_hosts')
summary_general <- dbGetQuery(db, 'SELECT * FROM summary')
sampled_infections <- dbGetQuery(db, 'SELECT * FROM sampled_infections')
# summary_general <- inner_join(sampled_hosts, summary_general)

# Main code ---------------------------------------------------------------
if (nrow(summary_general)%%30==1){# The beginning of the data set has a timestap of 0 and this line is unnecesary. So remove it
  summary_general <- summary_general[-1,]
}

# Prevalence
summary_general$prevalence <- summary_general$n_infected/10^4

#EIR
#EIR=biting_rate * prevalence
biting_rate <- get_biting_rate(parameter_file)
summary_general$b <- NA
summary_general$b[] <- biting_rate # The [] is for recycling the biting_rate
summary_general$EIR <- summary_general$prevalence*summary_general$b*30

# MOI
meanMOI <- sampled_infections %>% group_by(time, host_id) %>% summarise(MOI=length(strain_id)) %>% group_by(time) %>% summarise(meanMOI=mean(MOI))
summary_general <- inner_join(summary_general, meanMOI)

# Host age structure
hosts <- dbGetQuery(db, 'SELECT * FROM hosts')
names(hosts)[1] <- 'host_id'
hosts$lifespan <- round((hosts$death_time-hosts$birth_time)/30)
hosts <- subset(hosts, host_id%in%sampled_infections$host_id)
sampled_infections <- left_join(sampled_infections, hosts, by='host_id')
sampled_infections$host_age <- round((sampled_infections$time-sampled_infections$birth_time)/30)

sampled_infections %>% ggplot(aes(x=host_age))+geom_histogram() + labs(x='Infected host age (months)') # Plot the age structure of infected hosts

summary_general %>% 
  select(-n_infected) %>% 
  # filter(time>10000) %>%
  gather(variable, value, -time) %>% 
  ggplot(aes(time, value, color=variable))+
  geom_line()+
  facet_wrap(~variable, scales = 'free')+
  mytheme+theme(legend.position = 'none')

# Annual biting rate is given by taking an average over 12 months and multiplying by 30.
summary_general$year <- gl(n = max(summary_general$time)/360, length = nrow(summary_general), k = 1)
summary_general %>% group_by(year) %>% summarise(eir_y=mean(EIR)*30)

summary_general$month <- gl(n = 12, k = 1, length = nrow(summary_general),labels = c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'), ordered = F)
summary_general %>% 
  ggplot(aes(x=month,y=EIR))+
  geom_boxplot()+
  geom_point(stat='summary', fun.y=mean, color='red')+
  stat_summary(fun.y=mean, geom="line")+mytheme

# Create an object with the results
summary_general$PS <- parameter_space
summary_general$exp <- experiment
summary_general$scenario <- scenario
summary_general$run <- run
assign(paste('results_',base_name,sep=''), summary_general)

# Compare between experiments ---------------------------------------------

d <- rbind(results_PS03_S_E01_R1,
           # results_PS03_S_E02_R1,
           results_PS03_S_E03_R1,
           # results_PS03_S_E04_R1)
           results_PS03_S_E05_R1)

d <- rbind(results_PS03_N_E00_R1,results_PS03_N_E01_R1,results_PS03_N_E02_R1)
# mintime=d %>% group_by(exp) %>% summarise(m=max(time)) %>% summarise(min(m))
# mintime=mintime[1,1]
# pdf('seasonal_comparison.pdf',16,10)
time_range <- c(10800,14400)
d %>%
  select(-year, -month, -n_infected) %>% 
  filter(time>time_range[1]&time<time_range[2]) %>%
  gather(variable, value, -time, -exp, -PS, -scenario, -run) %>% 
  ggplot(aes(x=time, y=value, color=exp, group=exp))+
  geom_line()+
  # geom_vline(xintercept = c(21600,21960,22320,22680,23040,23400))+
  scale_x_continuous(breaks=pretty(x=subset(d, time>time_range[1]&time<time_range[2])$time,n=5))+
  facet_wrap(~variable, scales = 'free')+mytheme
d %>% 
  ggplot(aes(x=month,y=EIR, color=exp, group=exp))+
  # geom_boxplot()+
  geom_point(stat='summary', fun.y=mean)+
  scale_y_continuous(limits = c(0,10))+
  stat_summary(fun.y=mean, geom="line")+
  mytheme
dev.off()


png('~/Documents/malaria_interventions/burnin_test.png', 1200, 800)
d %>% filter(time>1440) %>% gather(variable, value, -time, -exp) %>% ggplot(aes(time, value, color=exp))+geom_line()+facet_wrap(~variable, scales = 'free')
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
calendar %>% filter(running_day==10930)
