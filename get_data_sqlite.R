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
  axis.text = element_text(color='black', family="Helvetica", size=19),
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
setwd('~/GitHub/')
PS <- '01'
scenario <- 'S'
exp <- '01' # 00 is for the checkpoint and control
run <- 1
base_name <- paste('PS',PS,'_',scenario,'_E',exp,sep='')
sqlite_file <- paste(base_name,'_R',run,'.sqlite',sep='')
parameter_file <- paste(base_name,'.py',sep='')

# Extract from sqlite -----------------------------------------------------
db <- dbConnect(SQLite(), dbname = sqlite_file)
sampled_hosts <- dbGetQuery(db, 'SELECT * FROM sampled_hosts')
summary_general <- dbGetQuery(db, 'SELECT * FROM summary')
summary_general <- inner_join(sampled_hosts, summary_general)

# Main code ---------------------------------------------------------------
summary_general <- summary_general[-1,]
# Prevalence
summary_general$prevalence <- summary_general$total_infected/10^4

#EIR
#EIR=biting_rate * prevalence
biting_rate <- get_biting_rate(parameter_file)
summary_general$b <- NA
summary_general$b[] <- biting_rate # The [] is for recycling the biting_rate
summary_general$EIR <- summary_general$prevalence*summary_general$b*30

# Burnin
# burnin_seq <- 1:10  # This is burnin in years, not days
# df <- map(burnin_seq, function(burnin){
#   x <- summary_general %>% filter(time>burnin*360)
#   x$burnin <- burnin
#   return(x)
# }) %>% bind_rows()
# df$burnin <- as.factor(df$burnin)
# df %>% ggplot(aes(time, prevalence, color=burnin))+geom_line()+facet_wrap(~burnin)
# df %>% ggplot(aes(time, EIR, color=burnin))+geom_line()+facet_wrap(~burnin)

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


summary_general$PS <- PS
summary_general$exp <- exp
summary_general$scenario <- scenario
summary_general$run <- run
assign(paste(base_name,'_R',run,'_results',sep=''), summary_general)

# Compare between experiments ---------------------------------------------

d <- rbind(PS01_S_E00_R1_results,PS01_S_E01_R1_results)
# mintime=d %>% group_by(exp) %>% summarise(m=max(time)) %>% summarise(min(m))
# mintime=mintime[1,1]
pdf('seasonal_comparison.pdf',16,10)
d %>%
  select(-year, -month, -n_infected) %>% 
  # filter(time>15000) %>% 
  gather(variable, value, -time, -exp, -PS, -scenario, -run) %>% 
  ggplot(aes(x=time, y=value, color=exp, group=exp))+
  geom_line()+
  scale_x_continuous(breaks=pretty(x=d$time,n=20))+
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
