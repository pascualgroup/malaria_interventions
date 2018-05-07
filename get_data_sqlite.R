library(sqldf)
library(tidyverse)


# Functions ---------------------------------------------------------------

chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 

# This function extracts the biting rates from the parameter file, which are in
# daily resolution and returns a vector of the averaged biting rates for a given
# period. (e.g. 30 would be biting rates averaged over a month and 1 will return
# the daily biting rates).
get_biting_rate <- function(parameter_file, sampling_period=30){
  x <- readLines(parameter_file)
  y <- x[grep('BITING_RATE_MEAN',x)[1]]
  BITING_RATE_MEAN <- as.numeric(str_extract(y, "[[:digit:]]"))
  y <- x[grep('DAILY_BITING_RATE_DISTRIBUTION',x)[1]]
  DAILY_BITING_RATE_DISTRIBUTION <- eval(parse(text=paste('c(',(str_sub(y, str_locate(y, '\\[(.*?)\\]')[1]+1, str_locate(y, '\\[(.*?)\\]')[2]-1)),')',sep='')))
  BITING_RATE <- BITING_RATE_MEAN*DAILY_BITING_RATE_DISTRIBUTION
  BITING_RATE <- chunk2(BITING_RATE, 360/sampling_period)
  sapply(BITING_RATE, mean)
}


# Arguments ---------------------------------------------------------------
setwd('~/Documents/malaria_interventions_sqlite')
exp <- 'test_05'
run <- 1
sqlite_file <- paste(exp,'_',run,'.sqlite',sep='')
parameter_file <- paste(exp,'.py',sep='')


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
summary_general$b <- biting_rate
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

summary_general %>% filter(time>1440) %>% gather(variable, value, -time) %>% ggplot(aes(time, value))+geom_line()+facet_wrap(~variable, scales = 'free')


# Annual biting rate is given by taking an average over 12 months and multiplying by 30.
summary_general$year <- gl(n = max(summary_general$time)/360, length = nrow(summary_general), k = 1)
summary_general %>% group_by(year) %>% summarise(eir_y=mean(EIR)*30)


summary_general_04$exp <- 'test4'
summary_general_05$exp <- 'test5'
d <- rbind(summary_general_04,summary_general_05)
d %>% filter(time>1440) %>% gather(variable, value, -time, -exp) %>% ggplot(aes(time, value, color=exp))+geom_line()+facet_wrap(~variable, scales = 'free')
