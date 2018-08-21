library(tidyverse)
library(magrittr)
library(sqldf)
setwd('~/Documents/malaria_interventions')

source('functions.R')

# Create parameter and job files -------------------------------------------

setwd('~/Documents/malaria_interventions_data/')

# Clear previous files if necessary
clear_previous_files(run = 6, exclude_sqlite = F, exclude_CP = F, exclude_control = F, test = T)
# Get data design 
design <- loadExperiments_GoogleSheets() 

# Create the reference experiments (checkpoint and control)
ps_range <- sprintf('%0.2d', 27:39)
exp_range <- sprintf('%0.3d', 0:4)
run_range <- 11:50
work_scenario <- 'N'
design_subset <- subset(design, PS %in% ps_range & scenario==work_scenario)
generate_files(row_range = 1:nrow(design_subset), run_range = run_range, experimental_design = design_subset)

# # If checkpoints already exist, can create the reference experiments (control only)
for (ps in ps_range){
  print(ps)
  design_control <- subset(design, PS==ps & scenario == work_scenario & exp=='001')
  # clear_previous_files(parameter_space = ps, scenario = work_scenario, exclude_sqlite = F, exclude_CP = T, exclude_control = F, test = F)
  seeds <- get_random_seed(PS = ps, scenario = work_scenario, run_range = run_range)
  generate_files(row_range = 1, run_range = run_range, experimental_design = design_control, random_seed = seeds)
}

# Create the corresponding IRS experiments  
for (ps in ps_range){
  design_irs <- create_intervention_scheme_IRS(PS_benchmark = ps, scenario_benchmark = work_scenario, IRS_START_TIMES = '29160', immigration_range=c(0), length_range=c(720,1800,3600), coverage_range=0.9, poolsize_bounce = 'False', write_to_file = F, design_ref=design)
  generate_files(row_range = 1:nrow(design_irs), run_range = run_range, 
                 random_seed = get_random_seed(ps, work_scenario, run_range = run_range), 
                 experimental_design=design_irs)
}

# ZIP all the PY and sbatch files
unlink('files_to_run.txt')
sink('files_to_run.txt', append = T)
# Add sbatch files
files <- list.files(path = '~/Documents/malaria_interventions_data/', pattern = 'sbatch', full.names = F) 
files <- files[str_detect(string = files, paste(work_scenario,'E',sep=''))]
files <- files[str_sub(files,3,4) %in% ps_range]
files <- files[str_sub(files,7,9) %in% exp_range]
for (i in 1:length(files)){
  cat(files[i]);cat('\n')
}
# Add py files
files <- list.files(path = '~/Documents/malaria_interventions_data/', pattern = '.py', full.names = F) 
files <- files[str_detect(string = files, paste('_',work_scenario,sep=''))]
files <- files[str_sub(files,3,4) %in% ps_range]
files <- files[str_sub(files,9,11) %in% exp_range]
files <- files[parse_number(str_sub(files,14,15)) %in% run_range]
for (i in 1:length(files)){
  cat(files[i]);cat('\n')
}
sink.reset()
# Create zip
unlink('files_to_run.zip')
system('zip files_to_run.zip -@ < files_to_run.txt')

# Copy the file to Midway and unzip it

# First run the checkpoints
paste("for i in ",paste(ps_range,collapse=' '),"; do sbatch 'PS'$i'",work_scenario,"E000.sbatch'; done;",sep='')

# Then run control and interventions
## Run in Midway terminal:
cat("rm job_ids.txt; sacct -u pilosofs --format=jobid,jobname --starttime 2018-08-07T14:48:00 --name=");cat(paste(paste(ps_range,work_scenario,'E000',sep=''),collapse = ','));cat(" >> 'job_ids.txt'")
## Copy file from Midway and run:
jobids <- read.table('job_ids.txt', header = F, skip=2) 
jobids <- na.omit(unique(parse_number(jobids$V1))) # 1 job id per PS
length(jobids)==length(ps_range)
## Copy the output of the following loop and paste in Midway
for (ps in ps_range){
  for (e in exp_range){
    cat(paste('sbatch -d afterok:',jobids[which(ps_range==ps)],' PS',ps,work_scenario,'E',e,'.sbatch',sep=''));cat('\n')
  }
}

# Or, if checkpoints are already finished:
unlink('jobs_to_run.sh')
sink('jobs_to_run.sh')
for (ps in ps_range){
  for (e in exp_range[-1]){
    cat(paste('sbatch PS',ps,work_scenario,'E',e,'.sbatch',sep=''));cat('\n')
  }
}
sink.reset()

# Zip files to reduce the clutter -----------------------------------------

system('rm *.sbatch')

scenario_range <- c('G')
exp_range <- sprintf('%0.3d', 0:4)
ps_range <- sprintf('%0.2d', 27:39)

for (ps in ps_range){
  for (scenario in scenario_range){
    unlink('files_tmp.txt')
    sink('files_tmp.txt', append = T)
    # Add py files
    files <- list.files(path = '~/Documents/malaria_interventions_data/', pattern = '.py', full.names = F) 
    files <- files[str_sub(files,3,4) %in% ps]
    files <- files[str_sub(files,6,6) %in% scenario]
    files <- files[str_sub(files,9,11) %in% exp_range]
    for (i in 1:length(files)){
      cat(files[i]);cat('\n')
    }
    sink.reset()
    # Create zip
    # unlink(paste(ps,'_sbatch_py.zip',sep=''))
    if (file.size('files_tmp.txt')>10){
      system(paste('zip ',ps,'_',scenario,'_py.zip -@ < files_tmp.txt -mu',sep=''))
    }
  }
}


# Verify files ------------------------------------------------------------

# --- sqlite files ---
files <- list.files(path = '~/Documents/malaria_interventions_data/', pattern = '\\.sqlite', full.names = F)
# Also consider adding following information: extracting random seeds; existence of corresponding parameter files
files_sqlite <- tibble(file_sqlite=files,
                    PS = sapply(str_split(files,'_'),function (x) parse_number(x[1])),
                    scenario=sapply(str_split(files,'_'),function (x) x[2]),
                    experiment=sapply(str_split(files,'_'),function (x) str_sub(x[3],2,4)),
                    run= sapply(str_split(files,'_'),function (x) parse_number(x[4])),
                    CP=sapply(str_split(files,'_'),function (x) str_detect(x[5],'CP')),
                    size=round(file.size(files)/1024^2,2)
)
files_sqlite$PS <- sprintf('%0.2d', files_sqlite$PS)

subset(files_sqlite, size<1)

files_sqlite$run_time <- unlist(lapply(files_sqlite$file_sqlite[is.na(files_sqlite$CP)], function (f) {
  print(f)
  db <- dbConnect(SQLite(), dbname = paste('~/Documents/malaria_interventions_data/',f,sep=''))
  x <- dbGetQuery(db, 'SELECT max(time) FROM summary')
  dbDisconnect(db)
  x
}))
unique(files_sqlite$run_time)
files_sqlite %>% filter(is.na(run_time))


files_sqlite %<>% mutate(scenario=factor(scenario, levels=c('S','N','G')))  %>% arrange(CP, PS, scenario, experiment, run)         
files_sqlite %>% filter(is.na(CP) & as.numeric(PS)>=27) %>% 
  filter(scenario=='S') %>% 
  group_by(PS,scenario,experiment) %>% 
  summarise(runs_completed=length(run)) %>% 
  filter(runs_completed<50) %>% print(n = Inf)
files_sqlite %>% group_by(scenario, PS, experiment) %>% summarise(x=sum(runs_completed)) %>% print(n = Inf)
files_sqlite %>% filter(scenario=='G') %>% print(n = Inf)

#--- PY files ---
files <- list.files(path = '~/Documents/malaria_interventions_data/', pattern = 'py.zip', full.names = T)
files <- map(files, function(f){
  unzip(f, list=T)
}) %>% bind_rows()

files_py <- tibble(file_py=files$Name,
                   PS = sapply(str_split(files$Name,'_'),function (x) parse_number(x[1])),
                   scenario=sapply(str_split(files$Name,'_'),function (x) x[2]),
                   experiment=sapply(str_split(files$Name,'_'),function (x) str_sub(x[3],2,4)),
                   run= sapply(str_split(files$Name,'_'),function (x) parse_number(x[4])),
                   size=round(files$Length,2),
                   date=files$Date
)
files_py$PS <- sprintf('%0.2d', files_py$PS)

files_py %<>% filter(as.numeric(PS)>=27) %>% mutate(scenario=factor(scenario, levels=c('S','N','G'))) %>% arrange(scenario, PS, experiment, run)
files_py %<>% filter(as.numeric(PS)>=27) %>% group_by(PS,scenario,experiment) %>% summarise(runs_count=length(run))
files_py %>% filter(scenario=='G') %>% filter(as.numeric(PS)>=27) %>% distinct(PS,run)

# --- Generate missing GI files --- For those runs that produced an error when
# fitting the curve, take curve-fit information from a random run where fit was
# successful.

existing_py_G <- files_py %>% filter(scenario=='G') %>% filter(as.numeric(PS)>=27 & experiment=='000') %>% distinct(PS,run)
py_files_G <- vector(mode = 'list', length = length(ps_range))
names(py_files_G) <- ps_range
for (ps in ps_range){
  existing <- subset(existing_py_G, PS==ps)$run
  py_files_G[[which(ps_range==ps)]]$existing <- existing
  py_files_G[[which(ps_range==ps)]]$missing <- setdiff(1:50,existing)
}


experiment <- '000'
for (parameter_space in ps_range){
  py_files_G_ps <- py_files_G[[which(ps_range==parameter_space)]]
  for (run in py_files_G_ps$missing){
    benchmark_file <- paste('PS',parameter_space,'_G_E',exp,'_R',sample(py_files_G_ps$existing,1),'.py',sep='')
    x <- readLines(benchmark_file)
  }
}



# --- CP sqlite files ---
# In Midway run: ls *.sqlite >> CP_files.txt and copy the file to local computer
files <- read.table('CP_files.txt',header = F, stringsAsFactors = F)$V1
files_CP <- tibble(file_CP=files,
          PS = sapply(str_split(files,'_'),function (x) parse_number(x[1])),
          scenario=sapply(str_split(files,'_'),function (x) x[2]),
          experiment=sapply(str_split(files,'_'),function (x) str_sub(x[3],2,4)),
          run= sapply(str_split(files,'_'),function (x) parse_number(x[4])),
          CP=sapply(str_split(files,'_'),function (x) str_detect(x[5],'CP'))
)
files_CP$PS <- sprintf('%0.2d', files_CP$PS)

missing_py_files <- full_join(files_CP, files_py, by = c("PS", "scenario", "experiment", "run")) %>% 
  filter(as.numeric(PS)>=27) %>% 
  select(scenario, PS, experiment, run, file_CP, file_py) %>% 
  mutate(scenario=factor(scenario, levels=c('S','N','G'))) %>% 
  arrange(PS, scenario, experiment, run)

missing_py_files %>% filter(is.na(file_py))

#--- Check seeds ---
for(i in 1:nrow(files_df)){
  print(i)
  unzip(paste(files_df$PS[i],'_',files_df$scenario[i],'_py.zip',sep=''), files = files_df$file_py[i])
  files_df$seed[i] <- get_random_seed(PS = files_df$PS[i], scenario = as.character(files_df$scenario[i]), experiment = files_df$experiment[i], run_range = files_df$run[i])
  unlink(files_df$file_py[i])
}

seed_check <- files_df %>% distinct(scenario, PS, run, seed) %>% group_by(scenario, PS) %>% summarise(length(run))


# Generate sbatch files to create the missing checkpoints -----------------

ps_scen_combinations <- as.tibble(files_CP_missing) %>% distinct(scenario,PS)

for (i in 1:nrow(ps_scen_combinations)){
  ps <- ps_scen_combinations$PS[i]
  scen <- ps_scen_combinations$scenario[i]

  unzip(paste(ps,'_',scen,'_py.zip',sep=''), files = subset(files_CP_missing, scenario==scen & PS==ps)$file_py)
  experiment <- '000' # Use 000 for checkpoints and 001 for control
  base_name <- paste('PS',ps,scen,'E',experiment,sep='')

  SLURM_ARRAY_RANGE <- collapse_with_commas(subset(files_CP_missing, scenario==scen & PS==ps)$run, use_brackets = F)

  experimental_design <- subset(design, PS==ps & scenario==scen & exp=='000')
  # Write the job file for the exepriment
  job_lines <- readLines('~/Documents/malaria_interventions/job_file_ref.sbatch')
  wall_time <-  experimental_design$wall_time
  mem_per_cpu <-  experimental_design$mem_per_cpu
  job_lines[2] <- paste('#SBATCH --job-name=',paste(ps,scen,'E',experiment,sep=''),sep='')
  job_lines[3] <- paste('#SBATCH --time=',wall_time,sep='')
  job_lines[4] <- paste('#SBATCH --output=slurm_output/',base_name,'_%A_%a.out',sep='')
  job_lines[5] <- paste('#SBATCH --error=slurm_output/',base_name,'_%A_%a.err',sep='')
  job_lines[6] <- paste('#SBATCH --array=',SLURM_ARRAY_RANGE,sep='')
  job_lines[9] <- paste('#SBATCH --mem-per-cpu=',mem_per_cpu,sep='')
  job_lines[19] <- paste("PS='",ps,"'",sep='')
  job_lines[20] <- paste("scenario='",scen,"'",sep='')
  job_lines[21] <- paste("exp='",experiment,"'",sep='')
  job_lines[26] <- "CHECKPOINT='create'"
  output_file=paste(base_name,'.sbatch',sep = '')
  write_lines(job_lines, output_file)
  cat(paste('sbatch ',base_name,'.sbatch',sep=''));cat('\n')
}


# Add a line in the parameter files ---------------------------------------
# This adds a line to the end of each py file
ps_range <- sprintf('%0.2d', 28:39)
work_scenario <- 'G'
exp_range <- sprintf('%0.3d', 2:4)
run_range <- 1:10
for (ps in ps_range){
  for (e in exp_range){
    for (r in run_range){
      file <- paste('PS',ps,'_',work_scenario,'_E',e,'_R',r,'.py',sep='')
      print(file)
      unzip(paste(ps,'_',work_scenario,'_py.zip',sep=''), file)
      if(file.exists(file)) {write('POOLSIZE_BOUNCE_BACK_AFTER_INTERVENTION=False',file, append=TRUE)}
    }
  }
}

for (ps in ps_range){
  generate_sbatch(ps, scen = 'N', experiment = '003', runs = 1:10, unzip_py_files = F)
}

# Or generate sbatch for particular runs, as in the generalized immunity
cases <- files_py %>% filter(scenario=='G') %>% filter(as.numeric(PS)>27) %>% distinct(PS,run)
for (ps in ps_range){
  generate_sbatch(ps, scen = 'G', experiment = '002', runs = subset(cases, PS==ps)$run, unzip_py_files = F)
  generate_sbatch(ps, scen = 'G', experiment = '003', runs = subset(cases, PS==ps)$run, unzip_py_files = F)
  generate_sbatch(ps, scen = 'G', experiment = '004', runs = subset(cases, PS==ps)$run, unzip_py_files = F)
}

# General stuff -----------------------------------------------------------

# Extract the random seeds of all experiments and runs
random_seeds <- c()
for (ps in ps_range){
  for (e in exp_range){
    for(r in 1:5){
      x <- data.frame(ps,e,r,seed=get_random_seed(PS = ps, scenario = 'N', run_range = r))
      random_seeds <- rbind(random_seeds, x)
    }
  }
}
y <- random_seeds %>% group_by(ps,r,seed) %>% summarise(s=length(seed))


for (ps in sprintf('%0.2d', 27:31)){
  for (e in sprintf('%0.3d', 1:4)){
    cat(paste('sbatch PS',ps,'NE',e,'.sbatch',sep=''));cat('\n')
  }
}



design_irs <- create_intervention_scheme_IRS(PS_benchmark = '09', scenario_benchmark = 'S', IRS_START_TIMES = '29160', immigration_range=c(0,0.001), length_range=c(1800,3600), coverage_range=0.9, write_to_file = T, design_ref=design)



# plots -------------------------------------------------------------------
setwd('~/Documents/malaria_interventions/')
seasonality <- read_csv('mosquito_population_seasonality.csv', col_names = c('day','num_mosquitos'))
seasonality$grp <- 'S'
seasonality2 <- NULL
for (i in 1:10){
  seasonality2 <- rbind(seasonality2,seasonality)
}
seasonality=seasonality2
seasonality$day <- 1:3600
irs01 <- read_csv('IRS01.csv', col_names = c('day','num_mosquitos'))
irs01$grp <- 'IRS01'
x <- rbind(seasonality,irs01)
x %>% ggplot(aes(day, num_mosquitos, color=grp))+geom_line()

