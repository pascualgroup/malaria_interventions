library(tidyverse)
library(magrittr)
library(sqldf)
library(rPython)
library(googlesheets)


source('~/Documents/malaria_interventions/functions.R')

# Create parameter and job files -------------------------------------------

setwd('/home/shai/Documents/malaria_interventions')
sqlite_path_global <- '/media/Data/PLOS_Biol/sqlite_S'
parameter_files_path_global <- '/media/Data/PLOS_Biol/parameter_files'

# Clear previous files if necessary
# clear_previous_files(run = 6, exclude_sqlite = F, exclude_CP = F, exclude_control = F, test = T)
# Get data design 
# design <- loadExperiments_GoogleSheets(local = T, workBookName = 'PLOS_Biol_design.csv') 
design <- loadExperiments_GoogleSheets(local = F, workBookName = 'PLOS_Biol_design', sheetID = 2) 

# Create the reference experiments (checkpoint and control)
ps_range <- sprintf('%0.2d', 1)
exp_range <- sprintf('%0.3d', 0:1)
run_range <- 1
work_scenario <- 'G'
# Generate 000 and 001 experiments
design_subset <- subset(design, PS %in% ps_range & scenario==work_scenario)
generate_files(row_range = 1:nrow(design_subset), run_range = run_range, experimental_design = design_subset, 
               biting_rate_file = 'fixed_biting_rates_1.csv', 
               target_folder = '/media/Data/PLOS_Biol/parameter_files')

# # If checkpoints already exist, can create the reference experiments (control only)
for (ps in ps_range){
  print(ps)
  design_control <- subset(design, PS==ps & scenario == work_scenario & exp=='001')
  # clear_previous_files(parameter_space = ps, scenario = work_scenario, exclude_sqlite = F, exclude_CP = T, exclude_control = F, test = F)
  seeds <- get_random_seed(PS = ps, scenario = work_scenario, run_range = run_range)
  generate_files(row_range = 1, run_range = run_range, experimental_design = design_control, random_seed = seeds)
}

# Create IRS experiments  
for (ps in ps_range){
  design_irs <- create_intervention_scheme_IRS(PS_benchmark = ps, scenario_benchmark = work_scenario, IRS_START_TIMES = '29160', immigration_range=c(0), length_range=c(720,1800,3600), coverage_range=0.9, poolsize_bounce = 'False', write_to_file = F, design_ref=design)
  generate_files(row_range = 1:nrow(design_irs), run_range = run_range, 
                 random_seed = get_random_seed(ps, work_scenario, run_range = run_range), 
                 experimental_design=design_irs)
}

setwd('/media/Data/PLOS_Biol/parameter_files/')
# ZIP all the PY and sbatch files
unlink('files_to_run.txt')
sink('files_to_run.txt', append = T)
# Add sbatch files
files <- list.files(path = '/media/Data/PLOS_Biol/parameter_files', pattern = 'sbatch', full.names = F) 
files <- files[str_detect(string = files, paste(work_scenario,'E',sep=''))]
files <- files[str_sub(files,3,4) %in% ps_range]
files <- files[str_sub(files,7,9) %in% exp_range]
for (i in 1:length(files)){
  cat(files[i]);cat('\n')
}
# Add py files
files <- list.files(path = '/media/Data/PLOS_Biol/parameter_files', pattern = '.py', full.names = F) 
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
files_S <- list.files(path = '/media/Data/PLOS_Biol/sqlite_S/', pattern = '\\.sqlite', full.names = F)
files_G <- list.files(path = '/media/Data/PLOS_Biol/sqlite_G/', pattern = '\\.sqlite', full.names = F)
files <- c(files_S, files_G)
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
  db <- dbConnect(SQLite(), dbname = paste('/media/Data/malaria_interventions_data/',f,sep=''))# see if folder is correct
  x <- dbGetQuery(db, 'SELECT max(time) FROM summary')
  dbDisconnect(db)
  x
}))
unique(files_sqlite$run_time)
files_sqlite %>% filter(is.na(run_time))


files_sqlite %<>% mutate(scenario=factor(scenario, levels=c('S','G')))  %>%
  arrange(CP, PS, scenario, experiment, run)         
files_sqlite %<>% filter(is.na(CP) & as.numeric(PS)>=27) %>% 
  # filter(scenario=='S') %>% 
  group_by(PS,scenario,experiment) %>% 
  summarise(runs_completed=length(run)) 
files_sqlite %>% 
  filter(runs_completed<50) %>% print(n = Inf)
files_sqlite %>% group_by(scenario, PS, experiment) %>% summarise(x=sum(runs_completed)) %>% print(n = Inf)
files_sqlite %>% filter(scenario=='G') %>% print(n = Inf)


# f <- read_csv('sqlite_files.txt', col_names = F)
# f$PS = sapply(str_split(f$X1,'_'),function (x) parse_number(x[1]))
# f$scenario=sapply(str_split(f$X1,'_'),function (x) x[2])
# f$experiment=sapply(str_split(f$X1,'_'),function (x) str_sub(x[3],2,4))
# f$run= sapply(str_split(f$X1,'_'),function (x) parse_number(x[4]))
# f$PS <- sprintf('%0.2d', f$PS)
# 
# f %>% select(PS,scenario,experiment,run) %>% full_join(files_sqlite, by=c('PS','scenario','experiment','run')) %>% select(PS,scenario,experiment,run)
# 
# missing=files_sqlite %>% anti_join(f, by=c('PS','scenario','experiment','run')) %>% select(PS,scenario,experiment,run,file_sqlite) %>% print(n=Inf)
# for (fl in missing$file_sqlite){
#   file.rename(paste('~/Documents/malaria_interventions_data/sqlite_G/',fl,sep=''), paste('~/Documents/malaria_interventions_data/sqlite_G/tocopy/',fl,sep=''))
# }





#--- PY files ---
files <- list.files(path = '/media/Data/malaria_interventions_data/parameter_files/', pattern = 'py.zip', full.names = T)
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
files_py %>% filter(as.numeric(PS)>=27) %>% group_by(PS,scenario,experiment) %>% summarise(runs_count=length(run))
files_py %>% filter(scenario=='G') %>% filter(as.numeric(PS)>=27) %>% distinct(PS,run)


# Generate missing GI parameter files ------------------------------------- 
# For those runs that produced an error when fitting the curve, take curve-fit
# information from a random run where fit was successful.

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
existing_py_G <- files_py %>% filter(scenario=='G') %>% filter(as.numeric(PS)>=27 & experiment=='000') %>% distinct(PS,run)
# A list of existing and missing runs
py_files_G <- vector(mode = 'list', length = length(ps_range))
names(py_files_G) <- ps_range
for (ps in ps_range){
  existing <- subset(existing_py_G, PS==ps)$run
  py_files_G[[which(ps_range==ps)]]$existing <- existing
  py_files_G[[which(ps_range==ps)]]$missing <- setdiff(1:50,existing)
}


for (ps in ps_range){
  py_files_G_ps <- py_files_G[[which(ps_range==ps)]]
  # Create 000 files for all runs
  for (missing_run in py_files_G_ps$missing){
    # RANDOMLY select EXISTING GI run from which to take the fitting parameters
    existing_run <- sample(py_files_G_ps$existing,1)
    existing_file <- paste('PS',ps,'_G_E000_R',existing_run,'.py',sep='')
    unzip(paste(ps,'_G_py.zip',sep=''), existing_file)
    x <- readLines(existing_file)
    params_GI <- get_generalized_immunity(ps, existing_run)
    # Get the seed from the corresponding selection scenario
    unzip(paste(ps,'_S_py.zip',sep=''), paste('PS',ps,'_S_E000_R',missing_run,'.py',sep=''))
    seed <- get_random_seed(PS = ps, scenario = 'S', experiment = '000', run_range = missing_run)
    unlink(paste('PS',ps,'_S_E000_R',missing_run,'.py',sep=''))
    # Use the parameters to generate the file
    design_subset <- subset(design, PS %in% ps & scenario=='G' & exp=='000')
    generate_files(row_range = 1, random_seed = seed, run_range = missing_run, experimental_design = design_subset, params_GI = params_GI)
    # Remove the existing file
    unlink(existing_file)
  }
  
  # Get the seeds used for the missing runs to generate the rest of the experiments
  seeds <- get_random_seed(PS = ps, scenario = 'G', run_range = py_files_G_ps$missing)
  
  # Generate control experimets
  design_subset <- subset(design, PS %in% ps & scenario=='G' & exp=='001')
  generate_files(row_range = 1, run_range = py_files_G_ps$missing, experimental_design = design_subset, random_seed = seeds)
  
  # Generate IRS experiments
  design_irs <- create_intervention_scheme_IRS(PS_benchmark = ps, scenario_benchmark = work_scenario, IRS_START_TIMES = '29160', immigration_range=c(0), length_range=c(720,1800,3600), coverage_range=0.9, poolsize_bounce = 'False', write_to_file = F, design_ref=design)
  generate_files(row_range = 1:nrow(design_irs), run_range = py_files_G_ps$missing,
                 random_seed = seeds, 
                 experimental_design=design_irs)
  
  # Generate sbatch files
  generate_sbatch(ps = ps, scen = 'G', runs = py_files_G_ps$missing, experiment = '000', unzip_py_files = F)
  generate_sbatch(ps = ps, scen = 'G', runs = py_files_G_ps$missing, experiment = '001', unzip_py_files = F)
  generate_sbatch(ps = ps, scen = 'G', runs = py_files_G_ps$missing, experiment = '002', unzip_py_files = F)
  generate_sbatch(ps = ps, scen = 'G', runs = py_files_G_ps$missing, experiment = '003', unzip_py_files = F)
  generate_sbatch(ps = ps, scen = 'G', runs = py_files_G_ps$missing, experiment = '004', unzip_py_files = F)
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

ps_scen_combinations <- expand.grid(PS=27:39,scenario='G')

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



# Generate sbatch files to extract data -----------------------------------
sbatch_arguments <- expand.grid(PS=sprintf('%0.2d', 4:6),
                                scen=c('S','N','G'), 
                                array='1', 
                                stringsAsFactors = F)
sbatch_arguments$cutoff_prob <- rep(c(0.25,0.7,0.9),3)
sbatch_arguments$mem_per_cpu <- c(4000,4000,8000,12000,16000,32000,4000,6000,10000)
sbatch_arguments$time <- c('02:00:00','02:00:00','02:00:00','10:00:00','10:00:00','10:00:00','02:00:00','02:00:00','02:00:00')
# sbatch_arguments$mem_per_cpu <- c(rep(4000,2),rep(8000,2),rep(12000,3),rep(16000,2),rep(32000,4),
#                                   rep(8000,2),rep(16000,2),rep(32000,8),64000)
# sbatch_arguments$time <- c(rep('01:00:00',4),rep('03:00:00',4),rep('06:00:00',5),
#                            rep('04:00:00',4),rep('08:00:00',4),rep('12:00:00',5)) # For G

ps_range <- sprintf('%0.2d', 4:6)
for (scenario in c('S','N','G')){
  for (ps in ps_range){
    x <- readLines('~/Documents/malaria_interventions/get_data_midway_plosbiol.sbatch')
    str_sub(x[2],22,24) <- paste(ps,scenario,sep='')
    str_sub(x[3],16,18) <- subset(sbatch_arguments, PS==ps & scen==scenario)$time
    str_sub(x[4],33,35) <- paste(ps,scenario,sep='')
    str_sub(x[5],32,34) <- paste(ps,scenario,sep='')
    str_sub(x[6],17,20) <- subset(sbatch_arguments, PS==ps & scen==scenario)$array
    str_sub(x[9],23,25) <- subset(sbatch_arguments, PS==ps & scen==scenario)$mem_per_cpu
    str_sub(x[19],5,7) <- ps
    str_sub(x[20],11,13) <- scenario
    str_sub(x[22],13,16) <- subset(sbatch_arguments, PS==ps & scen==scenario)$cutoff_prob
    writeLines(x, paste(parameter_files_path_global,'/','PS',ps,'_',scenario,'_get_data_midway.sbatch',sep=''))
  }
}
# for i in {27..39}; do sbatch 'PS'$i'_G_get_data_midway.sbatch'; done


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

