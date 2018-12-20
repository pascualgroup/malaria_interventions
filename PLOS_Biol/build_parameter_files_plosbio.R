source('~/Documents/malaria_interventions/functions.R')
prep.packages(c('tidyverse','magrittr','sqldf','rPython','googlesheets'))


# Create parameter and job files -------------------------------------------
if (detect_locale()=='Lab'){
  setwd('/home/shai/Documents/malaria_interventions')
  sqlite_path_global <- '/media/Data/PLOS_Biol/sqlite_S'
  parameter_files_path_global <- '/media/Data/PLOS_Biol/parameter_files'
}
if (detect_locale()=='Mac'){
  setwd('~/GitHub/malaria_interventions')
  sqlite_path_global <- '~/GitHub/PLOS_Biol/sqlite_S'
  parameter_files_path_global <- '~/GitHub/PLOS_Biol/parameter_files'
}

# Clear previous files if necessary
# clear_previous_files(run = 6, exclude_sqlite = F, exclude_CP = F, exclude_control = F, test = T)
# Get data design 
# design <- loadExperiments_GoogleSheets(local = T, workBookName = 'PLOS_Biol_design.csv') 
design <- loadExperiments_GoogleSheets(local = F, workBookName = 'PLOS_Biol_design', sheetID = 2) 

# Create the reference experiments (checkpoint and control)
ps_range <- sprintf('%03d', 500:599)
exp_range <- sprintf('%0.3d', 1)
run_range <- 1
work_scenario <- 'S'
# Generate 000 and 001 experiments
design_subset <- subset(design, PS %in% ps_range & scenario==work_scenario & exp %in% exp_range)
generate_files(row_range = 1:nrow(design_subset), run_range = run_range, 
               experimental_design = design_subset, 
               # The radom_seed is necessary if using CP for intervention when
               # not creating the intervention experiment parameter file at the
               # same time as the 000 file.
               #random_seed = get_random_seed(PS = 18, scenario = work_scenario, run_range = run_range, folder = '/media/Data/PLOS_Biol/parameter_files/'),
               biting_rate_file = design_subset$DAILY_BITING_RATE_DISTRIBUTION[1],
               target_folder = parameter_files_path_global)

# # If checkpoints already exist, can create the reference experiments (control only)
# for (ps in ps_range){
#   print(ps)
#   design_control <- subset(design, PS==ps & scenario == work_scenario & exp=='001')
#   # clear_previous_files(parameter_space = ps, scenario = work_scenario, exclude_sqlite = F, exclude_CP = T, exclude_control = F, test = F)
#   seeds <- get_random_seed(PS = ps, scenario = work_scenario, run_range = run_range)
#   generate_files(row_range = 1, run_range = run_range, experimental_design = design_control, random_seed = seeds)
# }


setwd(parameter_files_path_global)
# ZIP all the PY and sbatch files
unlink('files_to_run.txt')
sink('files_to_run.txt', append = T)
# Add sbatch files
files <- list.files(path = parameter_files_path_global, pattern = 'sbatch', full.names = F) 
files <- files[str_detect(string = files, paste(work_scenario,'E',sep=''))]
files <- files[str_sub(files,3,5) %in% ps_range]
files <- files[str_sub(files,8,10) %in% exp_range]
for (i in 1:length(files)){
  cat(files[i]);cat('\n')
}
# Add py files
files <- list.files(path = parameter_files_path_global, pattern = '.py', full.names = F) 
files <- files[str_detect(string = files, paste('_',work_scenario,sep=''))]
files <- files[str_sub(files,3,5) %in% ps_range]
files <- files[str_sub(files,10,12) %in% exp_range]
files <- files[parse_number(str_sub(files,14,17)) %in% run_range]
for (i in 1:length(files)){
  cat(files[i]);cat('\n')
}
sink.reset()

# Create zip
unlink('files_to_run.zip')
system('zip files_to_run.zip -@ < files_to_run.txt')

# Copy the file to Midway and unzip it, then run the run_E000.sh file for run the checkpoints:
sink('/media/Data/PLOS_Biol/parameter_files/run_E000.sh')
for (ps in ps_range){
  cat('sbatch PS',ps,work_scenario,'E000.sbatch',sep='');cat('\n')
}
sink.reset()

# Then run control and interventions
## Run in Midway terminal:
cat("rm job_ids.txt; sacct -u pilosofs --format=jobid,jobname --starttime 2018-12-18T11:13:00 --name=");cat(paste(paste(ps_range,work_scenario,'E000',sep=''),collapse = ','));cat(" >> 'job_ids.txt'")
## Copy file from Midway and run:
jobids <- read.table('job_ids.txt', header = F, skip=2) 
jobids <- na.omit(unique(parse_number(jobids$V1))) # 1 job id per PS
length(jobids)==length(ps_range)
exp_range='002'
sink('/media/Data/PLOS_Biol/parameter_files/run_E002.sh')
for (ps in ps_range){
  for (e in exp_range){
    cat(paste('sbatch -d afterok:',jobids[which(ps_range==ps)],' PS',ps,work_scenario,'E',e,'.sbatch',sep=''));cat('\n')
  }
}
sink.reset()




# Or, if checkpoints are already finished:
setwd('/media/Data/PLOS_Biol/parameter_files')
unlink('jobs_to_run.sh')
sink('jobs_to_run.sh')
for (ps in ps_range){
  for (e in exp_range[-1]){
    cat(paste('sbatch PS',ps,work_scenario,'E',e,'.sbatch',sep=''));cat('\n')
  }
}
sink.reset()

# Zip files to reduce the clutter -----------------------------------------
setwd('/media/Data/PLOS_Biol/parameter_files')
system('rm *.sbatch')

scenario_range <- c('S','N','G')
exp_range <- sprintf('%0.3d', 0:2)
ps_range <- sprintf('%0.3d', 500:599)

for (ps in ps_range){
  for (scenario in scenario_range){
    unlink('files_tmp.txt')
    sink('files_tmp.txt', append = T)
    # Add py files
    files <- list.files(path = '/media/Data/PLOS_Biol/parameter_files/', pattern = '\\.py', full.names = F) 
    if (str_length(ps)==2){
      files <- files[str_sub(files,3,4) %in% ps]
      files <- files[str_sub(files,6,6) %in% scenario]
      files <- files[str_sub(files,9,11) %in% exp_range]
    }
    if (str_length(ps)==3){
      files <- files[str_sub(files,3,5) %in% ps]
      files <- files[str_sub(files,7,7) %in% scenario]
      files <- files[str_sub(files,10,12) %in% exp_range]
    }
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
unlink('files_tmp.txt')


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
files_sqlite$PS <- sprintf('%0.3d', files_sqlite$PS)

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
files_py$PS <- sprintf('%0.3d', files_py$PS)

files_py %<>% filter(as.numeric(PS)>=27) %>% mutate(scenario=factor(scenario, levels=c('S','N','G'))) %>% arrange(scenario, PS, experiment, run)
files_py %>% filter(as.numeric(PS)>=27) %>% group_by(PS,scenario,experiment) %>% summarise(runs_count=length(run))
files_py %>% filter(scenario=='G') %>% filter(as.numeric(PS)>=27) %>% distinct(PS,run)


# Experiments for cutoff --------------------------------------------------

cutoff_design <- expand.grid(PS=sprintf('%0.2d', 5),
                             scenario=c('S','G','N'),
                             cutoff_prob=seq(0.25,0.95,0.05),
                             array='11-50',
                             calculate_mFst=F,
                             stringsAsFactors = F)
cutoff_design$job_name <- paste(cutoff_design$PS, cutoff_design$scenario,'_',cutoff_design$cutoff_prob,sep='')
cutoff_design$job_name <- gsub('\\.','',cutoff_design$job_name)

cutoff_design$mem_per_cpu <- 32000
for (i in 1:nrow(cutoff_design)){
  if(cutoff_design[i,'cutoff_prob']<0.35 & cutoff_design[i,'PS']=='06'){cutoff_design$mem_per_cpu[i] <- 58000}
  if(cutoff_design[i,'cutoff_prob']>0.8 & cutoff_design[i,'PS']!='06'){cutoff_design$mem_per_cpu[i] <- 16000}
}

cutoff_design$time <- '15:00:00'
for (i in 1:nrow(cutoff_design)){
  if(cutoff_design[i,'cutoff_prob']<0.35 & cutoff_design[i,'PS']=='06'){cutoff_design$time[i] <- '36:00:00'}
  if(cutoff_design[i,'cutoff_prob']>0.8 & cutoff_design[i,'PS']!='06'){cutoff_design$time[i] <- '10:00:00'}
}

#Use function make_sbatch_get_data()



# Generate sbatch files to extract data -----------------------------------

# Creating these files can be tricky. Pay carefult attention to combinations of
# PS, scenario, experiment and cutoffs in the resulting sbatch files!!!

sbatch_arguments <- expand.grid(PS=sprintf('%0.2d', 4:6),
                                scen=c('S','N','G'),
                                array='11-50', 
                                layers='1:300',
                                exp=c('001'),
                                stringsAsFactors = F)
sbatch_arguments$cutoff_prob <- rep(c(0.3,0.6,0.85),3)
sbatch_arguments$mem_per_cpu <- rep(c(6000,12000,32000),3)
sbatch_arguments$time <- rep(c('04:00:00','05:00:00','10:00:00'),3)

# sbatch_arguments <- subset(sbatch_arguments, scen=='G'&PS=='06')



# cal <- as.tibble(build_calendar(num_years = 25, 10))
# cal %>% filter(!is.na(survey)) %>% group_by(survey,layer) %>%
#   summarise(first_day=min(running_day),last_day=max(running_day), year=unique(year_sim)+2002, month=unique(month_sim))
sbatch_arguments <- expand.grid(PS=sprintf('%0.2d', 18),
                                scen=c('N','G'),
                                exp='002',
                                array='1',
                                layers='118,126,138,142,154,162',
                                # layers='1:300',
                                stringsAsFactors = F)
sbatch_arguments$cutoff_prob <- 0.85
sbatch_arguments$mem_per_cpu <- 8000
sbatch_arguments$time <- '01:00:00'

#Use function make_sbatch_get_data()


# Generate files for sensitivity analysis of main results ---------------------------------

# Generate 84 simulations for selection in high diversity with variation across 3 parameters:
# 1. Diversity 
# 2. Mean biting rate
# 3. Naive duration of infection (range 6-18 months)

# This code creates the design
diversity_range <- seq(10000,16000,1000) # Main analysis is 12000
biting_range <- seq(0.3,0.5,0.1) # Main analysis is 0.00040
naive_doi_range <- 1/(seq(180,540,120)/60) # Main analysis is 1/6
sensitivity_params <- expand.grid(N_GENES_INITIAL=diversity_range, BITING_RATE_MEAN=biting_range, TRANSITION_RATE_NOT_IMMUNE=naive_doi_range)
sensitivity_params <- as.tibble(sensitivity_params)
nrow(sensitivity_params)
sensitivity_params$PS <- str_pad((1:nrow(sensitivity_params))+99,width = 3, side = 'left', pad = '0')

design <- loadExperiments_GoogleSheets(local = F, workBookName = 'PLOS_Biol_design', sheetID = 2) 
design_seed_000 <- subset(design, PS=='06' & scenario=='S' & exp=='000')
design_seed_000 %<>% slice(rep(1:n(), each = nrow(sensitivity_params)))
design_seed_000$PS <- sensitivity_params$PS
design_seed_000$N_GENES_INITIAL <- sensitivity_params$N_GENES_INITIAL
design_seed_000$BITING_RATE_MEAN <- sensitivity_params$BITING_RATE_MEAN
design_seed_000$TRANSITION_RATE_NOT_IMMUNE <- sensitivity_params$TRANSITION_RATE_NOT_IMMUNE

design_seed_001 <- subset(design, PS=='06' & scenario=='S' & exp=='001')
design_seed_001 %<>% slice(rep(1:n(), each = nrow(sensitivity_params)))
design_seed_001$PS <- sensitivity_params$PS
design_seed_001$N_GENES_INITIAL <- sensitivity_params$N_GENES_INITIAL
design_seed_001$BITING_RATE_MEAN <- sensitivity_params$BITING_RATE_MEAN
design_seed_001$TRANSITION_RATE_NOT_IMMUNE <- sensitivity_params$TRANSITION_RATE_NOT_IMMUNE

design <- design_seed_000 %>% bind_rows(design_seed_001)
design$mem_per_cpu <- 32000
design$wall_time <- '20:00:00'

# Create files to get data after sqlites are created.
sbatch_arguments <- expand.grid(PS=sprintf('%0.3d', 100:183),
                                scen=c('S'),
                                array='1', 
                                layers='1:300',
                                exp=c('001'),
                                stringsAsFactors = F)
sbatch_arguments$cutoff_prob <-0.85
sbatch_arguments$mem_per_cpu <- 20000
sbatch_arguments$time <- '05:00:00'

#Use function make_sbatch_get_data()

# Generate files for empirical data comparison ----------------------------

# Generate 100 simulations per scenario for seasonality in high diversity with variation across 2 parameters:
# 1. Diversity 
# 2. Mean biting rate

# This code creates the design. Use that design to create the sqlite files.

diversity_range <- round(seq(30000,40000,length.out = 10)) # Main analysis is 35000
biting_range <- seq(0.0001,0.0003,length.out = 10) # Main analysis is 0.00020
empirical_comparison_params <- expand.grid(N_GENES_INITIAL=diversity_range, BITING_RATE_MEAN=biting_range)
empirical_comparison_params <- as.tibble(empirical_comparison_params)
nrow(empirical_comparison_params)
empirical_comparison_params$PS <- str_pad((1:nrow(empirical_comparison_params))+499,width = 3, side = 'left', pad = '0')

design <- loadExperiments_GoogleSheets(local = F, workBookName = 'PLOS_Biol_design', sheetID = 2) 
work_scenario <- 'G'

design_seed_000 <- subset(design, PS=='18' & scenario==work_scenario & exp=='000')
design_seed_000 %<>% slice(rep(1:n(), each = nrow(empirical_comparison_params)))
design_seed_000$PS <- empirical_comparison_params$PS
design_seed_000$N_GENES_INITIAL <- empirical_comparison_params$N_GENES_INITIAL
design_seed_000$BITING_RATE_MEAN <- empirical_comparison_params$BITING_RATE_MEAN

design_seed_001 <- subset(design, PS=='18' & scenario==work_scenario & exp=='001')
design_seed_001 %<>% slice(rep(1:n(), each = nrow(empirical_comparison_params)))
design_seed_001$PS <- empirical_comparison_params$PS
design_seed_001$N_GENES_INITIAL <- empirical_comparison_params$N_GENES_INITIAL
design_seed_001$BITING_RATE_MEAN <- empirical_comparison_params$BITING_RATE_MEAN

design_seed_002 <- subset(design, PS=='18' & scenario==work_scenario & exp=='002')
design_seed_002 %<>% slice(rep(1:n(), each = nrow(empirical_comparison_params)))
design_seed_002$PS <- empirical_comparison_params$PS
design_seed_002$N_GENES_INITIAL <- empirical_comparison_params$N_GENES_INITIAL
design_seed_002$BITING_RATE_MEAN <- empirical_comparison_params$BITING_RATE_MEAN

design <- design_seed_000 %>% bind_rows(design_seed_001) %>% bind_rows(design_seed_002)
design$mem_per_cpu <- 32000
design$wall_time <- '20:00:00'


if (detect_locale()=='Lab'){
  setwd('/home/shai/Documents/malaria_interventions')
  sqlite_path_global <- '/media/Data/PLOS_Biol/sqlite_S'
  parameter_files_path_global <- '/media/Data/PLOS_Biol/parameter_files'
}
if (detect_locale()=='Mac'){
  setwd('~/GitHub/malaria_interventions')
  sqlite_path_global <- '~/GitHub/PLOS_Biol/sqlite_S'
  parameter_files_path_global <- '~/GitHub/PLOS_Biol/parameter_files'
}


# Create the reference experiments (checkpoint and control)
ps_range <- sprintf('%03d', 500:599)
exp_range <- sprintf('%0.3d', 0:2)
run_range <- 1

for (ps in ps_range){
design_subset <- subset(design, PS %in% ps & scenario==work_scenario & exp %in% exp_range)
generate_files(row_range = 1:nrow(design_subset), run_range = run_range, 
               experimental_design = design_subset, 
               # The radom_seed is necessary if using CP for intervention when
               # not creating the intervention experiment parameter file at the
               # same time as the 000 file.
               #random_seed = get_random_seed(PS = 18, scenario = work_scenario, run_range = run_range, folder = '/media/Data/PLOS_Biol/parameter_files/'),
               biting_rate_file = design_subset$DAILY_BITING_RATE_DISTRIBUTION[1],
               target_folder = parameter_files_path_global)
}

# Copy the file to Midway and unzip it, then run the run_E000.sh file for run the checkpoints:
sink('/media/Data/PLOS_Biol/parameter_files/run_E000.sh')
for (ps in ps_range){
  cat('sbatch PS',ps,work_scenario,'E000.sbatch',sep='');cat('\n')
}
sink.reset()

jobids <- c(55805308:55805370,55805372:55805408) # from Midway after runnign run_E000.sh
sink('/media/Data/PLOS_Biol/parameter_files/run_E001.sh')
for (ps in ps_range){
  cat('sbatch -d afterok:',jobids[which(ps_range==ps)],' PS',ps,work_scenario,'E001.sbatch',sep='');cat('\n')
}
sink.reset()

sink('/media/Data/PLOS_Biol/parameter_files/run_E002.sh')
for (ps in ps_range){
  cat('sbatch -d afterok:',jobids[which(ps_range==ps)],' PS',ps,work_scenario,'E002.sbatch',sep='');cat('\n')
}
sink.reset()

# Create files to extract data
sbatch_arguments <- expand.grid(PS=sprintf('%0.3d', 500:599),
                                scen='N',
                                array='1', 
                                layers='1:100',
                                exp='002',
                                stringsAsFactors = F)
sbatch_arguments$cutoff_prob <- 0.85
sbatch_arguments$mem_per_cpu <- 20000
sbatch_arguments$time <- '02:00:00'
sbatch_arguments$after_job <- c(55803955:55804045,55804560,55804046,55804074,55804047:55804052)
# use make_sbatch_get_data here
make_sbatch_get_data(sbatch_arguments, make_networks = T)
make_sbatch_get_data(sbatch_arguments, repertoire_persistence = T) # For neutral simulations



#  Verify result file on Midway -------------------------------------------


# check for sqlite files
for i in '04' '05' '06'
do
ls sqlite/PS$i*_S_* -lv | wc -l
ls sqlite/PS$i*_N_* -lv | wc -l
ls sqlite/PS$i*_G_* -lv | wc -l
done
ls
# Check for result files
for i in 'sampled_alleles' 'sampled_infections' 'sampled_strains' 'summary_general' 'network_info' 'layer_summary' 'node_list' 'Infomap_multilayer.txt' 'Infomap_multilayer_expanded' 'modules' 'strain_sequences' 'unique_repertoires' 'temporal_diversity' 'mFst' 
do
echo $i
# ls Results/04_S/*_$i* -lv | wc -l
# ls Results/04_N/*_$i* -lv | wc -l
# ls Results/04_G/*_$i* -lv | wc -l
# ls Results/05_S/*_$i* -lv | wc -l
# ls Results/05_N/*_$i* -lv | wc -l
# ls Results/05_G/*_$i* -lv | wc -l
ls Results/06_S/*_$i* -lv | wc -l
ls Results/06_N/*_$i* -lv | wc -l
ls Results/06_G/*_$i* -lv | wc -l
done

# move result files to one folder
for i in 'S' 'G' 'N'
do
cd '/scratch/midway2/pilosofs/PLOS_Biol/Results/04_'$i
mv *_0.3_* ../cutoff_to_use/
cd '/scratch/midway2/pilosofs/PLOS_Biol/Results/05_'$i
mv *_0.6_* ../cutoff_to_use/
cd '/scratch/midway2/pilosofs/PLOS_Biol/Results/06_'$i
mv *_0.85_* ../cutoff_to_use/
done

# Check the cutoff folders for the files
for i in 'Infomap_multilayer' 'Infomap_multilayer_expanded' 'layer_summary' 'mFst' 'modules' 'network_info' 'node_list' 'sampled_alleles' 'sampled_infections' 'sampled_strains' 'strain_sequences' 'summary_general' 'temporal_diversity' 'unique_repertoires'
do
echo $i
  ls PS06_S*$i* -lv | wc -l
  ls PS06_G*$i* -lv | wc -l
  ls PS06_N*$i* -lv | wc -l
done


# Rename files and folders from 2-digit to 3-digit PS ---------------------------------
file_list <- list.files('/media/Data/PLOS_Biol/parameter_files/test', full.names = T, recursive = T, include.dirs = F)
for (f in file_list){
  folder <- str_sub(f,1,max(str_locate_all(f, '/')[[1]]))
  file_name <- str_sub(f,max(str_locate_all(f, '/')[[1]])+1,str_length(f))
  pos <- max(str_locate(file_name,'_'))
  pos_min <- ifelse(str_detect(file_name,'PS'),3,1)
  str_sub(file_name,pos_min,pos-1) <- str_pad(str_sub(file_name,pos_min,pos-1),3,'left','0')
  new_name <- paste(folder,file_name,sep='')
  file.rename(f,new_name)
}

file_list <- list.dirs('/media/Data/PLOS_Biol/parameter_files/test', full.names = T, recursive = T)
file_list <- file_list[2:length(file_list)]
for (f in file_list){
  folder <- str_sub(f,1,max(str_locate_all(f, '/')[[1]]))
  file_name <- str_sub(f,max(str_locate_all(f, '/')[[1]])+1,str_length(f))
  pos <- max(str_locate(file_name,'_'))
  pos_min <- ifelse(str_detect(file_name,'PS'),3,1)
  str_sub(file_name,pos_min,pos-1) <- str_pad(str_sub(file_name,pos_min,pos-1),3,'left','0')
  new_name <- paste(folder,file_name,sep='')
  system(paste('mv ',f,new_name))
}


# Calendar ----------------------------------------------------------------
cal <- as.tibble(build_calendar(num_years = 25, 10))
cal %>% filter(!is.na(survey)) %>% group_by(survey,layer) %>%
  summarise(first_day=min(running_day),last_day=max(running_day), year=unique(year_sim)+2002, month=unique(month_sim))

