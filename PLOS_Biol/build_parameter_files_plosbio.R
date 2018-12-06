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
ps_range <- sprintf('%0.2d', 14)
exp_range <- sprintf('%0.3d', 0:1)
run_range <- 1
work_scenario <- 'S'
# Generate 000 and 001 experiments
design_subset <- subset(design, PS %in% ps_range & scenario==work_scenario & exp %in% exp_range)
generate_files(row_range = 1:nrow(design_subset), run_range = run_range, 
               experimental_design = design_subset, 
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
files <- files[str_sub(files,3,4) %in% ps_range]
files <- files[str_sub(files,7,9) %in% exp_range]
for (i in 1:length(files)){
  cat(files[i]);cat('\n')
}
# Add py files
files <- list.files(path = parameter_files_path_global, pattern = '.py', full.names = F) 
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
exp_range <- sprintf('%0.3d', 0:1)
ps_range <- sprintf('%0.2d', 12:14)

for (ps in ps_range){
  for (scenario in scenario_range){
    unlink('files_tmp.txt')
    sink('files_tmp.txt', append = T)
    # Add py files
    files <- list.files(path = '/media/Data/PLOS_Biol/parameter_files/', pattern = '.py', full.names = F) 
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


# Experiments for cutoff --------------------------------------------------

cutoff_design <- expand.grid(PS=sprintf('%0.2d', 4:6),
                             scenario=c('S','G','N'),
                             cutoff_prob=seq(0.2,9,0.05),
                             array='4-10',
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

for (i in 1:nrow(cutoff_design)){
  x <- readLines('~/GitHub/malaria_interventions/get_data_midway_plosbiol.sbatch')
  ps <- cutoff_design$PS[i]
  scenario <- cutoff_design$scenario[i]
  cutoff_prob <- cutoff_design$cutoff_prob[i]
  x[2] <- paste('#SBATCH --job-name=',cutoff_design$job_name[i],sep='')
  str_sub(x[3],16,18) <- cutoff_design$time[i]
  x[4] <- paste('#SBATCH --output=slurm_output/',cutoff_design$job_name[i],'_%A_%a.out',sep='')
  x[5] <- paste('#SBATCH --error=slurm_output/',cutoff_design$job_name[i],'_%A_%a.err',sep='')
  str_sub(x[6],17,20) <- cutoff_design$array[i]
  str_sub(x[9],23,25) <- cutoff_design$mem_per_cpu[i]
  str_sub(x[19],5,7) <- ps
  str_sub(x[20],11,13) <- scenario
  str_sub(x[22],13,16) <- cutoff_prob
  writeLines(x, paste('~/GitHub/PLOS_Biol/Cutoff/','PS',ps,'_',scenario,'_',cutoff_prob,'_get_data_midway.sbatch',sep=''))
}

sink('~/GitHub/PLOS_Biol/Cutoff/run_cutoff_experiments.sh')
for (i in 1:nrow(cutoff_design)){
  ps <- cutoff_design$PS[i]
  scenario <- cutoff_design$scenario[i]
  cutoff_prob <- cutoff_design$cutoff_prob[i]
  cat('sbatch',paste('PS',ps,'_',scenario,'_',cutoff_prob,'_get_data_midway.sbatch',sep=''));cat('\n')
}
sink.reset()



# Generate sbatch files to extract data -----------------------------------
sbatch_arguments <- expand.grid(PS=sprintf('%0.2d', 4:6),
                                scen=c('S','N','G'),
                                array='6-50', 
                                stringsAsFactors = F)
sbatch_arguments$cutoff_prob <- rep(c(0.3,0.6,0.85),3)
sbatch_arguments$mem_per_cpu <- rep(c(6000,12000,32000),3)
sbatch_arguments$time <- rep(c('04:00:00','05:00:00','10:00:00'),3)

sbatch_arguments <- subset(sbatch_arguments, scen=='G'&PS=='06')

sbatch_arguments <- expand.grid(PS=sprintf('%0.2d', 11:12),
                                scen=c('S'),
                                array='1', 
                                stringsAsFactors = F)
sbatch_arguments$cutoff_prob <- 0.85
sbatch_arguments$mem_per_cpu <- 32000
sbatch_arguments$time <- '10:00:00'

for (scenario in unique(sbatch_arguments$scen)){
  for (ps in unique(sbatch_arguments$PS)){
    x <- readLines('~/Documents/malaria_interventions/PLOS_Biol/get_data_midway_plosbiol.sbatch')
    str_sub(x[2],22,24) <- paste(ps,scenario,sep='')
    str_sub(x[3],16,18) <- subset(sbatch_arguments, PS==ps & scen==scenario)$time
    str_sub(x[4],33,35) <- paste(ps,scenario,sep='')
    str_sub(x[5],32,34) <- paste(ps,scenario,sep='')
    str_sub(x[6],17,20) <- subset(sbatch_arguments, PS==ps & scen==scenario)$array
    str_sub(x[9],23,25) <- subset(sbatch_arguments, PS==ps & scen==scenario)$mem_per_cpu
    str_sub(x[19],5,7) <- ps
    str_sub(x[20],11,13) <- scenario
    str_sub(x[22],13,16) <- subset(sbatch_arguments, PS==ps & scen==scenario)$cutoff_prob
    writeLines(x, paste(parameter_files_path_global,'/','PS',ps,'_',scenario,'_',subset(sbatch_arguments, PS==ps & scen==scenario)$cutoff_prob,'_get_data_midway.sbatch',sep=''))
  }
}
# for i in {27..39}; do sbatch 'PS'$i'_G_get_data_midway.sbatch'; done
sink('/media/Data/PLOS_Biol/parameter_files/run_experiments.sh')
for (i in 1:nrow(sbatch_arguments)){
  ps <- sbatch_arguments$PS[i]
  scenario <- sbatch_arguments$scen[i]
  cutoff_prob <- sbatch_arguments$cutoff_prob[i]
  cat('sbatch ',paste('PS',ps,'_',scenario,'_',cutoff_prob,'_get_data_midway.sbatch',sep=''));cat('\n')
}
sink.reset()



# Generate files for sensitivity analysis ---------------------------------

# Generate 100 simulations per scenario in high diversity with variation across 3 parameters:
# 1. Diversity (vary between 10000 and 20000)
# 2. Mean biting rate (0.35-0.5)
# 3. Naive duration of infection (range 6-18 months)

diversity_range <- seq(10000,16000,1000) # Main analysis is 12000
biting_range <- seq(0.4,0.5,0.05) # Main analysis is 0.5
naive_doi_range <- 1/(seq(180,540,60)/60) # Main analysis is 1/6
sensitivity_params <- expand.grid(N_GENES_INITIAL=diversity_range, BITING_RATE_MEAN=biting_range, TRANSITION_RATE_NOT_IMMUNE=naive_doi_range)
nrow(sensitivity_params)


#  Verify result file on Midway -------------------------------------------


# check for sqlite files
for i in '04' '05' '06'
do
ls sqlite/PS$i*_S_* -lv | wc -l
ls sqlite/PS$i*_N_* -lv | wc -l
ls sqlite/PS$i*_G_* -lv | wc -l
done

# Check for result files
for i in 'sampled_alleles' 'sampled_infections' 'sampled_strains' 'summary_general' 'network_info' 'layer_summary' 'node_list' 'Infomap_multilayer' 'Infomap_multilayer_expanded' 'modules' 'strain_sequences' 'unique_repertoires' 'temporal_diversity' 'mFst' 
do
echo $i
# ls Results/04_S/*_0.3_$i* -lv | wc -l
# ls Results/04_N/*_0.3_$i* -lv | wc -l
# ls Results/04_G/*_0.3_$i* -lv | wc -l
# ls Results/05_S/*_0.6_$i* -lv | wc -l
# ls Results/05_N/*_0.6_$i* -lv | wc -l
# ls Results/05_G/*_0.6_$i* -lv | wc -l
# ls Results/06_S/*_0.85_$i* -lv | wc -l
ls Results/06_N/*_0.85_$i* -lv | wc -l
# ls Results/06_G/*_0.85_$i* -lv | wc -l
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
