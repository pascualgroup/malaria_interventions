# Initialize --------------------------------------------------------------
source('functions.R')
prep.packages(c('dplyr','readr','stringr','tidyr','tibble','magrittr','purrr','sqldf','rPython'))


if (length(commandArgs(trailingOnly=TRUE))==0) {
  args <- c('N','000','1:10')
} else {
  args <- commandArgs(trailingOnly=TRUE)
}
# job_ps <- sprintf('%0.2d',as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID')))
job_scenario <- as.character(args[1]) 
job_exp <- as.character(args[2]) 
job_run_range <- as.character(args[3]) 
job_run_range <- eval(parse(text=job_run_range))

# Create parameter and job files -------------------------------------------

# Get data design 
design <- readr::read_csv('malaria_interventions_design.csv') 
design$BITING_RATE_MEAN <- as.numeric(format(design$BITING_RATE_MEAN, scientific = FALSE))

if (job_exp=='000'){
  design_subset <- subset(design, PS == job_ps & scenario==job_scenario & exp==job_exp)
  if (job_scenario=='S'){
    generate_files(row_range = 1, run_range = job_run_range, random_seed = NULL, experimental_design = design_subset)
  } else {
    seeds <- get_random_seed(PS = job_ps, scenario = 'S', run_range = job_run_range)
    generate_files(row_range = 1, run_range = job_run_range, random_seed = seeds, experimental_design = design_subset)
  }
}

if (job_exp=='001'){
  design_subset <- subset(design, PS == job_ps & scenario==job_scenario & exp==job_exp)
  seeds <- get_random_seed(PS = job_ps, scenario = job_scenario, run_range = job_run_range)
  generate_files(row_range = 1, run_range = job_run_range, experimental_design = design_subset, random_seed = seeds)
}

# Create the corresponding IRS experiments. This scheme can change, or it can be
# moved to the google sheet for a different workflow in which experiments are
# specified directly instead of bening generated.
if (parse_number(job_exp)>1){
  design_irs <- create_intervention_scheme_IRS(PS_benchmark = job_ps, 
                                               scenario_benchmark = job_scenario,
                                               IRS_START_TIMES = '29160', 
                                               immigration_range=c(0), 
                                               length_range=c(720,1800,3600), 
                                               coverage_range=0.9, 
                                               poolsize_bounce = 'False', 
                                               write_to_file = F, 
                                               design_ref=design)
  seeds <- get_random_seed(job_ps, job_scenario, run_range = job_run_range)
  generate_files(row_range = 1:3, run_range = job_run_range, 
                 random_seed = seeds, 
                 experimental_design=design_irs)
}
