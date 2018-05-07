RANDOM_SEED = 1

N_GENES_INITIAL = 12000 # pool
N_GENES_PER_STRAIN = 60
N_LOCI = 2
N_ALLELES_INITIAL = N_LOCI * [1200]

SELECTION_MODE = 'SPECIFIC_IMMUNITY'
N_INFECTIONS_FOR_GENERAL_IMMUNITY = 0
GENERAL_IMMUNITY_PARAMS = [-1, -1, -1, -1]
CLEARANCE_RATE_IMMUNE = 0.041628833086078301 # For generalized immunity

TRANSITION_RATE_NOT_IMMUNE = 1.0 / 6.0 # switching rate. every 6 days
TRANSITION_RATE_IMMUNE = 1000

PRINT_FUNCTION_TRACE = False
PRINT_DEBUG_LEVEL = 0

T_YEAR = 360.0
T_BURNIN = 0 # in days
T_END = 50.0 * T_YEAR

SAMPLE_DB_FILENAME = '"test_02.sqlite"'

PRINT_INFO_PERIOD = 120.0

VERIFICATION_ON = True
VERIFICATION_PERIOD = 50.0 * T_YEAR

SAVE_TO_CHECKPOINT = False
CHECKPOINT_SAVE_FILENAME = '""'
CHECKPOINT_SAVE_PERIOD = T_YEAR * 20

LOAD_FROM_CHECKPOINT = False
CHECKPOINT_LOAD_FILENAME = '""'

OUTPUT_HOSTS = True
OUTPUT_GENES = False # Set to false to save space
OUTPUT_STRAINS = False # Set to false to save space

HOST_SAMPLING_ON = True
HOST_SAMPLING_PERIOD = 30.0
HOST_SAMPLE_SIZE = 100

GENE_TRANSMISSIBILITY = 0.5
COINFECTION_REDUCES_TRANSMISSION = True

ECTOPIC_RECOMBINATION_RATE = 1.8e-07
P_ECTOPIC_RECOMBINATION_IS_CONVERSION = 0

IMMUNITY_LOSS_RATE = 0.001

MUTATION_RATE = 1.42e-08

T_LIVER_STAGE = 14.0


MEAN_HOST_LIFETIME = 30.0 * T_YEAR
MAX_HOST_LIFETIME = 80.0 * T_YEAR

N_POPULATIONS = 1

N_HOSTS = [10000]
N_INITIAL_INFECTIONS = [20]

BITING_RATE_MEAN = [1]
BITING_RATE_RELATIVE_AMPLITUDE = [0.0]
BITING_RATE_PEAK_PHASE = [0.0]
DAILY_BITING_RATE_DISTRIBUTION = [0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5] # the "raw" values of mosquito numbers from mathematica. needs to be of length 360. Otherwise will run without it

IRS_ON = False
IRS_START_TIMES = [100,200]
BITING_RATE_FACTORS = [[0.05,0.06,0.07],[0.05,0.06,0.07]] #each vector is the "raw" values of mosquito numbers from mathematica. The ABM will use this vector to multiply the BITING_RATE_MEAN until it ends, and then go back to the baseline DAILY_BITING_RATE_DISTRIBUTION
IRS_IMMIGRATION_RATE_FACTORS = [0.2,0.6]

MDA_ON = False
MDA_START_TIMES = [250,300]
HOST_FAIL_RATE = [0.2,0.2] # % of host sthat did not take the drug
DRUG_EFF_DURATION = [15,15] # How long the drug remains effective in the body
MDA_IMMIGRATION_RATE_FACTORS = [0.2,0.6]


IMMIGRATION_ON = True
IMMIGRATION_RATE = [1.0]
P_IMMIGRATION_INCLUDES_NEW_GENES = 0.5
N_IMMIGRATION_NEW_GENES = 0
