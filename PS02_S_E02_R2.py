RANDOM_SEED=868381
N_GENES_INITIAL=12000
N_GENES_PER_STRAIN=60
N_LOCI=2
N_ALLELES_INITIAL=N_LOCI*[1200]
SELECTION_MODE='SPECIFIC_IMMUNITY'
N_INFECTIONS_FOR_GENERAL_IMMUNITY=0
GENERAL_IMMUNITY_PARAMS=[-1, -1, -1, -1]
CLEARANCE_RATE_IMMUNE=0.041628833086078301 # For generalized immunity
TRANSITION_RATE_NOT_IMMUNE=1.0 / 6.0 # switching rate. every 6 days
TRANSITION_RATE_IMMUNE=1000
PRINT_FUNCTION_TRACE=False
PRINT_DEBUG_LEVEL=0
T_YEAR=360.0
T_BURNIN=2880
T_END=36000
SAMPLE_DB_FILENAME='"PS02_S_E02_R2.sqlite"'
PRINT_INFO_PERIOD=120.0
VERIFICATION_ON=False
VERIFICATION_PERIOD=36000
SAVE_TO_CHECKPOINT=False
CHECKPOINT_SAVE_FILENAME='""'
CHECKPOINT_SAVE_PERIOD=0
LOAD_FROM_CHECKPOINT=True
CHECKPOINT_LOAD_FILENAME='"PS02_S_E00_R2_CP.sqlite"'
OUTPUT_HOSTS=True
OUTPUT_GENES=False # Set to false to save space
OUTPUT_STRAINS=False # Set to false to save space
HOST_SAMPLING_ON=True
HOST_SAMPLING_PERIOD=30.0
HOST_SAMPLE_SIZE=100
GENE_TRANSMISSIBILITY=0.5
COINFECTION_REDUCES_TRANSMISSION=True
ECTOPIC_RECOMBINATION_RATE=1.8e-07
P_ECTOPIC_RECOMBINATION_IS_CONVERSION=0
IMMUNITY_LOSS_RATE=0.001
MUTATION_RATE=1.42e-08
T_LIVER_STAGE=14.0
MEAN_HOST_LIFETIME=30.0 * T_YEAR
MAX_HOST_LIFETIME=80.0 * T_YEAR
N_POPULATIONS=1
N_HOSTS=[10000]
N_INITIAL_INFECTIONS=[20]
BITING_RATE_MEAN=[5e-05]
BITING_RATE_RELATIVE_AMPLITUDE=[0.0]
BITING_RATE_PEAK_PHASE=[0.0]
DAILY_BITING_RATE_DISTRIBUTION=[184.982060636238,200.255404185268,190.091657432075,176.936673923227,164.343031052803,152.793760078616,142.280309311918,132.721990757893,124.033625709693,116.13590946434,108.956491569867,102.429591067436,96.4954777083097,91.0998775063187,86.1935081934926,81.7316148582299,77.6735473752113,73.9823941296359,70.624631409777,67.569823745941,64.7903263402798,62.261040675501,59.9591757095815,57.8640360180417,55.9568281889751,54.2204871570074,52.6395139643939,51.1998352707491,49.8886665452749,48.6943983409682,47.6064821131221,46.6153343553555,45.7122530802748,44.889322061416,44.13935303385,43.4558081247075,42.8327483785961,42.2647680297853,41.7469538312647,41.2748368576276,40.8443501562017,40.4517970480741,40.0938070759306,39.7673175268291,39.4695383049637,39.1979296104896,38.9501783120656,38.7241775187336,38.5180090184007,38.3299237028482,38.1583293639103,38.0017737769016,37.8589341338113,37.7286047949206,37.609686843514,37.50117832764,37.4021653181591,37.3118147288408,37.2293672736381,37.1541298330001,37.0854713133807,37.0228146684226,36.9656352249567,36.9134526965649,36.8658304874747,36.8223689332749,36.7827047634528,36.7465053716365,36.7134680370484,36.683316586589,36.6557984022823,36.6306826034548,36.6077597837153,36.586838896635,36.5677447120854,36.5503175349244,36.5344117986761,36.5198946799105,36.5066448234627,36.4945517554252,36.4835144475775,36.4734406657602,36.4642461883492,36.4558541345434,36.4481944766219,36.4412033547379,36.4348224858472,36.4289985845382,36.423682949872,36.4188312037933,36.4144029157693,36.4103610895911,36.4066719764741,36.4033048155977,36.4002316266374,36.3974265546817,36.3948664137349,36.3925295614191,36.390396747356,36.3884499478837,36.3866731020119,36.3850512618356,36.383570977466,36.3822198706082,36.3809866531596,36.3798611099218,36.3788337913255,36.3778961228354,36.3770403186877,36.376259127581,36.3755461792137,36.3748953467457,36.3743013727248,36.3737591652055,36.3732643065013,36.3728126044509,36.3724003301053,36.3720240316874,36.3716805650306,36.3713671195684,36.3710809645286,36.3708198308674,36.3705814728698,36.3703638797406,36.370165316919,36.3699840324805,36.369818578311,36.3696675718257,36.3695297102287,36.3694039040808,36.3692890594064,36.3691842391427,36.3690885778915,36.3690012449681,36.368921567071,36.3688488589252,36.3687824611981,36.3687218882497,36.3686666214112,36.368616137142,36.368570066055,36.3685280375853,36.3684896516034,36.3684545980557,36.3684226195349,36.3683934248545,36.3683667563531,36.3683424240487,36.3683202157831,36.3682999305956,36.368281422106,36.3682645324877,36.368249105087,36.3682350279662,36.3682221849002,36.3682104538153,36.3681997446467,36.3681899755132,36.3681810545177,36.3681729058745,36.3681654734027,36.3681586903366,36.3681524946063,36.3681468412579,36.3681416820438,36.3681369698813,36.3681326707489,36.368128748638,36.3681251665486,36.3681218969029,36.3681189135428,36.3681161882129,36.3681137006326,36.3681114340905,36.3681093667845,36.3681074769326,36.36810574951,36.3681041736667,36.3681027356301,36.3681014214769,36.3824169711611,36.5868636264479,37.3137292989518,38.8896588849721,41.5969269867187,45.6776517011925,51.3462962768688,58.8024246075306,68.2404532572648,79.85474880501,93.8398405839524,110.386742354164,129.677078750806,151.876740164501,177.130350470818,205.557250316746,237.249171883638,272.269439482822,310.65335757952,352.409425053749,397.521022857601,445.948360485904,497.63044888562,552.487019808868,610.420329955446,671.316794398049,735.04845978698,801.474377036833,870.441800378375,941.787281377538,1015.33771227084,1090.91125222515,1168.31823275306,1247.36198486262,1327.83965484193,1409.54300728571,1492.25915115092,1575.77130548484,1659.8595289313,1744.30143733407,1828.87291466673,1913.34881406379,1997.50365029742,2081.11228344749,2163.95059344204,2245.79614473898,2326.42884041966,2405.6315648126,2483.19081324801,2558.89730838636,2632.54660134387,2703.93965664881,2772.88342010047,2839.19136746662,2902.68403370367,2963.18952097676,3020.54398416973,3074.59209310879,3125.18747033453,3172.19310299439,3215.48172820662,3254.93619083601,3290.44977281168,3321.92638760145,3349.2811705967,3372.44043808459,3391.34189360059,3405.93488836559,3416.18054736254,3422.05188380067,3423.53387682192,3420.623514707,3413.3298065236,3401.67376335901,3385.68834152934,3365.41835574442,3340.92036440204,3312.26251407319,3279.524351653,3242.79662352343,3202.18101552305,3157.78988588425,3109.7459624656,3058.18199001872,3003.24039148679,2945.07291252369,2883.84014358727,2819.71112706063,2752.86287733781,2683.47992901588,2611.75372944546,2537.88232202856,2462.06973786712,2384.52534631951,2305.46335094028,2225.10224033472,2143.66415945827,2061.37433812571,1978.46049021498,1895.15220628465,1811.68034630028,1728.27643230333,1645.17204331725,1562.59821799246,1480.78485370311,1399.96012162063,1320.34992419284,1242.1773279901,1165.66198459152,1091.01964102924,1018.46169386935,948.194705729685,880.41998007034,815.333164886091,753.1238763191,693.975388265441,638.06417359528,585.55935270184,536.621744821678,491.402242115739,450.01571640615,412.35452014996,378.127761460069,347.029170321546,318.7736375522,293.101017375666,269.774799848469,248.580096363242,229.321704355466,211.822250033048,195.920630103989,181.4704890447,168.338856724313,156.404958764567,145.559075293999,135.701539412099,126.741814829866,118.59766678999,111.194393294844,104.464143350716,98.3452956785572,92.7818696933424,87.7230381579216,83.1226257338353,78.938706103822,75.1332009977282,71.6715344021282,68.5223056833404,65.6570004450704,63.049727737086,60.676978167339,58.5174022114277,56.5516126906658,54.7620077186872,53.1325998586964,51.6488713786716,50.2976368482979,49.0669189753349,47.945843498837,46.9245231865883,45.9939792650736,45.1460537475325,44.3733299143677,43.6690703446293,43.0271456468427,42.4419851037608,41.9085232715463,41.4221516208221,40.9786797313279,40.5742921640239,40.2055199166146,39.8692027697277,39.562466316178,39.2826923331566,39.0274969459958,38.7947100974123,38.5823535152005,38.3886258017695,38.2118850795468,38.050636299575]
IRS_ON=False
IRS_START_TIMES=[100,200]
BITING_RATE_FACTORS=[[0.05,0.06,0.07],[0.05,0.06,0.07]] #each vector is the "raw" values of mosquito numbers from mathematica. The ABM will use this vector to multiply the BITING_RATE_MEAN until it ends, and then go back to the baseline DAILY_BITING_RATE_DISTRIBUTION
IRS_IMMIGRATION_RATE_FACTORS=[0.2,0.6]
MDA_ON=False
MDA_START_TIMES=[250,300]
HOST_FAIL_RATE=[0.2,0.2] # % of host sthat did not take the drug
DRUG_EFF_DURATION=[15,15] # How long the drug remains effective in the body
MDA_IMMIGRATION_RATE_FACTORS=[0.2,0.6]
IMMIGRATION_ON=True
IMMIGRATION_RATE=[1.0]
P_IMMIGRATION_INCLUDES_NEW_GENES=0.5
N_IMMIGRATION_NEW_GENES=0
