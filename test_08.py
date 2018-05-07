RANDOM_SEED=1
N_GENES_INITIAL=12000 # pool
N_GENES_PER_STRAIN=60
N_LOCI=2
N_ALLELES_INITIAL=N_LOCI * [1200]
SELECTION_MODE='SPECIFIC_IMMUNITY'
N_INFECTIONS_FOR_GENERAL_IMMUNITY=0
GENERAL_IMMUNITY_PARAMS=[-1, -1, -1, -1]
CLEARANCE_RATE_IMMUNE=0.041628833086078301 # For generalized immunity
TRANSITION_RATE_NOT_IMMUNE=1.0 / 6.0 # switching rate. every 6 days
TRANSITION_RATE_IMMUNE=1000
PRINT_FUNCTION_TRACE=False
PRINT_DEBUG_LEVEL=0
T_YEAR=360.0
T_BURNIN=0 # in days
T_END=50.0 * T_YEAR
SAMPLE_DB_FILENAME='"test_08.sqlite"'
PRINT_INFO_PERIOD=120.0
VERIFICATION_ON=True
VERIFICATION_PERIOD=50.0 * T_YEAR
SAVE_TO_CHECKPOINT=False
CHECKPOINT_SAVE_FILENAME='""'
CHECKPOINT_SAVE_PERIOD=T_YEAR * 20
LOAD_FROM_CHECKPOINT=False
CHECKPOINT_LOAD_FILENAME='""'
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
BITING_RATE_MEAN=[0.00005]
BITING_RATE_RELATIVE_AMPLITUDE=[0.0]
BITING_RATE_PEAK_PHASE=[0.0]
DAILY_BITING_RATE_DISTRIBUTION=[184.982060636238,200.255404185268,190.091657432075,176.936673923227,164.343031052803,152.793760078616,142.280309311918,132.721990757893,124.033625709693,116.13590946434,108.956491569867,102.429591067436,96.4954777083097,91.0998775063187,86.1935081934926,81.7316148582299,77.6735473752113,73.9823941296359,70.624631409777,67.569823745941,64.7903263402798,62.261040675501,59.9591757095815,57.8640360180417,55.9568281889751,54.2204871570074,52.6395139643939,51.1998352707491,49.8886665452749,48.6943983409682,47.6064821131221,46.6153343553555,45.7122530802748,44.889322061416,44.13935303385,43.4558081247075,42.8327483785961,42.2647680297853,41.7469538312647,41.2748368576276,40.8443501562017,40.4517970480741,40.0938070759306,39.7673175268291,39.4695383049637,39.1979296104896,38.9501783120656,38.7241775187336,38.5180090184007,38.3299237028482,38.1583293639103,38.0017737769016,37.8589341338113,37.7286047949206,37.609686843514,37.50117832764,37.4021653181591,37.3118147288408,37.2293672736381,37.1541298330001,37.0854713133807,37.0228146684226,36.9656352249567,36.9134526965649,36.8658304874747,36.8223689332749,36.7827047634528,36.7465053716365,36.7134680370484,36.683316586589,36.6557984022823,36.6306826034548,36.6077597837153,36.586838896635,36.5677447120854,36.5503175349244,36.5344117986761,36.5198946799105,36.5066448234627,36.4945517554252,36.4835144475775,36.4734406657602,36.4642461883492,36.4558541345434,36.4481944766219,36.4412033547379,36.4348224858472,36.4289985845382,36.423682949872,36.4188312037933,36.4144029157693,36.4103610895911,36.4066719764741,36.4033048155977,36.4002316266374,36.3974265546817,36.3948664137349,36.3925295614191,36.390396747356,36.3884499478837,36.3866731020119,36.3850512618356,36.383570977466,36.3822198706082,36.3809866531596,36.3798611099218,36.3788337913255,36.3778961228354,36.3770403186877,36.376259127581,36.3755461792137,36.3748953467457,36.3743013727248,36.3737591652055,36.3732643065013,36.3728126044509,36.3724003301053,36.3720240316874,36.3716805650306,36.3713671195684,36.3710809645286,36.3708198308674,36.3705814728698,36.3703638797406,36.370165316919,36.3699840324805,36.369818578311,36.3696675718257,36.3695297102287,36.3694039040808,36.3692890594064,36.3691842391427,36.3690885778915,36.3690012449681,36.368921567071,36.3688488589252,36.3687824611981,36.3687218882497,36.3686666214112,36.368616137142,36.368570066055,36.3685280375853,36.3684896516034,36.3684545980557,36.3684226195349,36.3683934248545,36.3683667563531,36.3683424240487,36.3683202157831,36.3682999373544,17.2685264724421,9.39143102937512,6.10950248045012,4.70717192297586,4.0747762400795,3.76073447724819,3.58240003975983,3.46624060717364,3.38241725884589,3.31828744574698,3.26788862282171,3.22788396299796,3.19605010611872,3.17072433788052,3.15059560995322,3.13461266998481,3.12193095566676,3.11187380618734,3.10390077244086,3.09758139889294,3.09257342902415,3.08860509388228,3.08546078301834,3.08296948220987,3.08099558676188,3.07943175165972,3.07819276441813,3.07721116038011,3.0764334351585,3.07581728827643,3.07874393310733,3.13264439369241,3.31887711564604,3.67986783658507,4.21967320565986,4.93157564143459,5.81541736584223,6.88233793543282,8.15306501128148,9.65455442368255,11.4170821612436,13.4721705348628,15.851067729003,18.5834486215355,21.6961549195511,25.2119488030417,29.1483498995162,33.5166609996754,38.3212851724897,43.5594092524015,49.2210847141386,55.2896847072645,61.7426792516797,68.5526366881019,75.6883415547682,83.1159359298789,90.7999757079674,98.70436299445,106.793097489269,115.030854429546,123.383379319793,131.817746549569,140.302476854765,148.807565911519,157.304448571381,165.765931249986,174.166083064684,182.480129029785,190.684355256758,198.756007061152,206.673215274601,214.414927512243,221.96086329307,229.291486522242,236.387975130425,243.232231968568,249.80687686919,256.09526300482,262.081498800008,267.750475165406,273.087890850223,278.080282749366,282.715066329811,286.980561296233,290.866028924552,294.361699349385,297.458810956704,300.149626676449,302.427467556521,304.286730241992,536.431832936876,751.177441888436,955.365750060786,1152.2808028855,1341.1804093601,1519.90154057349,1686.58626257878,1840.11634774943,1980.00873226163,2106.21817990318,2218.97601548017,2318.67908921722,2405.81776303892,2480.93058498128,2544.57620063128,2597.31628939643,2639.70543550747,2672.28535074762,2695.58182449233,2710.10338518468,2716.34102679394,2714.76856927755,2705.84342146309,2690.00756071496,2667.68860174271,2639.30090054034,2605.24664775496,2565.91687685858,2521.69243022394,2472.94483809627,2420.03707284887,2363.32428398479,2303.15432876616,2239.8683188472,2173.80102256161,2105.28120030327,2034.63186036811,1962.17044335603,1888.20894032843,1813.0539505115,1737.00668436291,1660.36291780757,1583.41290342138,1506.44124429877,1429.72673631915,1353.54218449747,1278.15419912759,1203.82297742381,1130.80207643222,1059.33818296517,989.670886214962,922.032458289431,856.647646820967,793.733481116684,733.499087122366,676.145492532286,621.865373255395,570.842625910249,523.251537628624,479.255258325911,438.979845250971,402.327973564152,369.018109172077,338.752442598622,311.253562483666,286.268307168292,263.56650642208,242.939030230084,224.195896529601,207.164518788144,191.688099874069,177.624202893236,164.843443046594,153.228274001279,142.671918677829,133.077395723445,124.356589076513,116.429478028152,109.223359648572,102.672218036656,96.7160712721448,91.3004624721247,86.3759126104703,81.8975044301678,77.8244289898049,74.1196414396585,70.7494888719584,67.6834217972237,64.8936916553918,62.3551070185338,60.0447892777967,57.9419654904442,56.0277720488413,54.2850784210725,52.6983296856962,51.2533968523437,49.9374505264446,48.7388346032204,47.646963580732,46.6522186260868,45.7458597041134,44.9199488866357,44.1672653909468,43.4812506614521,42.8559400038718,42.2859099931714,41.7662296615629,41.2924121083905,40.8603770081757,40.4664110318973,40.1071353627688,39.7794739917963,39.4806260073334,39.2080441033523,38.9594049651866,38.732594133751,38.5256869653966,38.3369281956453,38.1647194546882,38.0076035166529]
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
