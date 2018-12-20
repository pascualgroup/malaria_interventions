# Initialize --------------------------------------------------------------
source('/home/shai/Documents/malaria_interventions/functions.R')
prep.packages(c("stringr","ggplot2","igraph","Matrix","dplyr","cowplot"))

setwd('/media/Data/PLOS_Biol/empirical')

plotSurveyLayer <- function(x, zeroDiag=T, cutoff_g=NULL){
  if(zeroDiag){diag(x) <- 0}
  g <- graph.adjacency(x, mode = 'directed', weighted = T, diag = F)
  if(!is.null(cutoff_g)){g <- delete_edges(g, which(E(g)$weight<quantile(E(g)$weight, cutoff_g)))} # remove all edges smaller than the cutoff
  plot(g, 
       vertex.color='#8F25DD',
       vertex.size=6,
       vertex.label=NA,
       # edge.arrow.mode='-', 
       edge.arrow.width=1,
       edge.arrow.size=0.2,
       edge.curved=0.5, 
       edge.width=E(g)$weight*10)  
}



# Get and clean empirical data ---------------------------------------------------

mainData <- data.table::fread('S1ToS6_all_renamed_otuTable_binary.txt',sep = '\t',header = T)
names(mainData)[1] <- 'OTU_ID'
anyDuplicated(mainData$OTU_ID)
anyDuplicated(colnames(mainData))
varGroupings <- data.table::fread('S1ToS6_all_renamed_centroids_DBLaUpsType.csv', header = T)
varGroupings <- as.tibble(varGroupings)
varGroupings %<>% separate(col=type, into=c('var_type','sample','size'),sep=';')

setequal(mainData$OTU_ID,varGroupings$var_type)
setdiff(mainData$OTU_ID,varGroupings$var_type)
setdiff(varGroupings$var_type,mainData$OTU_ID)

mainData <- subset(mainData, OTU_ID%in%subset(varGroupings, Grouping=='BC')$var_type)
OTU_ID <- mainData$OTU_ID
mainData <- mainData[,-1]
mainData <- as.matrix(mainData)
rownames(mainData) <- OTU_ID
mainData[1:5,1:5]


# Get surveys for isolates
isolates_df <- tibble(isolate_code=colnames(mainData), isolate_name=str_sub(isolates_all,3,9))
isolates_df %<>% mutate(survey=str_sub(isolate_code, 1, str_locate(isolates_all, 'MRS')[,1]-1))
isolates_df %>% count(survey)
isolates_df %<>% 
  filter(isolate_name!='MRS1011') %>% # Remove the single individual with a chronic infection
  mutate(survey=ifelse(survey=='R2S1','S1',survey)) %>% 
  mutate(survey=ifelse(survey=='RS1','S1',survey)) %>% 
  mutate(survey=ifelse(survey=='RS4','S4',survey))
nrow(isolates_df)
isolates_df %<>% distinct()
nrow(isolates_df)

cleanSurveyData <- function(main_data, isolates_df, Survey, min.var.per.isolate=40, max.var.per.isolate=60, plotit=F){
  
  x <- isolates_df %>% filter(survey==Survey) %>% select(isolate_code)
  x <- x$isolate_code
  
  data_survey <- main_data[,x]
  data_survey <- bipartite::empty(data_survey)
  
  if(plotit){
    require(ggplot2)
    require(cowplot)
    x <- tibble(vars_in_pop=rowSums(data_survey))
    p1 <- ggplot(data = x, aes(x=vars_in_pop))+geom_histogram(fill='#1E9B95')+
      labs(x='var types in population', title='Number of times a var type\noccurs in the host population')+
      manuscript_theme
    x <- tibble(vars_in_isolate=colSums(data_survey))
    p2 <- ggplot(data = x, aes(x=vars_in_isolate))+geom_histogram(fill='#2C69C2')+
      labs(x='var types in isolate', title='Number of var types\nin an isolate')+
      geom_vline(xintercept = c(min.var.per.isolate,max.var.per.isolate), color='red')+
      manuscript_theme
    title <- ggdraw() + draw_label(Survey, size=20)
    p <- plot_grid(p1, p2,
                   labels = c("A", "B"),
                   nrow = 1, ncol = 2)
    p <- plot_grid(title, p, nrow=2, rel_heights = c(0.096,0.904))
    print(p)
  }
  data_survey <- data_survey[, which(colSums(data_survey)>=min.var.per.isolate & colSums(data_survey)<=max.var.per.isolate)]  # Remove isolates by constrains on number of var genes
  data_survey <- bipartite::empty(data_survey)
  dim(data_survey)
  # df <- df[which(rowSums(df)>1),] # Remove var types that appear only once
  # df <- bipartite::empty(df)
  cat(nrow(data_survey),' var types in ',ncol(data_survey),' isolates')
  return(data_survey)
}

# Limit to 60 vars (MOI=1)
maxVarIsolate <- 55
Data_S1 <- cleanSurveyData(main_data = mainData, isolates_df, Survey = 'S1', min.var.per.isolate = 40, max.var.per.isolate = maxVarIsolate, plotit = F)
Data_S2 <- cleanSurveyData(main_data = mainData, isolates_df, Survey = 'S2', min.var.per.isolate = 40, max.var.per.isolate = maxVarIsolate, plotit = F)
Data_S3 <- cleanSurveyData(main_data = mainData, isolates_df, Survey = 'S3', min.var.per.isolate = 40, max.var.per.isolate = maxVarIsolate, plotit = F)
Data_S4 <- cleanSurveyData(main_data = mainData, isolates_df, Survey = 'S4', min.var.per.isolate = 40, max.var.per.isolate = maxVarIsolate, plotit = F)
Data_S5 <- cleanSurveyData(main_data = mainData, isolates_df, Survey = 'S5', min.var.per.isolate = 40, max.var.per.isolate = maxVarIsolate, plotit = F)
Data_S6 <- cleanSurveyData(main_data = mainData, isolates_df, Survey = 'S6', min.var.per.isolate = 40, max.var.per.isolate = maxVarIsolate, plotit = F)

# isoaltes in rows, var genes in columns
Data_S1 <- t(Data_S1)
Data_S2 <- t(Data_S2)
Data_S3 <- t(Data_S3)
Data_S4 <- t(Data_S4)
Data_S5 <- t(Data_S5)
Data_S6 <- t(Data_S6)

# isolates_to_include <- unique(c(rownames(Data_S1),rownames(Data_S2),rownames(Data_S3),rownames(Data_S4),rownames(Data_S5),rownames(Data_S6)))
# vars_to_include <- unique(c(colnames(Data_S1),colnames(Data_S2),colnames(Data_S3),colnames(Data_S4),colnames(Data_S5),colnames(Data_S6)))
# mainDataClean <- mainData[vars_to_include,isolates_to_include]
# dim(mainDataClean)


# Define Layers -----------------------------------------------------------

empiricalLayer_1 <- overlapAlleleAdj(Data_S1)
empiricalLayer_2 <- overlapAlleleAdj(Data_S2)
empiricalLayer_3 <- overlapAlleleAdj(Data_S3)
empiricalLayer_4 <- overlapAlleleAdj(Data_S4)
empiricalLayer_5 <- overlapAlleleAdj(Data_S5)
empiricalLayer_6 <- overlapAlleleAdj(Data_S6)

diag(empiricalLayer_1) <- 0
diag(empiricalLayer_2) <- 0
diag(empiricalLayer_3) <- 0
diag(empiricalLayer_4) <- 0
diag(empiricalLayer_5) <- 0
diag(empiricalLayer_6) <- 0

intralayer_matrices_empirical <- list(empiricalLayer_1, empiricalLayer_2, empiricalLayer_3, empiricalLayer_4,empiricalLayer_5,empiricalLayer_6)

plotSurveyLayer(intralayer_matrices_empirical[[3]])

# Interlayer edges --------------------------------------------------------

# Build "interlayer networks"
print('Building inter-layer networks...')
interlayer_matrices_empirical <- list()
for (current_layer in 1:5){
  next_layer <- current_layer+1

  strain_copies_t <- rownames(intralayer_matrices_empirical[[current_layer]]) # repertoires at time t
  strain_copies_t1 <- rownames(intralayer_matrices_empirical[[next_layer]]) # repertoires at time t+1
  # need minimum of 2 strains in t and t+1 to build a matrix
  if (length(strain_copies_t)<2 | length(strain_copies_t1)<2){
    print(paste('No interlayer edges between layers ',current_layer,' and ',next_layer,' because there are < 2 repertoires in one of them.',sep=''))
    return(NULL)
  } 
  
  data_tmp <- bipartite::empty(mainData[,union(strain_copies_t,strain_copies_t1)])
  x <- overlapAlleleAdj(t(data_tmp)) # This is the similarity matrix for all the repertoires in both layers.
  inter_layer_edges_matrix <- x[strain_copies_t,strain_copies_t1]# Pull only the similarity values between the repertoires from the correct layers (i.e. create a bipartite)
  interlayer_matrices_empirical[[current_layer]] <- inter_layer_edges_matrix
  print(paste('Built interlayer edges for layers: ',current_layer,' (',length(strain_copies_t),' isolates with MOI=1)',' --> ',next_layer,' (',length(strain_copies_t1),' isolates with MOI=1)',sep=''))
}



# Get edge weight distributions from simulations --------------------------
edge_weights_simulated <- c()
for (ps in 500:599){
  print(ps)
  x <- get_edge_disributions(PS = ps, scenario = 'S', exp = '001', run = 1,cutoff_prob = 0.85, get_inter = F)
  x <- subset(x, value!=0)
  edge_weights_simulated <- rbind(edge_weights_simulated, x)
}

edge_weights_simulated <- as.tibble(edge_weights_simulated)
edge_weights_simulated$grp <- 'Simulated'


# Define cutoff -----------------------------------------------------------

# # Limit the surveys
# intralayer_matrices_empirical <- intralayer_matrices_empirical[1:4]
# interlayer_matrices_empirical <- interlayer_matrices_empirical[1:4]

cutoff_prob_empirical <- 0.85

# Get empirical edge weights
intralayer_edges <- unlist(sapply(intralayer_matrices_empirical, as.vector))
interlayer_edges <- unlist(sapply(interlayer_matrices_empirical, as.vector))
edges <- c(intralayer_edges,interlayer_edges)

#Create a data frame for empirical edge weights
edges_empirical <- tibble(value=intralayer_edges,grp='Empirical')
edges_empirical %<>% filter(value!=0) # remove instances with no edges

# Calculate cutoff based on empirical data
cutoff_value <- quantile(edges_empirical$value, probs = cutoff_prob_empirical)

# Plot distributions
edge_weights_simulated %>% 
  bind_rows(edges_empirical) %>%
  ggplot(aes(x=value,fill=grp))+
  geom_density(aes(y=..scaled..))+
  scale_fill_manual(values=c('#8F25DD','#10A4EF'))+
  geom_vline(xintercept = cutoff_value, color='red')

# Apply cutoff to layers
print('Applying cutoff to layers...')
for (i in 1:length(intralayer_matrices_empirical)){
  # print(i)
  x <- intralayer_matrices_empirical[[i]]
  x[x<cutoff_value] <- 0
  intralayer_matrices_empirical[[i]] <- x
}


# Distribution of repertoire persistence ----------------------------------
# In the neutral scenario, where repertoire persistence is purely a result of
# population dynamics, find the probability of persisting EXACTLY t layers. This
# is different than calculating the probablity of persisting at least t layers.
# to calculate the later, we need to look at some cumulative sum of the
# probablities: cumulative=1-cumsum(prob)

get_repertoire_persistence_without_modules <- function(PS,scenario,exp,run,cutoff_prob,folder='/media/Data/PLOS_Biol/Results/'){
  file <- paste(folder,PS,'_',scenario,'/','PS',PS,'_',scenario,'_E',exp,'_R',run,'_',cutoff_prob,'_repertoire_persistence_without_modules.txt',sep='')
  if(file.exists(file)){
    print(paste(PS,scenario,exp,run,cutoff_prob,sep=' | '))
    x <- read_delim(file, delim=',', col_types = 'ciii')
    x$PS <- PS
    x$scenario <- scenario
    x$exp <- exp
    x$run <- run
    return(x)
  } else {
    print(paste('File does not exist: ',file,sep=''))
    return()
  }
}

repertoire_persistence <- c()
for (ps in 500:599){
  x <- get_repertoire_persistence_without_modules(ps,'N','002',1,0.85)
  repertoire_persistence <- rbind(x,repertoire_persistence)
}

repertoire_persistence_prob <- repertoire_persistence %>%
  count(persistence) %>%
  mutate(prob=n/sum(n))
# Calculate the prob of persisting for at least a given amount of layers
repertoire_persistence_prob$cum_prob <- 1-(c(0,cumsum(repertoire_persistence_prob$prob)[-nrow(repertoire_persistence_prob)]))
repertoire_persistence_prob %>% print(n=Inf)

ggplot(repertoire_persistence_prob)+
  geom_line(aes(x=persistence, y=cum_prob),size=1)+ 
  scale_y_continuous(limits = c(0,1))+
  scale_x_continuous(limits = c(0,15), breaks = 0:15)


# Build Infomap objects -----------------------------------------------
network_object <- vector(mode = 'list', length = 2)
names(network_object) <- c('intralayer_matrices','interlayer_matrices')
network_object$intralayer_matrices <- intralayer_matrices_empirical
network_object$interlayer_matrices <- interlayer_matrices_empirical
infomap_empirical <- build_infomap_objects(network_object = network_object, write_to_infomap_file = T,
                                           infomap_file_name = '/media/Data/PLOS_Biol/empirical/infomap_empirical.txt', 
                                          return_objects = T,repertoire_persistence_prob = repertoire_persistence_prob)

infomap_empirical$infomap_interlayer %>% ggplot()+
  geom_density(aes(x=w),fill='purple',alpha=0.6)+
  geom_density(aes(x=w_rescaled),fill='gray',alpha=0.6)

# Run Infomap -------------------------------------------------------------

system("./Infomap_v01926 infomap_empirical.txt . -i multilayer -d -N 50 --rawdir --two-level --tree --expanded")


# Read Infomap results ----------------------------------------------------
node_list <- infomap_empirical$nodeList
print ('Reading infomap file...')
infomap_file <- '/media/Data/PLOS_Biol/empirical/infomap_empirical_expanded.tree'
lines <- readLines(infomap_file)
cat(lines[1]);cat('\n')
# x <- fread(infomap_file, skip = 2, stringsAsFactors = F) # Read results of Infomap
x <- read_delim(infomap_file, 
                delim = ' ',
                col_types = list(col_character(), col_double(), col_character(), col_integer(), col_integer()), 
                col_names = c('path', 'flow', 'name', 'layer', 'node'),
                skip = 2) # Read results of Infomap
print(x)

# Create a data frame to store results
modules <- tibble(module=rep(NA,nrow(x)),
                  nodeID=rep(NA,nrow(x)),
                  layer=rep(NA,nrow(x)),
                  path=x$path)

print('Creating module data frame...')
modules$module <- as.numeric(str_split(string = modules$path, pattern = ':', simplify = T)[,1])
modules$nodeID <- str_trim(str_split(string = x$name, pattern = '\\|', simplify = T)[,1])
modules$layer <- as.numeric(str_trim(str_split(string = x$name, pattern = '\\|', simplify = T)[,2])) # can also use x$layer

# modules %>%  ggplot(aes(x=layer,y=module))+geom_point()

# Rename modules because Infomap gives random names
print('Adding information on strains...')
modules2 <- modules %>% 
  distinct(module,layer) %>% 
  arrange(module,layer)
x <- c(1,table(modules2$module))
module_birth_layers <- modules2 %>% slice(cumsum(x)) %>% arrange(layer,module)
module_renaming <- data.frame(module=module_birth_layers$module, module_renamed = 1:max(module_birth_layers$module)) 
modules2 %<>% left_join(module_renaming)
modules2 %<>% full_join(modules) 
modules2 %<>% select(-module, -path)
names(modules2)[2] <- 'module'
modules2 %<>% arrange(module, layer, nodeID)

# Change node IDs to repertoire names
modules2$nodeID <- as.integer(modules2$nodeID)
node_list$nodeID <- as.integer(node_list$nodeID)
print(paste('Same strains in module and the node_list strain data frames?',setequal(modules2$nodeID,node_list$nodeID)))
modules2 %<>% left_join(node_list) %>% 
  rename(strain_unique=nodeLabel) %>% 
  mutate(strain_id=str_split(strain_unique,'_', simplify = T)[,1])

# Add information
modules2$PS <- 'Empirical'
modules2$scenario <- 'S'
modules2$exp <- '002'
modules2$run <- 1
modules2$cutoff_prob <- cutoff_prob_empirical

modules_empirical <- modules2


# Module analysis ---------------------------------------------------------

modules_empirical %>% ggplot(aes(x=layer,y=module))+geom_point()

modules_empirical %>% 
  group_by(module) %>% mutate(layer_birth=min(layer),layer_death=max(layer),persistence=layer_death-layer_birth+1) %>% 
  distinct(module,persistence)

modules_empirical %>% group_by(module) %>% summarise(numStrains=length(strain_unique))



# Proportiion of new moduels appearing after intervention

# Module/Repertoire distribution across isolates
# Because we work with MOI=1 then the distribution of repertoires in isoaltes 
# will always produce maximum entropy, since an isolate is a repertoire. So it 
# remains to see how modules are distributed. The number of cases is actually
# the number of strains (or isolates); again, becaues of MOI=1.


# Empirical ---------------------------------------------------------------

IRS_S <- get_data(parameter_space = '18', scenario = 'S', experiment = '002', run = 1, cutoff_prob = 0.85, use_sqlite = T, tables_to_get = 'summary_general')[[1]]
IRS_G <- get_data(parameter_space = '18', scenario = 'G', experiment = '002', run = 1, cutoff_prob = 0.85, use_sqlite = T, tables_to_get = 'summary_general')[[1]]
IRS_N <- get_data(parameter_space = '18', scenario = 'N', experiment = '002', run = 1, cutoff_prob = 0.85, use_sqlite = T, tables_to_get = 'summary_general')[[1]]

IRS_S$layer <- 1:300
IRS_G$layer <- 1:300
IRS_N$layer <- 1:300
IRS_S %>% bind_rows(IRS_G) %>% bind_rows(IRS_N) %>% 
  mutate(scenario=factor(scenario, levels=c('S','G','N'))) %>% 
  filter(layer%in%100:200) %>% 
  select(-year, -month, -n_infected) %>% 
  gather(variable, value, -pop_id, -time, -exp, -PS, -scenario, -run, -layer) %>% 
  group_by(pop_id, time, layer, exp, PS, scenario, variable) %>%
  summarise(value_mean=mean(value), value_sd=sd(value)) %>% # Need to average across runs
  filter(variable %in% monitored_variables) %>%
  ggplot()+
  geom_line(aes(x=layer, y=value_mean, color=scenario))+
  geom_errorbar(aes(ymin=value_mean-value_sd,ymax=value_mean+value_sd,x=layer, group=scenario),color='gray',width=0.001, alpha=0.3)+
  geom_vline(xintercept = c(131,149),color-'black')+
  scale_color_manual(values = scenario_cols)+
  facet_wrap(~variable, scales='free')+
  mytheme




# Analyzing simulations ---------------------------------------------------
