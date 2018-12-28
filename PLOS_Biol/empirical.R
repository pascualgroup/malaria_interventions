# Initialize --------------------------------------------------------------
source('/home/shai/Documents/malaria_interventions/functions.R')
prep.packages(c("tidyverse","magrittr","data.table","igraph","Matrix","dplyr","cowplot",'sqldf'))

setwd('/media/Data/PLOS_Biol/empirical')

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
dim(mainData)
mainData[1:5,1:5]


# Diversity of var genes##################
dim(mainData)
var_numbers <- rowSums(mainData)
hist(var_numbers)
vegan::diversity(var_numbers)
vegan::diversity(rep(1,9000))
##########################################

isolates_all <- colnames(mainData)

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
# Run once, and then just read results from file:
# unlink('/media/Data/PLOS_Biol/empirical/edge_weights_simulated.csv')
# edge_weights_simulated <- c()
# for (ps in c(520,540,560,580,599)){
#   print(paste('[',Sys.time(),']', ps))
#   x <- get_edge_disributions(PS = ps, scenario = 'S', exp = '001', run = 1,cutoff_prob = 0.95, get_inter = F) # The cutoff_prob here is just for file names. it does not really matter because these ar eedge weights before cutoff is applied
#   # x <- subset(x, value!=0)
#   x <- x[sample(nrow(x),round(nrow(x)*0.4),F),] # Randomly sample 20% of the edge values because there are just too many.
#   edge_weights_simulated <- rbind(edge_weights_simulated, x)
#   if (ps%%20==0){
#     fwrite(edge_weights_simulated[,-c(6,7)], '/media/Data/PLOS_Biol/empirical/edge_weights_simulated.csv', append = T)
#     edge_weights_simulated <- c()
#   }
# }
edge_weights_simulated <- fread('/media/Data/PLOS_Biol/empirical/edge_weights_simulated.csv')
edge_weights_simulated <- as.tibble(edge_weights_simulated)
edge_weights_simulated$grp <- 'Simulated'

# See where the quantile passes 0.
quantile(edge_weights_simulated$value, probs = seq(0.85,1,0.005))

ggplot(edge_weights_simulated, aes(x=value))+
  scale_x_continuous(breaks=seq(0,1,0.05))+
  geom_density()


# Explore simulated data ------------------------------------------

monitored_variables <- c('prevalence', 'meanMOI','n_circulating_strains', 'n_circulating_genes', 'n_alleles', 'n_total_bites')

PS_range <- as.character(500:599)
# cases <- expand.grid(ps=PS_range, scenario='S', exp=c('001','002'), run=1)
# cases$cutoff_prob <- 0.85 # This is just for fime names. there is no cutoff because there are no modules here.
# exploratory <- c()
# for (i in 1:nrow(cases)){
#   print(paste('PS: ',cases$ps[i],' | Scenario: ',cases$scenario[i],' | exp: ',cases$exp[i], ' | run: ',cases$run[i],sep=''))
#   tmp <- get_data(parameter_space = cases$ps[i], scenario = cases$scenario[i], experiment = cases$exp[i], run = cases$run[i], cutoff_prob = cases$cutoff_prob[i], use_sqlite = T, tables_to_get = 'summary_general')[[1]]
#   exploratory <- rbind(exploratory, tmp)
# }
# write_csv(exploratory, '/media/Data/PLOS_Biol/Results/exploratory.csv')
exploratory <- read_csv('/media/Data/PLOS_Biol/Results/exploratory.csv')
exploratory$month <- factor(exploratory$month, levels=c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'))
# EIR
exploratory %>% 
  ggplot(aes(x=month, y=EIR))+
  geom_boxplot()+
  stat_summary(aes(group=exp), fun.y=mean, geom="line", size=1)+
  facet_wrap(~exp)+
  scale_y_continuous(breaks=seq(0,20,2))+
  mytheme


layers_df <- exploratory %>% distinct(time,year,month) %>% mutate(layer=1:300) %>% mutate(mylabels=paste(layer,year,month,sep='_'))

exploratory_summary <- exploratory %>%
  left_join(layers_df) %>% 
  select(-n_infected) %>% 
  # filter(time>time_range[1]&time<time_range[2]) %>%
  gather(variable, value, -pop_id, -year, -month, -time, -layer, -exp, -PS, -scenario, -run, -mylabels) %>% 
  group_by(time, layer, year, month, mylabels, exp, scenario, variable) %>%
  summarise(value_mean=mean(value), value_sd=sd(value)) %>% # Need to average across PS
  filter(variable %in% monitored_variables)
# First make sure that the location of the surveys and IRS in the simulations is correct
exploratory_summary %>% 
  filter(layer %in% 110:180) %>%
  filter(variable %in% c('n_total_bites')) %>% 
  ggplot()+
  geom_line(aes(x=layer, y=value_mean, color=exp))+
  geom_errorbar(aes(x=layer, ymin=value_mean-value_sd,ymax=value_mean+value_sd, color=exp),width=0.2, alpha=0.3)+
  scale_color_manual(values = c('gray50','red'))+
  scale_x_continuous(breaks=110:180, labels = unique(subset(exploratory_summary, layer%in%110:180)$mylabels))+
  geom_vline(xintercept = c(118,126,138,142,154,162), color='black', size=1)+
  geom_vline(xintercept = c(131,138,145), color='black', linetype='dotted', size=1)+ # IRS
  facet_wrap(~variable, scales='free')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size=8))
# Now plot the variables
exploratory_summary %>% 
  filter(layer %in% 110:180) %>%
  filter(variable!='n_total_bites') %>% 
  ggplot()+
  geom_line(aes(x=layer, y=value_mean, color=exp))+
  geom_errorbar(aes(x=layer, ymin=value_mean-value_sd,ymax=value_mean+value_sd, color=exp),width=0.2, alpha=0.3)+
  scale_color_manual(values = c('gray50','red'))+
  geom_vline(xintercept = c(131,138,145), color='black', linetype='dotted', size=0.5)+ # IRS
  geom_vline(xintercept = c(118,126,138,142,154,162), color='black', size=0.5)+
  facet_wrap(~variable, scales='free')+
  mytheme
 
exploratory_summary %>% 
  filter(layer %in% c(118,126,138,142,154,162)) %>% 
  filter(exp=='002') %>% 
  filter(variable=='n_circulating_genes')

# Number of genes per survey in the simulated vs. empirical data. This is only
# for the subsampled repertoires Do that by going directly to the node list, to
# ensure that the genes I count belong to strains that were fed into Infomap.
num_genes <- c()
for (ps in 500:599){
  print(ps)
  sampled_strains <- read_csv(paste('/media/Data/PLOS_Biol/Results/',ps,'_S/PS',ps,'_S_E002_R1_0.85_sampled_strains.csv',sep=''))
  sampled_strains$strain_id <- as.character(sampled_strains$strain_id)
  sampled_alleles <- read_csv(paste('/media/Data/PLOS_Biol/Results/',ps,'_S/PS',ps,'_S_E002_R1_0.85_sampled_alleles.csv',sep=''))
  node_list <- read_csv(paste('/media/Data/PLOS_Biol/Results/',ps,'_S/PS',ps,'_S_E002_R1_0.85_node_list.csv',sep=''))
  node_list %<>% mutate(strain_id=str_split(nodeLabel,'_')[[1]][1])
  x <- node_list %>% 
    left_join(sampled_strains) %>% 
    left_join(sampled_alleles) %>%
    distinct(strain_id,gene_id)
  num_genes <- c(num_genes, length(unique(x$gene_id)))
}

empirical_data_genes <- c(colnames(Data_S1),
                          colnames(Data_S2),
                          colnames(Data_S3),
                          colnames(Data_S4),
                          colnames(Data_S5),
                          colnames(Data_S6))
empirical_data_genes <- unique(empirical_data_genes)                    

mean(num_genes)
length(empirical_data_genes)

# Cutoff sensitivity in simulated data ------------------------------------

 
# Create files to run tests on Midway
sbatch_arguments <- expand.grid(PS=sprintf('%0.3d', 500:599),
                              scen=c('S','G','N'),
                              array='1', 
                              layers='118,126,138,142,154,162',
                              exp=c('002'),
                              cutoff_prob=0.85,
                              time = '01:00:00',
                              mem_per_cpu = 4000,
                              stringsAsFactors = F)
nrow(sbatch_arguments)
make_sbatch_get_data(sbatch_arguments = sbatch_arguments,
                     make_networks = F,
                     prepare_infomap = T, # make network is nested in this, so no need to call it
                     run_Infomap = T,
                     read_infomap_results = T,
                     temporal_diversity = F,
                     module_Fst = F)



experiments <- expand.grid(PS=sprintf('%0.3d', 500:599),
                           scen=c('S','G','N'),
                           exp=c('001','002'),
                           # cutoff_prob=seq(0.95,0.97,0.005),
                           cutoff_prob=0.85,
                           stringsAsFactors = F)

# factorial_design_test <- tibble(option=1:4,unit=c('genes','alleles','alleles','genes'),MOI=c('1','1','All','All'),cutoff_prob=c(0.95,0.85,0.85,0.95))

module_results_simulations <- c()
# interlayer_edge_weights <- c()
for (i in 1:nrow(experiments)){
  print(i)
  ps <- experiments$PS[i]
  scenario <- experiments$scen[i]
  cutoff_prob <- experiments$cutoff_prob[i]
  exp <- experiments$exp[i]
  x <- get_modularity_results(ps,scenario,exp,1,cutoff_prob,folder = paste('/media/Data/PLOS_Biol/Results/',ps,'_',scenario,'_3','/',sep=''))
  # y <- readLines(paste('/media/Data/PLOS_Biol/Results/',ps,'_',scenario,'_3','/PS',ps,'_',scenario,'_E',exp,'_R1_',cutoff_prob,'_network_info.csv',sep=''))[6]
  # x$cutoff_value <- as.numeric(y)
  module_results_simulations <- rbind(module_results_simulations,x)
  
  # x <- read_delim(paste('/media/Data/PLOS_Biol/Results/',ps,'_S/PS',ps,'_S_E',exp,'_R1_',cutoff_prob,'_Infomap_multilayer.txt',sep=''), delim=' ', col_names = c('layer_s','node_s','layer_t','node_t','w'))
  # x$cutoff_prob <- cutoff_prob
  # x$exp <- exp
  # interlayer_edge_weights <- rbind(interlayer_edge_weights,x)
}
module_results_simulations <- as.tibble(module_results_simulations)
module_results_simulations %>% filter(exp=='002') %>%  distinct(scenario, PS, exp)

ggplot(interlayer_edge_weights, aes(w, fill=exp))+geom_density()+facet_wrap(~cutoff_prob, scales='free')
interlayer_edge_weights %>% group_by(cutoff_prob,exp) %>% summarise(min_w=min(w),max_w=max(w))


modsize <- module_results_simulations %>% 
  group_by(scenario,PS,exp,cutoff_prob,module) %>% summarise(size=n())
module_results_simulations %>% 
  # filter(exp=='002') %>% 
  filter(PS==599) %>%
  # filter(scenario=='S') %>% 
  distinct(PS,scenario,cutoff_prob,module, layer) %>% 
  group_by(PS,scenario,cutoff_prob,module) %>% 
  summarise(birth_layer=min(layer),death_layer=max(layer)+1) %>% 
  left_join(modsize) %>% 
  ggplot(aes(xmin=birth_layer, xmax=death_layer, ymin=module, ymax=module,color=size))+
  geom_rect(size=2)+
  # geom_rect(size=2, color=scenario_cols[1])+
  scale_x_continuous(breaks=1:6)+
  scale_color_viridis_c()+
  facet_grid(exp~scenario)


modsize <- module_results_simulations %>% 
  group_by(scenario,PS,exp,cutoff_prob,layer,module) %>% summarise(size=n())
modsize %>%
  filter(PS==520) %>%
  # filter(exp=='002') %>%
  ggplot(aes(x=layer,y=size))+
  geom_bar(aes(fill=as.factor(module)),stat = "identity",position='stack', color='black')+
  # geom_area(aes(color=as.factor(module), fill=as.factor(module)),position='stack')+
  facet_grid(exp~scenario)+
  mytheme+theme(legend.position='none')

# cutoffs
# cutoffs <- module_results_simulations %>% distinct(scenario, PS, exp, cutoff_prob, cutoff_value) %>% print(n=Inf)

# Persistence
module_persistence <- module_results_simulations %>% 
  select(scenario, exp, PS, run, cutoff_prob, layer, module) %>% 
  group_by(scenario, exp, PS,run,cutoff_prob,module) %>% 
  summarise(birth_layer=min(layer), death_layer=max(layer), persistence=death_layer-birth_layer+1)

total_modules <- module_persistence %>% 
  group_by(scenario, exp, PS, cutoff_prob) %>% summarise(n_modules=length(unique(module)))

module_persistence %>% 
  group_by(scenario, exp, PS, cutoff_prob) %>% 
  count(persistence) %>% left_join(total_modules) %>% mutate(prop=n/n_modules) %>% 
  # filter(cutoff_prob==0.95) %>% 
  ggplot(aes(y=prop, x=persistence,fill=scenario))+
  facet_wrap(~exp)+
  geom_col(position='dodge')+
  scale_fill_manual(values=c('blue','orange','red'))

# Proportion of modules persisting up to S3
module_persistence %>% 
  filter(death_layer==3) %>% 
  group_by(scenario, exp, PS, cutoff_prob) %>% summarise(n_1to3=n()) %>% 
  left_join(total_modules) %>% mutate(prop=n_1to3/n_modules) %>% 
  ggplot(aes(x=scenario,y=prop,fill=scenario))+
  geom_boxplot()+
  facet_wrap(~exp)+
  scale_fill_manual(values=c('blue','orange','red'))

# Proportion of modules born after layer 3
module_persistence %>% 
  filter(birth_layer==3) %>% 
  group_by(scenario, exp, PS, cutoff_prob) %>% summarise(n_after_3=n()) %>% 
  left_join(total_modules) %>% mutate(prop=n_after_3/n_modules) %>% 
  ggplot(aes(x=scenario,y=prop,fill=scenario))+
  geom_boxplot()+
  facet_wrap(~exp)+
  scale_fill_manual(values=c('blue','orange','red'))

# Proportion of modules that passed layer 3
module_persistence %>% 
  filter(birth_layer<3, death_layer>3) %>% 
  group_by(scenario, exp, PS, cutoff_prob) %>% summarise(n_after_3=n()) %>% 
  left_join(total_modules) %>% mutate(prop=n_after_3/n_modules) %>% 
  ggplot(aes(x=scenario,y=prop,fill=scenario))+
  geom_boxplot()+
  facet_wrap(~exp)+
  scale_fill_manual(values=c('blue','orange','red'))

# Survival analysis
library(survival)
library(survminer)
x <- module_persistence
x <- subset(x, x$birth_layer!=x$death_layer)
x <- subset(x, birth_layer==1)
unique(x$scenario)
unique(x$exp)
x$event <- ifelse(x$death_layer==6,0,1)

fit <- with(x, survfit(Surv(birth_layer, death_layer, event)~scenario+exp))
ggsurvplot_facet(fit, data=x, facet.by = 'exp', conf.int = TRUE)
# +scale_color_manual(values=c('blue','orange','red'))


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

# Copy that file to Midway so it can be used for the simulations
write_csv(repertoire_persistence_prob,'/media/Data/PLOS_Biol/empirical/repertoire_persistence_prob.csv')


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

system("./Infomap_v01926 infomap_empirical.txt . -i multilayer -d -N 10 --rawdir --two-level --tree --expanded")


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

# Temporal diversity
# 
# max_layer_to_persist <- max(modules_empirical$layer)-min(modules_empirical$layer)+1 # This is important for analyses that do not have sequential number of layers (1:300), like in the empirical IRS.
# 
# module_persistence <- modules_empirical %>% 
#   select(layer, module) %>% 
#   group_by(module) %>% 
#   summarise(birth_layer=min(layer), death_layer=max(layer), persistence=death_layer-birth_layer+1) %>% 
#   mutate(relative_persistence=persistence/(max_layer_to_persist-birth_layer+1))
# 
# sampled_strains <- read_csv(file_strains, col_types = 'ccc')
# sampled_strains <-  sampled_strains[,-3]
# suppressMessages(modules %<>% select(scenario, PS, scenario, exp, run, cutoff_prob, module, strain_id) %>% left_join(sampled_strains))
# allele_freq <- xtabs(~module+allele_locus, modules)
# module_diversity <- vegan::diversity(allele_freq)/log(ncol(allele_freq))
# 
# module_persistence$D <- module_diversity
# module_persistence$statistic <- module_diversity*module_persistence$relative_persistence




# Read results of simulated data ------------------------------------------

module_results_simulations <- temporal_diversity_simulations <- mFst_simulations<- c()
for (ps in 500:599){
  for(scen in c('S','G','N')){
    x <- get_modularity_results(ps,scen,'002',1,0.85,folder = paste('/media/Data/PLOS_Biol/Results/',ps,'_',scen,'/',sep=''))
    module_results_simulations <- rbind(module_results_simulations,x)
    
    x <- get_temporal_diversity(ps,scen,'002',1,0.85,folder = paste('/media/Data/PLOS_Biol/Results/',ps,'_',scen,'/',sep=''))
    temporal_diversity_simulations <- rbind(temporal_diversity_simulations,x)
    
    x <- get_mFst(ps,scen,'002',1,0.85,folder = paste('/media/Data/PLOS_Biol/Results/',ps,'_',scen,'/',sep=''))
    mFst_simulations <- rbind(mFst_simulations,x)
  }
}
module_results_simulations <- as.tibble(module_results_simulations)
temporal_diversity_simulations <- as.tibble(temporal_diversity_simulations)
mFst_simulations <- as.tibble(mFst_simulations)

module_results_simulations %>% 
  filter(PS==500) %>%
  filter(scenario=='N') %>% 
  distinct(module, layer) %>% 
  group_by(module) %>% summarise(birth_layer=min(layer),death_layer=max(layer)+1) %>% 
  ggplot(aes(xmin=birth_layer, xmax=death_layer, ymin=module, ymax=module))+
  geom_rect(size=2, color=scenario_cols[1])+
  scale_x_continuous(breaks=1:6)+
  manuscript_theme+theme(axis.title = element_blank())


# Relative persistence
module_persistence <- module_results_simulations %>% 
  select(scenario, PS, run, cutoff_prob, layer, module) %>% 
  group_by(scenario, PS,run,cutoff_prob,module) %>% 
  summarise(birth_layer=min(layer), death_layer=max(layer), persistence=death_layer-birth_layer+1) %>% 
  mutate(relative_persistence=persistence/(300-birth_layer+1)) %>% 
  mutate(type='Module') %>% 
  rename(id=module) %>% mutate(id=as.character(id))

ggplot(module_persistence)+
  geom_density(aes(persistence,fill=scenario, y=..scaled..))
  # geom_histogram(aes(persistence,fill=scenario), position='dodge')


temporal_diversity_simulations %>%
  ggplot(aes(statistic))+geom_density()







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
