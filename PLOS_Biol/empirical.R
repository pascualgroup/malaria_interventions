# Initialize --------------------------------------------------------------
## @knitr Initialize
source('~/Documents/malaria_interventions/functions.R')
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

## @knitr END

# Get and clean empirical data (GENES) ---------------------------------------------------
data_genes <- data.table::fread('~/Dropbox/Qixin_Shai_Malaria/S1-S6/S1ToS6_all_renamed_otuTable_binary.txt',sep = '\t',header = T)
names(data_genes)[1] <- 'OTU_ID'
anyDuplicated(data_genes$OTU_ID)
anyDuplicated(colnames(data_genes))
varGroupings <- data.table::fread('~/Dropbox/Qixin_Shai_Malaria/S1-S6/S1ToS6_all_renamed_centroids_DBLaUpsType.csv', header = T)
varGroupings <- as.tibble(varGroupings)
varGroupings %<>% separate(col=type, into=c('var_type','sample','size'),sep=';')

setequal(data_genes$OTU_ID,varGroupings$var_type)
setdiff(data_genes$OTU_ID,varGroupings$var_type)
setdiff(varGroupings$var_type,data_genes$OTU_ID)

data_genes <- subset(data_genes, OTU_ID%in%subset(varGroupings, Grouping=='BC')$var_type)
OTU_ID <- data_genes$OTU_ID
data_genes <- data_genes[,-1]
data_genes <- as.matrix(data_genes)
rownames(data_genes) <- OTU_ID
dim(data_genes)
data_genes[1:5,1:5]

isolates_all <- colnames(data_genes)

# Get surveys for isolates
isolates_df <- tibble(isolate_code=colnames(data_genes), isolate_name=str_sub(isolates_all,3,9))
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

# Diversity of var genes
var_gene_richness_empirical <- c()
var_gene_diversity_empirical <- c()
for (s in paste('S',1:6,sep='')){
  x <- data_genes[,subset(isolates_df, survey==s)$isolate_code]
  x <- bipartite::empty(x)
  var_gene_richness_empirical <- c(var_gene_richness_empirical,nrow(x))
  var_gene_diversity_empirical <- c(var_gene_diversity_empirical,vegan::diversity(rowSums(x)))
}


# Limit to 60 vars (MOI=1)
maxVarIsolate <- 55
Data_S1_genes <- cleanSurveyData(main_data = data_genes, isolates_df, Survey = 'S1', min.var.per.isolate = 40, max.var.per.isolate = maxVarIsolate, plotit = F)
Data_S2_genes <- cleanSurveyData(main_data = data_genes, isolates_df, Survey = 'S2', min.var.per.isolate = 40, max.var.per.isolate = maxVarIsolate, plotit = F)
Data_S3_genes <- cleanSurveyData(main_data = data_genes, isolates_df, Survey = 'S3', min.var.per.isolate = 40, max.var.per.isolate = maxVarIsolate, plotit = F)
Data_S4_genes <- cleanSurveyData(main_data = data_genes, isolates_df, Survey = 'S4', min.var.per.isolate = 40, max.var.per.isolate = maxVarIsolate, plotit = F)
Data_S5_genes <- cleanSurveyData(main_data = data_genes, isolates_df, Survey = 'S5', min.var.per.isolate = 40, max.var.per.isolate = maxVarIsolate, plotit = F)
Data_S6_genes <- cleanSurveyData(main_data = data_genes, isolates_df, Survey = 'S6', min.var.per.isolate = 40, max.var.per.isolate = maxVarIsolate, plotit = F)

data_moi1_genes <- list(Data_S1_genes,Data_S2_genes,Data_S3_genes,Data_S4_genes)
sapply(data_moi1_genes,dim)
# data_moi1_genes <- list(Data_S1,Data_S2,Data_S3,Data_S4,Data_S5,Data_S6)

# isolates_to_include <- unique(c(rownames(Data_S1),rownames(Data_S2),rownames(Data_S3),rownames(Data_S4),rownames(Data_S5),rownames(Data_S6)))
# vars_to_include <- unique(c(colnames(Data_S1),colnames(Data_S2),colnames(Data_S3),colnames(Data_S4),colnames(Data_S5),colnames(Data_S6)))
# data_genesClean <- data_genes[vars_to_include,isolates_to_include]
# dim(data_genesClean)

# Get and clean empirical data (ALLELES) ---------------------------------------------------

## @knitr prepare_data
data_allele_1 <- data.table::fread('~/Dropbox/Qixin_Shai_Malaria/S1-S6/S1ToS6_VR1_clusterMap_groupBC_table.txt',sep = '\t',header = T)
# anyDuplicated(data_allele_1$OTU_ID)
# anyDuplicated(colnames(data_allele_1))
OTU_ID <- data_allele_1$OTU_ID
data_allele_1 <- data_allele_1[,-1]
data_allele_1 <- as.matrix(data_allele_1)
rownames(data_allele_1) <- OTU_ID
# dim(data_allele_1)
# data_allele_1[1:5,1:5]
data_allele_2 <- data.table::fread('~/Dropbox/Qixin_Shai_Malaria/S1-S6/S1ToS6_VR2_clusterMap_groupBC_table.txt',sep = '\t',header = T)
# anyDuplicated(data_allele_2$OTU_ID)
# anyDuplicated(colnames(data_allele_2))
OTU_ID <- data_allele_2$OTU_ID
data_allele_2 <- data_allele_2[,-1]
data_allele_2 <- as.matrix(data_allele_2)
rownames(data_allele_2) <- OTU_ID
# dim(data_allele_2)
# data_allele_2[1:5,1:5]
setequal(colnames(data_allele_1),colnames(data_allele_2))

data_allleles <- rbind(data_allele_1,data_allele_2)
isolates_all_allleles <- colnames(data_allleles)

# Get surveys for isolates
isolates_df_alleles <- tibble(isolate_code=colnames(data_allleles), isolate_name=str_sub(isolates_all_allleles,3,9))
isolates_df_alleles %<>% mutate(survey=str_sub(isolate_code, 1, str_locate(isolates_all_allleles, 'MRS')[,1]-1))
isolates_df_alleles %>% count(survey)
isolates_df_alleles %<>% 
  filter(isolate_name!='MRS1011') %>% # Remove the single individual with a chronic infection
  mutate(survey=ifelse(survey=='R2S1','S1',survey)) %>% 
  mutate(survey=ifelse(survey=='RS1','S1',survey)) %>% 
  mutate(survey=ifelse(survey=='RS4','S4',survey))
nrow(isolates_df_alleles)
isolates_df_alleles %<>% distinct()
nrow(isolates_df_alleles)

# Limit to 60 vars (MOI=1)
maxVarIsolate <- 110
Data_S1_alleles <- cleanSurveyData(main_data = data_allleles, isolates_df_alleles, Survey = 'S1', min.var.per.isolate = 80, max.var.per.isolate = maxVarIsolate, plotit = F)
Data_S2_alleles <- cleanSurveyData(main_data = data_allleles, isolates_df_alleles, Survey = 'S2', min.var.per.isolate = 80, max.var.per.isolate = maxVarIsolate, plotit = F)
Data_S3_alleles <- cleanSurveyData(main_data = data_allleles, isolates_df_alleles, Survey = 'S3', min.var.per.isolate = 80, max.var.per.isolate = maxVarIsolate, plotit = F)
Data_S4_alleles <- cleanSurveyData(main_data = data_allleles, isolates_df_alleles, Survey = 'S4', min.var.per.isolate = 80, max.var.per.isolate = maxVarIsolate, plotit = F)
Data_S5_alleles <- cleanSurveyData(main_data = data_allleles, isolates_df_alleles, Survey = 'S5', min.var.per.isolate = 80, max.var.per.isolate = maxVarIsolate, plotit = F)
Data_S6_alleles <- cleanSurveyData(main_data = data_allleles, isolates_df_alleles, Survey = 'S6', min.var.per.isolate = 80, max.var.per.isolate = maxVarIsolate, plotit = F)

# isoaltes in rows, var genes in columns
Data_S1_alleles <- t(Data_S1_alleles)
Data_S2_alleles <- t(Data_S2_alleles)
Data_S3_alleles <- t(Data_S3_alleles)
Data_S4_alleles <- t(Data_S4_alleles)
Data_S5_alleles <- t(Data_S5_alleles)
Data_S6_alleles <- t(Data_S6_alleles)

data_moi1_alleles <- list(Data_S1_alleles,Data_S2_alleles,Data_S3_alleles,Data_S4_alleles,Data_S5_alleles,Data_S6_alleles)

# dimensions of data and freq. of alleles
x <- data_allleles[unique(unlist(lapply(data_moi1_alleles,colnames))),unique(unlist(lapply(data_moi1_alleles,rownames)))]
dim(x)
x <- tibble(freq=rowSums(x))
ggplot(x, aes(freq))+geom_density()

## @knitr END


# Compare gene to allele data
sapply(data_moi1_genes,dim)
sapply(data_moi1_alleles,dim)
isolates_genes <- lapply(data_moi1_genes,colnames)
isolates_alleles <- lapply(data_moi1_alleles,colnames)
for (i in 1:6){
  cat('# shared isolates in data sets, layer ');cat(i);cat(' ');
  cat(length(base::intersect(isolates_genes[[i]],isolates_alleles[[i]])))
  cat('\n')
}

# Define Layers -----------------------------------------------------------

## @ knitr build_network
empiricalLayer_1 <- overlapAlleleAdj(t(Data_S1_genes))
empiricalLayer_2 <- overlapAlleleAdj(t(Data_S2_genes))
empiricalLayer_3 <- overlapAlleleAdj(t(Data_S3_genes))
empiricalLayer_4 <- overlapAlleleAdj(t(Data_S4_genes))
empiricalLayer_5 <- overlapAlleleAdj(t(Data_S5_genes))
empiricalLayer_6 <- overlapAlleleAdj(t(Data_S6_genes))

diag(empiricalLayer_1) <- 0
diag(empiricalLayer_2) <- 0
diag(empiricalLayer_3) <- 0
diag(empiricalLayer_4) <- 0
diag(empiricalLayer_5) <- 0
diag(empiricalLayer_6) <- 0

intralayer_matrices_empirical <- list(empiricalLayer_1, empiricalLayer_2, empiricalLayer_3, empiricalLayer_4,empiricalLayer_5,empiricalLayer_6)
# intralayer_matrices_empirical <- list(empiricalLayer_1, empiricalLayer_2, empiricalLayer_3, empiricalLayer_4)



# plotSurveyLayer(intralayer_matrices_empirical[[3]])

# Interlayer edges --------------------------------------------------------

# Build "interlayer networks"
print('Building inter-layer networks...')
interlayer_matrices_empirical <- list()
data_for_interlayer <- data_genes
for (current_layer in 1:(length(intralayer_matrices_empirical)-1)){
  next_layer <- current_layer+1
  
  strain_copies_t <- rownames(intralayer_matrices_empirical[[current_layer]]) # repertoires at time t
  strain_copies_t1 <- rownames(intralayer_matrices_empirical[[next_layer]]) # repertoires at time t+1
  # need minimum of 2 strains in t and t+1 to build a matrix
  if (length(strain_copies_t)<2 | length(strain_copies_t1)<2){
    print(paste('No interlayer edges between layers ',current_layer,' and ',next_layer,' because there are < 2 repertoires in one of them.',sep=''))
    return(NULL)
  } 
  
  data_tmp <- bipartite::empty(data_for_interlayer[,union(strain_copies_t,strain_copies_t1)])
  x <- overlapAlleleAdj(t(data_tmp)) # This is the similarity matrix for all the repertoires in both layers.
  inter_layer_edges_matrix <- x[strain_copies_t,strain_copies_t1]# Pull only the similarity values between the repertoires from the correct layers (i.e. create a bipartite)
  interlayer_matrices_empirical[[current_layer]] <- inter_layer_edges_matrix
  print(paste('Built interlayer edges for layers: ',current_layer,' (',length(strain_copies_t),' isolates with MOI=1)',' --> ',next_layer,' (',length(strain_copies_t1),' isolates with MOI=1)',sep=''))
}

## @knitr END

# Define cutoff -----------------------------------------------------------

# # Limit the surveys
# intralayer_matrices_empirical <- intralayer_matrices_empirical[1:4]
# interlayer_matrices_empirical <- interlayer_matrices_empirical[1:4]

## @knitr set_cutoff

cutoff_prob_empirical <- 0.9

# Get empirical edge weights
intralayer_edges <- unlist(sapply(intralayer_matrices_empirical, as.vector))
quantile(intralayer_edges, probs = seq(0.85,1,0.005))
cutoff_value_intra <- quantile(intralayer_edges, probs = cutoff_prob_empirical)

interlayer_edges <- unlist(sapply(interlayer_matrices_empirical, as.vector))
quantile(interlayer_edges, probs = seq(0.85,1,0.005))
cutoff_value_inter <- quantile(interlayer_edges, probs = cutoff_prob_empirical)

# Plot distributions before applying cutoff

png('empirical_data_edge_weights.png')
intralayer_edges <- unlist(sapply(intralayer_matrices_empirical, as.vector))
interlayer_edges <- unlist(sapply(interlayer_matrices_empirical, as.vector))
tibble(w=intralayer_edges,grp='intra') %>%
  bind_rows(tibble(w=interlayer_edges,grp='inter')) %>% 
  ggplot(aes(x=w,fill=grp))+
  geom_density()+
  scale_fill_manual(values=c('#8F25DD','#10A4EF'))+
  geom_vline(xintercept = c(cutoff_value_intra,cutoff_value_inter), color=c('#10A4EF','#8F25DD'))
dev.off()

# Apply cutoff to layers
print('Applying cutoff to layers...')
for (i in 1:length(intralayer_matrices_empirical)){
  # print(i)
  x <- intralayer_matrices_empirical[[i]]
  x[x<cutoff_value_intra] <- 0
  intralayer_matrices_empirical[[i]] <- x
}
# Apply cutoff to inter-layer blocks
print('Applying cutoff to layers...')
for (i in 1:length(interlayer_matrices_empirical)){
  # print(i)
  x <- interlayer_matrices_empirical[[i]]
  x[x<cutoff_value_inter] <- 0
  interlayer_matrices_empirical[[i]] <- x
}


## @knitr END


# Build Infomap objects -----------------------------------------------
setwd('/media/Data/PLOS_Biol/empirical/')
network_object <- vector(mode = 'list', length = 2)
names(network_object) <- c('intralayer_matrices','interlayer_matrices')
network_object$intralayer_matrices <- intralayer_matrices_empirical
network_object$interlayer_matrices <- interlayer_matrices_empirical

# Get the survival probability from the 002 exp in S

## @knitr interlayer_edge_rescaling
# layers_to_include <- c(118,126,138,142)
layers_to_include <- c(118,126,138,142,154,162)
# surv_prob_S_003 <- NULL
# for (ps in 500:599){
#   print(ps)
#   x <- calculate_rep_survival(ps = ps, scenario = 'S', exp = '003', run = 1, cutoff_prob = 0.85, layers_to_include = layers_to_include)
#   x <- tibble(ps=ps,layer=layers_to_include[1:(length(layers_to_include)-1)],surv_prob=x)
#   # x <- tibble(ps=ps,layer=layers_to_include[1:5],surv_prob=x)
#   surv_prob_S_003 <- rbind(surv_prob_S_003,x)
# }
# png('repertoire_survival_prob.png')
# surv_prob_S_003 %>% ggplot(aes(x=as.factor(layer), group=layer, y=surv_prob))+geom_boxplot()
# dev.off()
# surv_prob_S_003 %<>% group_by(layer) %>% summarise(surv_prob=mean(surv_prob))
# repertoire_survival_prob <- surv_prob_S_003$surv_prob
infomap_empirical <- build_infomap_objects(network_object = network_object,  
                                           write_to_infomap_file = F,
                                           infomap_file_name = 'infomap_empirical.txt', 
                                           return_objects = T,
                                           rescale_by_survival_prob = F,
                                           repertoire_survival_prob = NULL)


infomap_empirical$infomap_interlayer <- rescale_by_divergence('569','S','003',1,layers_to_include,interlayer_edges=infomap_empirical$infomap_interlayer)

quantile(infomap_empirical$infomap_interlayer$w_rescaled,probs = seq(0.8,0.99,0.005))
cutoff_value_inter <- quantile(infomap_empirical$infomap_interlayer$w_rescaled,probs = 0.85)
infomap_empirical$infomap_interlayer %>% ggplot()+
  geom_density(aes(x=w_rescaled),fill='purple',alpha=0.6)+
  geom_density(aes(x=w),fill='gray',alpha=0.6)+
  geom_vline(xintercept = cutoff_value_inter)
infomap_empirical$infomap_interlayer %<>% 
  filter(w_rescaled>=cutoff_value_inter)

# facet_wrap(~layer_s, scale='free')

if (file.exists('infomap_empirical.txt')){unlink('infomap_empirical.txt')}
infomap_empirical$infomap_interlayer %<>% select(layer_s,node_s,layer_t,node_t,w_rescaled) %>% rename(w=w_rescaled)
edges_to_write <- infomap_empirical$infomap_intralayer %>% bind_rows(infomap_empirical$infomap_interlayer)
fwrite(edges_to_write, 'infomap_empirical.txt', sep=' ', col.names = F)





# Run Infomap -------------------------------------------------------------
setwd('/media/Data/PLOS_Biol/empirical/')
system("./Infomap_v01926 infomap_empirical.txt . -i multilayer -d -N 10 --two-level --rawdir --tree --expanded")


# Read Infomap results ----------------------------------------------------
node_list <- infomap_empirical$nodeList
print ('Reading infomap file...')
infomap_file <- 'infomap_empirical_expanded.tree'
lines <- readLines(infomap_file)
cat(lines[1]);cat('\n')
# x <- fread(infomap_file, skip = 2, stringsAsFactors = F) # Read results of Infomap
x <- read_delim(infomap_file, 
                delim = ' ',
                col_types = list(col_character(), col_double(), col_character(), col_integer(), col_integer()), 
                col_names = c('path', 'flow', 'name', 'layer', 'node'),
                skip = 2) # Read results of Infomap
print(x)

x %<>% filter(layer<=length(intralayer_matrices_empirical))

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
modules2$exp <- '003'
modules2$run <- 1
modules2$cutoff_prob <- cutoff_prob_empirical

modules_empirical <- modules2


# Module analysis ---------------------------------------------------------
modsize <- modules_empirical %>% group_by(module) %>% summarise(numStrains=length(strain_unique))

modules_empirical %>% left_join(modsize) %>% ggplot(aes(x=layer,y=module,size=numStrains,color=numStrains))+
  geom_point()

modules_empirical %>% 
  group_by(module) %>% mutate(layer_birth=min(layer),layer_death=max(layer),persistence=layer_death-layer_birth+1) %>% 
  distinct(module,persistence)


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





# Compare edge weight distributions from simulations to empirical --------------------------
# Run once, and then just read results from file:
unlink('/media/Data/PLOS_Biol/empirical/edge_weights_simulated.csv')
edge_weights_simulated <- c()
for (ps in 750:769){
  print(paste('[',Sys.time(),']', ps))
  x <- get_edge_disributions(PS = ps, scenario = 'S', exp = '001', run = 1,cutoff_prob = 0.85, get_inter = F) # The cutoff_prob here is just for file names. it does not really matter because these ar eedge weights before cutoff is applied
  # x <- subset(x, value!=0)
  x <- x[sample(nrow(x),round(nrow(x)*0.2),F),] # Randomly sample 20% of the edge values because there are just too many.
  edge_weights_simulated <- rbind(edge_weights_simulated, x)
  if (ps%%20==0){
    fwrite(edge_weights_simulated[,-c(6,7)], '/media/Data/PLOS_Biol/empirical/edge_weights_simulated_700.csv', append = T)
    edge_weights_simulated <- c()
  }
}
edge_weights_simulated <- fread('/media/Data/PLOS_Biol/empirical/edge_weights_simulated.csv')
edge_weights_simulated <- as.tibble(edge_weights_simulated)
edge_weights_simulated$grp <- 'Simulated'

# See where the quantile passes 0.
quantile(edge_weights_simulated$value, probs = seq(0.85,1,0.005))

png('edge_weights_simulated.png')
ggplot(edge_weights_simulated, aes(x=value))+
  scale_x_continuous(breaks=seq(0,1,0.05))+
  geom_density()
dev.off()

# Work on that:
# edge_weights_simulated %>% 
#   bind_rows(edges_empirical) %>%
#   ggplot(aes(x=value,fill=grp))+
#   geom_density(aes(y=..scaled..))+
#   scale_fill_manual(values=c('#8F25DD','#10A4EF'))+
#   geom_vline(xintercept = cutoff_value, color='red')


# Explore simulated data ------------------------------------------

monitored_variables <- c('prevalence', 'meanMOI','n_circulating_strains', 'n_circulating_genes', 'n_alleles', 'n_total_bites')
exp_colors <- c('gray50','#873600','#16A085')

PS_range <- as.character(500:599)
cases <- expand.grid(ps=PS_range, scenario='S', exp=c('001','002','003'), run=1)
cases$cutoff_prob <- 0.85 # This is just for fime names. there is no cutoff because there are no modules here.
exploratory <- c()
for (i in 1:nrow(cases)){
  print(paste('PS: ',cases$ps[i],' | Scenario: ',cases$scenario[i],' | exp: ',cases$exp[i], ' | run: ',cases$run[i],sep=''))
  tmp <- get_data(parameter_space = cases$ps[i], scenario = cases$scenario[i], experiment = cases$exp[i], run = cases$run[i], cutoff_prob = cases$cutoff_prob[i], use_sqlite = T, tables_to_get = 'summary_general')[[1]]
  exploratory <- rbind(exploratory, tmp)
}
# write_csv(exploratory, '/media/Data/PLOS_Biol/Results/exploratory.csv')
exploratory <- read_csv('/media/Data/PLOS_Biol/Results/exploratory.csv')
exploratory$month <- factor(exploratory$month, levels=c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'))
# exploratory$diversity_grp <- ifelse(exploratory$PS%in%(500:599),'H','L')
layers_df <- exploratory %>% 
  distinct(time,year,month) %>% 
  mutate(layer=1:300) %>% 
  mutate(mylabels=paste(layer,year,month,sep='_'))
exploratory %<>% left_join(layers_df)

# EIR
png('EIR_simulations.png',1600,1000,res = 150)
exploratory %>% 
  # filter(layer%in%layers_to_include) %>% 
  ggplot(aes(x=month, y=EIR, color=exp))+
  scale_color_manual(values = exp_colors)+
  geom_boxplot()+
  stat_summary(aes(group=exp,color=exp), fun.y=mean, geom="line", size=1)+
  facet_grid(~exp)+
  # scale_y_continuous(breaks=seq(0,20,2))+
  mytheme
dev.off()

exploratory %>% 
  filter(layer%in%layers_to_include) %>%
  ggplot(aes(x=layer, y=EIR, color=exp))+
  scale_color_manual(values = exp_colors)+
  geom_boxplot(aes(group=layer))+
  stat_summary(aes(group=exp,color=exp), fun.y=mean, geom="line", size=1)+
  facet_grid(~exp)+
  geom_vline(xintercept = c(131,138.2,145), color='black', linetype='dotted', size=1)+ # IRS. Put 138.2 instead of 138 so lines will not overlap
  scale_x_continuous(breaks = layers_to_include, labels=unique(subset(exploratory_summary, layer%in%layers_to_include)$mylabels))+
  # scale_x_continuous(breaks = 118:162, labels=unique(subset(exploratory_summary, layer%in%118:162)$mylabels))+
  mytheme_no_legend+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size=8),panel.grid = element_blank())


exploratory_summary <- exploratory %>%
  select(-n_infected) %>% 
  # filter(time>time_range[1]&time<time_range[2]) %>%
  gather(variable, value, -pop_id, -year, -month, -time, -layer, -exp, -PS, -scenario, -run, -mylabels) %>% 
  group_by(time, layer, year, month, mylabels, exp, scenario, variable) %>%
  summarise(value_mean=mean(value), value_sd=sd(value)) %>% # Need to average across PS
  filter(variable %in% monitored_variables)
# First make sure that the location of the surveys and IRS in the simulations is correct
png('simulations_bites.png',1600,1000,res = 150)
exploratory_summary %>% 
  filter(layer %in% 110:180) %>%
  filter(variable %in% c('n_total_bites')) %>% 
  ggplot()+
  geom_line(aes(x=layer, y=value_mean, color=exp), size=1)+
  geom_errorbar(aes(x=layer, ymin=value_mean-value_sd,ymax=value_mean+value_sd, color=exp),width=0.2, alpha=0.3)+
  scale_color_manual(values = exp_colors)+
  scale_x_continuous(breaks=110:180, labels = unique(subset(exploratory_summary, layer%in%110:180)$mylabels))+
  geom_vline(xintercept = c(118,126,138,142,154,162), color='black', size=1)+
  geom_vline(xintercept = c(131,138.2,145), color='black', linetype='dotted', size=1)+ # IRS. Put 138.2 instead of 138 so lines will not overlap
  # facet_wrap(~diversity_grp, scales='free')+
  manuscript_theme+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size=8),panel.grid = element_blank())
dev.off()
# Now plot the variables
png('simulations_variables_time_series.png',1600,1000,res = 150)
exploratory_summary %>% 
  filter(layer %in% 110:180) %>%
  filter(variable!='n_total_bites') %>% 
  ggplot()+
  geom_line(aes(x=layer, y=value_mean, color=exp), size=1)+
  geom_errorbar(aes(x=layer, ymin=value_mean-value_sd,ymax=value_mean+value_sd, color=exp),width=0.2, alpha=0.3)+
  scale_color_manual(values = exp_colors)+
  geom_vline(xintercept = c(131,138,145), color='black', linetype='dotted', size=0.5)+ # IRS
  geom_vline(xintercept = c(118,126,138,142,154,162), color='black', size=0.5)+
  facet_wrap(~variable, scales='free')+
  mytheme_no_legend+
  theme(panel.grid = element_blank())
dev.off()


# Number of genes per survey in the simulated vs. empirical data. This is only
# for the subsampled repertoires Do that by going directly to the node list, to
# ensure that the genes I count belong to strains that were fed into Infomap.


var_gene_richness <- exploratory_summary %>% 
  filter(layer %in% c(118,126,138,142,154,162)) %>% 
  filter(exp=='003') %>% 
  filter(variable=='n_circulating_genes') %>% 
  bind_rows(tibble(layer=layers_to_include,variable='gene_richness_empirical',value_mean=var_gene_richness_empirical)) %>% 
  select(layer,variable,value_mean,value_sd) %>% 
  mutate(variable=ifelse(variable=='n_circulating_genes','genes_simulated','genes_empirical'))

ggplot(var_gene_richness, aes(x=as.factor(layer), y=value_mean, group=variable, fill=variable))+
  geom_col(position = 'dodge')+
  geom_errorbar(aes(x=as.factor(layer),ymin=value_mean-value_sd,ymax=value_mean+value_sd),position='dodge')

num_genes <- num_alleles <- c()
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
    distinct(strain_id,gene_id,allele_locus)
  num_genes <- c(num_genes, length(unique(x$gene_id)))
  num_alleles <- c(num_alleles, length(unique(x$allele_locus)))
}



# Cutoff sensitivity in simulated data ------------------------------------


# Create files to run tests on Midway
experiments <- expand.grid(PS=sprintf('%0.3d', 550:569),
                                scen=c('S'),
                                array='1', 
                                layers='118,126,138,142,154,162',
                                exp=c('001','003'),
                                cutoff_prob=0.9,
                                modularity_exp=c(1,4),
                                time = '00:30:00',
                                mem_per_cpu = 4000,
                                stringsAsFactors = F)
nrow(experiments)
make_sbatch_get_data(sbatch_arguments = experiments,
                     make_networks = F,
                     prepare_infomap = T, # make network is nested in this, so no need to call it
                     run_Infomap = T,
                     read_infomap_results = T,
                     temporal_diversity = F,
                     module_Fst = F)


# factorial_design_test <- tibble(option=1:4,unit=c('genes','alleles','alleles','genes'),MOI=c('1','1','All','All'),cutoff_prob=c(0.95,0.85,0.85,0.95))

module_results_simulations <- c()
simulated_edge_weights <- c()
for (i in 1:nrow(experiments)){
  print(paste(i,'/',nrow(experiments)))
  ps <- experiments$PS[i]
  scenario <- experiments$scen[i]
  cutoff_prob <- experiments$cutoff_prob[i]
  exp <- experiments$exp[i]
  modularity_exp <- experiments$modularity_exp[i]
  x <- get_modularity_results(ps,scenario,exp,1,cutoff_prob,folder = paste('/media/Data/PLOS_Biol/Results/',ps,'_',scenario,'_',modularity_exp,'/',sep=''))
  x$modularity_exp <- modularity_exp
  # y <- readLines(paste('/media/Data/PLOS_Biol/Results/',ps,'_',scenario,'_3','/PS',ps,'_',scenario,'_E',exp,'_R1_',cutoff_prob,'_network_info.csv',sep=''))[6]
  # x$cutoff_value <- as.numeric(y)
  module_results_simulations <- rbind(module_results_simulations,x)
  
  x <- read_delim(paste('/media/Data/PLOS_Biol/Results/',ps,'_',scenario,'_',modularity_exp,'/PS',ps,'_S_E',exp,'_R1_',cutoff_prob,'_Infomap_multilayer.txt',sep=''), delim=' ', col_names = c('layer_s','node_s','layer_t','node_t','w'))
  x %<>% filter(layer_s!=layer_t)
  x$PS <- ps
  x$scenario <- scenario
  x$exp <- exp
  x$cutoff_prob <- cutoff_prob
  x$modularity_exp <- modularity_exp
  simulated_edge_weights <- rbind(simulated_edge_weights,x)
}
module_results_simulations <- as.tibble(module_results_simulations)
module_results_simulations %>% group_by(scenario,PS,exp) %>% summarise(n=n())


## Compare edge weight distributions
simulated_edge_weights$grp <- 'Simulated'
names(simulated_edge_weights)[5] <- 'w_rescaled'
empirical_edge_weights <- infomap_empirical$infomap_interlayer
empirical_edge_weights$grp <- 'Empirical'

simulated_edge_weights %>% 
  filter(exp=='003') %>% 
  filter(modularity_exp==1) %>% 
  bind_rows(empirical_edge_weights) %>% 
  filter (layer_s!=layer_t) %>% 
  ggplot(aes(w_rescaled, fill=grp))+
  geom_density()+
  facet_wrap(~layer_s+grp, scales='free')





modsize <- module_results_simulations %>% 
  group_by(scenario,PS,exp,cutoff_prob,module) %>% summarise(size=n())
module_results_simulations %>% 
  filter(exp=='003') %>%
  distinct(PS,scenario,cutoff_prob,module, layer) %>% 
  group_by(PS,scenario,cutoff_prob,module) %>% 
  summarise(birth_layer=min(layer),death_layer=max(layer)+1) %>% 
  left_join(modsize) %>% 
  ggplot(aes(xmin=birth_layer, xmax=death_layer, ymin=module, ymax=module,color=size))+
  geom_rect(size=2)+
  scale_x_continuous(breaks=1:6)+
  scale_color_viridis_c()+ 
  facet_wrap(~cutoff_prob)+
  mytheme

png('/media/Data/PLOS_Biol/empirical/simulated_module_example.png',1600,1000,res = 150)
module_results_simulations %>% 
  # filter(exp=='002') %>% 
  filter(PS==550) %>%
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
  facet_grid(exp~scenario)+
  mytheme
dev.off()


module_results_simulations %>% 
  filter(exp=='003') %>%
  distinct(PS,scenario,cutoff_prob,module, layer) %>% 
  group_by(PS,scenario,cutoff_prob,module) %>% 
  summarise(birth_layer=min(layer),death_layer=max(layer)+1) %>% 
  left_join(modsize) %>% 
  ggplot(aes(xmin=birth_layer, xmax=death_layer, ymin=module, ymax=module,color=size))+
  geom_rect(size=2)+
  scale_x_continuous(breaks=1:6)+
  scale_color_viridis_c()+ 
  facet_wrap(~PS)+
  manuscript_theme


modsize <- module_results_simulations %>% 
  group_by(scenario,PS,exp,cutoff_prob,layer,module) %>% summarise(size=n())
modsize %>%
  filter(PS==550) %>%
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

png('/media/Data/PLOS_Biol/empirical/module_persistence_by_scenario.png',1600,1000,res = 150)
module_persistence %>% 
  group_by(scenario, exp, PS, cutoff_prob) %>% 
  count(persistence) %>% left_join(total_modules) %>% mutate(prop=n/n_modules) %>% 
  # filter(cutoff_prob==0.95) %>% 
  ggplot(aes(y=prop, x=persistence,fill=scenario))+
  facet_wrap(~exp)+
  geom_col(position='dodge')+
  scale_fill_manual(values=c('blue','orange','red'))+
  mytheme
dev.off()

png('/media/Data/PLOS_Biol/empirical/module_persistence_by_exp.png',1600,1000,res = 150)
module_persistence %>% 
  group_by(scenario, exp, PS, cutoff_prob) %>% 
  count(persistence) %>% left_join(total_modules) %>% mutate(prop=n/n_modules) %>% 
  # filter(cutoff_prob==0.95) %>% 
  ggplot(aes(y=prop, x=persistence,fill=exp))+
  facet_wrap(~scenario)+
  geom_col(position='dodge')+
  scale_fill_manual(values=exp_colors)+
  mytheme
dev.off()


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

module_persistence_empirical <- modules_empirical %>% 
  group_by(module) %>% mutate(birth_layer=min(layer),death_layer=max(layer),persistence=death_layer-birth_layer+1) %>% 
  distinct(scenario,exp,PS,module,birth_layer,death_layer,persistence)

module_persistence_empirical$scenario <- 'E'

x <- module_persistence 
# x <- module_persistence %>% bind_rows(module_persistence_empirical)
x <- subset(x, x$birth_layer!=x$death_layer)
x <- subset(x, birth_layer==1)
unique(x$scenario)
unique(x$exp)
x$event <- ifelse(x$death_layer==6,0,1)

fit <- with(x, survfit(Surv(birth_layer, death_layer, event)~scenario+exp))
png('/media/Data/PLOS_Biol/empirical/simulated_module_survival_by_scenario.png',1600,1000,res = 150)
ggsurvplot_facet(fit, data=x, facet.by = 'scenario', conf.int = TRUE)
dev.off()
png('/media/Data/PLOS_Biol/empirical/simulated_module_survival_by_exp.png',1600,1000,res = 150)
ggsurvplot_facet(fit, data=x, facet.by = 'exp', conf.int = TRUE)
dev.off()

# +scale_color_manual(values=c('blue','orange','red'))


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
for (ps in 750:769){
  x <- get_repertoire_persistence_without_modules(ps,'N','001',1,0.85)
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
write_csv(repertoire_persistence_prob,'/media/Data/PLOS_Biol/empirical/repertoire_persistence_prob_700.csv')




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

