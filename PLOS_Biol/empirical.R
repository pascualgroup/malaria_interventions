
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



# Distribution of repertoire persistence ----------------------------------
# In the neutral scenario, where repertoire persistence is purely a result of
# population dynamics, find the probability of persisting EXACTLY t layers. This
# is different than calculating the probablity of persisting at least t layers.
# to calculate the later, we need to look at some cumulative sum of the
# probablities: cumulative=1-cumsum(prob)
persistence_df %>%
  group_by(PS,scenario,exp,run,type) %>%
  count(persistence) %>% 
  mutate(prob=n/sum(n)) %>% 
  ggplot()+
  geom_line(aes(x=persistence, y=prob, color=type, linetype=type),size=1)+ geom_point(aes(x=persistence, y=prob),color='red')+scale_y_continuous(limits = c(0,1))



# Analyzing simulations ---------------------------------------------------
# !!! STILL OLD CODE !!!
# Randomly select individuals to match to the data. This is where we limit to MOI=1
empirical_number_individuals <- c(97,67,66,50)
for (i in 1:length(layersFull)){
  x <-layersFull[[i]]
  if (nrow(x) < empirical_number_individuals[i]){next}
  sampled_individuals <- sample(rownames(x), empirical_number_individuals[i], replace = F)
  x <- x[sampled_individuals,sampled_individuals]
  layersFull[[i]] <- x
}
print('MATCHED layer sizes:')
sapply(layersFull, nrow)

print('layer density:')
sapply(layersFull, function (x) gdensity(x,F,T))

# Extract the similarity matrix among all the repertoires in the 4 layers
similarityMatrix <- overlapAlleleAdj(table(strainComposition_reduced$strainId, strainComposition_reduced$geneId))
nodeLabel <- sort(unique(unlist(lapply(layersFull,rownames))))
similarityMatrix <- similarityMatrix[nodeLabel, nodeLabel]
diag(similarityMatrix) <- 0
write.csv(similarityMatrix, paste('../Results/',filename_base,'_similarityMatrix_no_cutoff.csv',sep='')) # save this so to be able tro plt edge weight distributons fast

# Apply cutoff: percentile of strongest links -------------------------------------------------
# calculates a cutoff value which will be used across all layers, and interlayer edges

if (cutoffMethod=='CP'){
  message('Applying cutoff by percentile...')
  cutoffValue <- quantile(as.vector(similarityMatrix[similarityMatrix!=0]), probs = cutoffPercentile)
  # ggHistogram(as.vector(similarityMatrix[similarityMatrix!=0]), xlab = 'Non-zero similarity values')+geom_vline(xintercept = cutoffValue, color='red')
  similarityMatrix[similarityMatrix<cutoffValue] <- 0 # This is to be used later for the interlayer edges
  # Apply cutoff on the layers
  for (i in 1:length(layersFull)){
    x <-layersFull[[i]]
    x[x<cutoffValue] <- 0
    layersFull[[i]] <- x
  }
  
  print('Layer density by applying a cutoff percentile:')
  sapply(layersFull, function (x) gdensity(x,F,T))
  write.csv(sapply(layersFull, function (x) gdensity(x,F,T)), paste('../Results/',filename_base,'_layerdensity.csv',sep=''))
}

# Apply cutoff: match density ---------------------------------------------
if (cutoffMethod=='MD'){
  mprint('Applying cut off by matching densities')
  
  density_vector_intra <- str_replace(str_replace(args[5],pattern = '\\[',''),'\\]','')
  density_vector_intra <- as.numeric(unlist(strsplit(density_vector_intra, "[,]")))
  density_vector_inter <- str_replace(str_replace(args[6],pattern = '\\[',''),'\\]','')
  density_vector_inter <- as.numeric(unlist(strsplit(density_vector_inter, "[,]")))
  
  for (i in 1:length(layersFull)){
    x <- match_density(layersFull[[i]],density_vector_intra[i])
    layersFull[[i]] <- x
  }
  
  print('Matched intra-layer densities:')
  print(sapply(layersFull, function (x) gdensity(x,F,T)))
  write.csv(sapply(layersFull, function (x) gdensity(x,F,T)), paste('../Results/',filename_base,'_layerdensity.csv',sep=''))
}


# Write files for Infomap -------------------------------------------------
message('Building Infomap files...')
nodeLabel <- sort(unique(unlist(lapply(layersFull,rownames))))
nodeList <- data.frame(nodeID=1:length(nodeLabel), nodeLabel)
# Build interlayer edges
ile_matrix <- list()
ile_edgelist <- list()

# Get the probabilities of persistence for 4,8 and 12 months
persistence_dist <-  data.table::fread(paste('neutral_',simNum,'_',cutoffPercentile,'_',cutoffMethod,'_repertoire_persistence.csv',sep=''))

get_persistence <- function(q){
  if (q%in%persistence_dist$persistence){
    rep_persistence <- subset(persistence_dist, persistence == q)$cumulative
  } else { # If q is not in the least take the closesty value
    mprint('Cannot find persistence value!')
    x <- persistence_dist$persistence
    dt = data.table(x, val = x) # you'll see why val is needed in a sec
    setattr(dt, "sorted", "x")  # let data.table know that w is sorted
    setkey(dt, x) # sorts the data
    # binary search and "roll" to the nearest neighbour
    # In the final expression the val column will have the you're looking for.
    closest <- dt[J(q), roll = "nearest"]$val
    rep_persistence <- subset(persistence_dist, persistence == closest)$cumulative
  }
  return(rep_persistence)
}


for (t in 1:(length(timeSlices)-1)){
  cat('Building interlayer edges for layers: ',t,'-->',t+1,'\n')
  strainCopies_t <- rownames(layersFull[[t]])
  strainCopies_t1 <- rownames(layersFull[[t+1]])
  if (length(strainCopies_t)<2 | length(strainCopies_t1)<2){
    print('need minimum of 2 strains in t and t+1 to build a matrix. SKIPPING layer')
    next
  } # need minimum of 2 strains in t and t+1 to build a matrix
  
  m <- similarityMatrix[strainCopies_t,strainCopies_t1]
  print(paste('density of interlayer edges layer',t,'to',t+1,':',gdensity(m,F,T)))
  if (cutoffMethod=='MD'){
    m <- match_density(m,density_vector_inter[t],bipartite = T)
    print(paste('MATCHED density of interlayer edges layer',t,'to',t+1,':',gdensity(m,F,T)))
  }
  
  # Divide edge weights by the probability of persistence. those that persisted for longer will have stronger values
  if (t==1){m <- m/get_persistence(8)} # 8 months between S1 and S2
  if (t==2){m <- m/get_persistence(12)} # 12 months between S2 and S3
  if (t==3){m <- m/get_persistence(4)} # 4 months between S3 and S4
  
  if(sum(m)==0){
    print('No interactions, SKIPPING layer')
    next
  }
  ile_matrix[[t]] <- m
  
  NZ <- which(m!=0, arr.ind = T) # non-zero elements
  
  edges_interlayer <- data.frame(layer_s=rep(t,nrow(NZ)),
                                 node_s=rep(NA,nrow(NZ)),
                                 layer_t=rep(t+1,nrow(NZ)),
                                 node_t=rep(NA,nrow(NZ)),
                                 w=rep(NA,nrow(NZ)))
  
  edges_interlayer$node_s <- rownames(m)[NZ[,1]]
  edges_interlayer$node_t <- colnames(m)[NZ[,2]]
  for (n in 1:nrow(edges_interlayer)){
    edges_interlayer[n,'w'] <- m[as.character(edges_interlayer[n,'node_s']),as.character(edges_interlayer[n,'node_t'])]
  }
  
  edges_interlayer$node_s <- nodeList$nodeID[match(edges_interlayer$node_s,nodeList$nodeLabel)]
  edges_interlayer$node_t <- nodeList$nodeID[match(edges_interlayer$node_t,nodeList$nodeLabel)]
  
  ile_edgelist[[t]] <- edges_interlayer
}
inter_edges_simulated <- do.call(rbind.data.frame,ile_edgelist)
intra_edges_simulated <- infomap_makeIntralayerEdges(layersFull,nodeList)

print('inter-layer densities:')
print(sapply(ile_matrix, dim))
print(sapply(ile_matrix, function (x) gdensity(x,T)))
write.csv(sapply(ile_matrix, function (x) gdensity(x,T)), paste('../Results/',filename_base,'_inter-layerdensity.csv',sep=''))


## Write file for infomap
message('Writing Infomap files')
file <- paste(filename_base,'_Infomap_multilayer','.txt',sep='')
print(paste('Infomap file:',file))
if (file.exists(file)){unlink(file)}
sink(file, append = T)
cat("# A network in a general multiplex format");cat('\n')
cat(paste("*Vertices",nrow(nodeList)));cat('\n')
write.table(nodeList, file, append = T,sep=' ', quote = T, row.names = F, col.names = F)
cat("*Multiplex");cat('\n')
cat("# layer node layer node [weight]");cat('\n')
cat("# Intralayer edges");cat('\n')
write.table(intra_edges_simulated, file, sep = ' ', row.names = F, col.names = F, quote = F, append = T)
cat("# Interlayer edges");cat('\n')
write.table(inter_edges_simulated, file, sep = ' ', row.names = F, col.names = F, quote = F, append = T)
sink.reset()


message('Running Infomap...')
#runInfomapCommand <- paste('./Infomap ',file,' . -2 -i multiplex --multiplex-relax-rate -1 -d -N 10 --rawdir --tree --expanded',sep='')
runInfomapCommand <- paste('./Infomap_v01914 ',file,' ',args[7],' . -2 -i multilayer -d -N 20 --rawdir --tree --expanded',sep='')

system(runInfomapCommand)

message('Reading Infomap results...')

modules_simulated <- infomap_readTreeFile(paste(filename_base,'_Infomap_multilayer','_expanded.tree',sep=''), reorganize_modules = T, remove_buggy_instances = T,max_layers = 4)

modules_simulated$scenario <- scenario
modules_simulated$simNum <- simNum

write.csv(modules_simulated, paste('../Results/',filename_base,'_modules.csv',sep=''))
