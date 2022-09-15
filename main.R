library("rjson")
library("dplyr")
library("igraph")
library("rlist")
library("tidyverse")
library("useful")
source("config.R")
source("tool.R")

# read from input file and data file
ListOfMassFeatures <- read.table("testdata0710.txt",header = TRUE)
human.metabolic <- fromJSON(file = "jsonformatter.json")$human_model_mfn

ListOfMassFeatures <- tibble(ListOfMassFeatures)
ListOfMassFeatures <- ListOfMassFeatures %>% mutate(raw_row_number = row_number()) %>% filter(m.z > MASS_RANGE[1] & m.z < MASS_RANGE[2])

# get max retention time and mz
max_retention_time = max(ListOfMassFeatures$retention_time)
max_mz = max(ListOfMassFeatures$m.z)

# add retention rank mark the outflow time
ListOfMassFeatures <-ListOfMassFeatures %>% arrange(retention_time) %>% mutate(retention_rank = row_number())

ListOfMassFeatures <-ListOfMassFeatures %>% arrange(raw_row_number)

# find the significant set
p_hotspots <- c(0.2, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0001)
N_hotspots <- sapply(p_hotspots,function(x){sum(ListOfMassFeatures$p.value<x,na.rm = TRUE)})

N_optimum = 300
N_minimum = 30
N_quantile = nrow(ListOfMassFeatures) / 4
chosen = 9999

for(i in length(N_hotspots)){
  if (N_hotspots[i] < N_quantile & N_hotspots[i] > N_minimum)
    chosen = i
}

ListOfMassFeatures <- ListOfMassFeatures %>% mutate(is_significant = ifelse(p.value<p_hotspots[chosen],TRUE,FALSE))

# List of row number
input_feature_list = (ListOfMassFeatures %>% filter(is_significant == TRUE))$raw_row_number

# tackle with data file
# human_data_cpd = as.undirected(read_graph("graph",format = "graphml"))
human_data_cpd = as.undirected(read_graph("graph_gml",format = "gml"))
total_cpd_list = get.vertex.attribute(human_data_cpd)
V(human_data_cpd)$name = get.vertex.attribute(human_data_cpd)$label

# Data Clean
rtime_tolerance <- max_retention_time * RETENTION_TIME_TOLERANCE_FRAC
rtime_tolerance_rank <- nrow(ListOfMassFeatures) * RETENTION_TIME_TOLERANCE_FRAC

wanted_ions <- pos_default



compounds <- human.metabolic$Compounds

# ion_cpd_tree
# initialization
ion_cpd_tree<-vector("list",2000)

for(i in 1:length(compounds)){
  if (compounds[[i]]$mw){
    for(j in 1:length(compounds[[i]]$adducts)){
      ion = names(compounds[[i]]$adducts)[j]
      mass = as.double(compounds[[i]]$adducts[names(compounds[[i]]$adducts)[j]])
      if(ion %in% pos_default & mass > MASS_RANGE[1] & mass < MASS_RANGE[2]) {
        ion_cpd_tree[[as.integer(mass)]][[length(ion_cpd_tree[[as.integer(mass)]]) + 1]] <- c(names(compounds)[i],ion,mass)
      }
    }
  }
}


# row dict
row_dict <- list()
a = 0
for (i in 1:nrow(ListOfMassFeatures)) {
  row_dict[[ListOfMassFeatures$raw_row_number[i]]] = ListOfMassFeatures[i,]
}


### List of Empirical Compound
# __match_to_mzFeatures__()
mztol = ListOfMassFeatures$m.z * 0.000010 # mz tolerance
mzint = as.integer(ListOfMassFeatures$m.z)
ListOfMassFeatures <-  ListOfMassFeatures %>% add_column(matched = NA)
for(k in 1:nrow(ListOfMassFeatures)) {
  match <- list()
  for (i in c(mzint[k] - 1, mzint[k], mzint[k] + 1)){
    for (j in ion_cpd_tree[[i]]){
      if (is.null(j)) next # if ion tree is null
      
      if(abs(as.double(j[[3]]) - as.double(ListOfMassFeatures$m.z[k])) < mztol[k]){
        match[[length(match) + 1]] = j
      }
    }
  }
  if(!is_empty(match)) ListOfMassFeatures$matched[k] = list(match) # List 
}


# index_Compounds_to_mzFeatures(self):
# cpds_to_mz_features [cpds: name,mz,retention time, raw row number]
# one element missing !!!!!!!!!!!!!!!!!!!!!!!1
cpds_to_mz_features = data.frame(cpds_id=character(),
                                 adduct=character(), 
                                 mz=double(), 
                                 retention_time=double(),
                                 raw_row_number=integer())
cpds_to_mz_features <- tibble(cpds_to_mz_features)

for(i in 1:nrow(ListOfMassFeatures)){
  if (!is.na(ListOfMassFeatures$matched[[i]][1])){
    for(j in 1:length(ListOfMassFeatures$matched[[i]])){
      cpds_id = ListOfMassFeatures$matched[[i]][[j]][1]
      adduct = ListOfMassFeatures$matched[[i]][[j]][2]
      mz = as.double(ListOfMassFeatures$matched[[i]][[j]][3])
      retention_time = as.double(ListOfMassFeatures$retention_time[i])
      raw_row_number = as.integer(i)
      # cpds_to_mz_features[[cpds]][[length(cpds_to_mz_features[[cpds]]) + 1]] <- 
      # c(adduct,mz,ListOfMassFeatures$retention_time[i],ListOfMassFeatures$raw_row_number[i])
      cpds_to_mz_features <- cpds_to_mz_features %>% add_row(cpds_id,adduct,mz,retention_time,raw_row_number)
      
    }
    
  }
}

# group adduct of a certain compound by their retentioin time
# __split_Compound__(self, compoundID, list_match_mzFeatures)
group_name = cpds_to_mz_features %>% group_by(cpds_id) %>% summarise()





#Ecompounds (get_user_data line 235)
ecompounds <- data.frame(matrix(ncol = 5, nrow = 0))
# #provide column names
colnames(ecompounds) <- c('str_row_ion', 'compounds','massfeature_rows','ions','coelute_group')

ecompounds$str_row_ion <- as.character(ecompounds$str_row_ion)
ecompounds$compounds <- as.character(ecompounds$compounds)
ecompounds$massfeature_rows <- as.integer(ecompounds$massfeature_rows)
ecompounds$ions <- as.character(ecompounds$ions)
ecompounds$coelute_group <- as.integer(ecompounds$coelute_group)
# ecompounds = data.frame(str_row_ion=character(),
#                         compounds=list(),
#                         massfeature_rows=list(),
#                         ions = list())
ecompounds <- tibble(ecompounds)


make_str_ion_row <- function(x){
  res_make_str_ion_row = character()
  # for(i in 1:nrow(x)){
  #   if (i == 1) {
  res_make_str_ion_row <- paste(x$raw_row_number[1],x$adduct[1],sep = "_")
  #   } else{
  #     res_make_str_ion_row <- paste(res_make_str_ion_row, paste(x$raw_row_number[i],x$adduct[i],sep = "_"), sep = ";")
  #   }
  # }
  return(res_make_str_ion_row)
}

# (ecompounds,a line of adducts, coelute_group_count)
add_to_ecompounds <- function(x,y,z){
  x <- x %>% add_row(str_row_ion = make_str_ion_row(y),
                                       compounds = y$cpds_id,
                                       massfeature_rows = y$raw_row_number,
                                       ions = y$adduct,
                                       coelute_group = z)
  return(x)
}

# add 2 str_row_ion!!!!!!!
coelute_group_count = 1
for(i in 1:nrow(group_name)) {
  adducts <- cpds_to_mz_features %>% filter(cpds_id == toString(group_name[i,])) %>% arrange(retention_time)
  # tmp <- adducts[1,]
  ecompounds <- add_to_ecompounds(ecompounds,adducts[1,],coelute_group_count)
  j = 1
  while(j < nrow(adducts)){
    a = ListOfMassFeatures[adducts$raw_row_number[j],]
    b = ListOfMassFeatures[adducts$raw_row_number[j + 1],]

    # group adducts with similar retention time(coelution) into a tmp
    if((abs(a$retention_time - b$retention_time) < rtime_tolerance) | (abs(a$retention_rank - b$retention_rank) < rtime_tolerance_rank)){
      ecompounds <- add_to_ecompounds(ecompounds,adducts[j+1,],coelute_group_count)    
      # tmp <- tmp %>% add_row(adducts[j + 1,])  
    } else {
      # ecompounds <- ecompounds %>% add_row(str_row_ion = make_str_ion_row(tmp),compounds = tmp$cpds_id,massfeature_rows = list(tmp$raw_row_number),ions = list(tmp$adduct))
      # tmp <- adducts[j + 1,]
      coelute_group_count = coelute_group_count + 1
      ecompounds <- add_to_ecompounds(ecompounds,adducts[j+1,],coelute_group_count)
    }
    j = j + 1
  }
  coelute_group_count = coelute_group_count + 1
  # ecompounds <- ecompounds %>% add_row(str_row_ion = make_str_ion_row(tmp),compounds = list(unique(tmp$cpds_id)),massfeature_rows = list(tmp$raw_row_number),ions = list(tmp$adduct))
}
ecompounds %>% group_by(coelute_group) %>% summarise()

# merge the ecompounds

# ecompounds %>% group_by(str_row_ion) %>% summarise()
# 
# ecompounds <- ecompounds %>% ungroup()

# get the empirical compound
ecompounds <- ecompounds %>% add_column(str_row_ions = NA)
 
group_sum <- ecompounds %>% group_by(coelute_group) %>% summarise(sris = toString(list(str_row_ion)))

for(i in 1:nrow(ecompounds)){
  ecompounds$str_row_ions[i] <- (group_sum %>% filter(coelute_group == ecompounds$coelute_group[i]))$sris
}


ecompounds <- ecompounds %>% add_column(EID = NA)
ecompounds_names = ecompounds %>% group_by(coelute_group) %>% group_by(str_row_ions)  %>% summarise()
for(i in 1:nrow(ecompounds_names)){
  ecompounds[which(ecompounds$str_row_ions == as.character(ecompounds_names[i,])),]$EID = i
}

cal_ev_score <- function(x){
  temp_score = 0
  for(i in 1:length(x)) {
    temp_score = temp_score + dict_weight_adduct[x[[i]]]
  }
  return(temp_score)
}
# add ion column
# ecompounds_v2 <- ecompounds %>% rowwise %>% group_by(EID) %>%  summarize(evidence_score = cal_ev_score(ions),
#                                                              primary_ion_present = ifelse(as.logical(length(intersect(ions,primary_ions))),TRUE,FALSE),
#                                                              ions = list(unique(ions)),
#                                                              row = list(unique(massfeature_rows)),
# )%>% ungroup()
ecompounds_v2 <- ecompounds %>% group_by(EID) %>%  summarize(evidence_score = cal_ev_score(ions),
                                                                           primary_ion_present = ifelse(as.logical(length(intersect(ions,primary_ions))),TRUE,FALSE),
                                                                          ions = ions,
                                                                          row = massfeature_rows,
                                                                          compounds = compounds
                                                                          )%>% ungroup()
ecompounds_v2 <- ecompounds_v2 %>% filter(primary_ion_present == TRUE) %>% select(-primary_ion_present)
### reference list
ref_rows = ListOfMassFeatures$raw_row_number

### significant list
sig_rows = (ListOfMassFeatures %>% filter(is_significant == TRUE))$raw_row_number

trio <- ecompounds_v2 %>% filter(row %in% intersect(ecompounds_v2$row, sig_rows))  



# =========================================================================================================
# work on module analysis

## Get Modules
seeds = trio$compounds
seed_cpds = trio$compounds



modules <- tibble()
modules <- modules %>% mutate(subgraph = NA, activity_score = double())
modules$subgraph <- as.list(modules$subgraph)


for (j in 1:SEARCH_STEPS){
  cpds_edges = data.frame(temp = as_ids(E(human_data_cpd)[from(seeds)]))
  cpds_edges = as_tibble(cpds_edges)
  cpds_edges = separate(data = cpds_edges, col = temp, into = c("V1", "V2"), sep = "\\|")
  cpds_edges <- cpds_edges %>% mutate(exist_or_not = NA)
  if (j == 1){
    for(i in 1:nrow(cpds_edges)){
      cpds_edges$exist_or_not[i] = (as.logical(length(intersect(cpds_edges$V1[i],seeds))) & as.logical(length(intersect(cpds_edges$V2[i],seeds))) ) 
      sig_edges <- cpds_edges %>% filter(exist_or_not == TRUE) %>% select(-exist_or_not)
      sig_graph <- as.undirected(graph_from_edgelist(as.matrix(sig_edges)))
    }
  } else {
    sig_graph <- as.undirected(graph_from_edgelist(as.matrix(cpds_edges %>% select(-exist_or_not))))
    seeds = V(sig_graph)$name
  }
  for(i in decompose.graph(sig_graph)){
    if (length(V(i)) > 3 & length(V(i)) < MODULE_SIZE_LIMIT) {
      
      # shave
      excessive = names(which(degree(i,setdiff(V(i)$name,seed_cpds)) == 1))
      while(length(excessive)) {
        i = delete_vertices(i, excessive)
        excessive = names(which(degree(i,setdiff(V(i)$name,seed_cpds)) == 1))
      }

      # calculate activity score
      # get empirical cpds number
      # just use EID
      # Ns
      emp_cpds_num = nrow(unique(trio %>% filter(compounds %in% intersect(V(i)$name,seed_cpds)) %>% select(EID)))
      #Nm
      nodes_num = as.double(length(V(i)))
      Q = 0
      if(nodes_num > 0){
        # m 
        ref_edge_num = length(E(human_data_cpd))
        expected = 0
        for (ii in V(i)$name){
          for (jj in V(i)$name){
            if (ii != jj) {
              expected = expected + as.numeric(degree(human_data_cpd,ii)) * as.numeric(degree(human_data_cpd,jj))
            } 
          }
        }
        expected = expected / (4.0 * ref_edge_num)
        Q = (length(E(i)) - expected) / ref_edge_num
        as = sqrt(length(seed_cpds) / nodes_num) * Q * (emp_cpds_num / nodes_num) * 100
      }else{
        as = 0
        Q = 0
      }
      # filter duplicates
      exist_in_module = FALSE
      if(nrow(modules) > 0){
        for (md_row in 1:nrow(modules)){
          if ((length(intersect(V(i)$name, V(modules$subgraph[[md_row]])$name)) == length(V(i)$name)) & (length(V(i)$name) == length(V(modules$subgraph[[md_row]])$name)))
            exist_in_module = TRUE
        }
      }
      if(!exist_in_module)
        modules <- modules %>% add_row(subgraph = list(i), activity_score = as)
    }
  }
  # if (j == 1) {
  #   seeds1 = V(sig_graph)$name
  # } else if (j == 2){
  #   seeds2 = V(sig_graph)$name
  # } else if (j == 3) {
  #   seeds3 = V(sig_graph)$name
  # } else{
  #   seeds4 = V(sig_graph)$name
  # }
}

# model2 use Newman's spectral split method
dQs = list()
for(i in 1:nrow(modules)){
  j = modules$subgraph[[i]]
  if(length(V(j)) > 5){
    num_edges = length(E(j))
    num_vertex = length(V(j))
    
    # make adjacency matrix
    am = matrix(0,length(V(j)),length(V(j)))
    md2_edgelist = get.edgelist(j)
    for(e in 1:nrow(md2_edgelist)){
      am[which(V(j)$name == md2_edgelist[e,1]),which(V(j)$name == md2_edgelist[e,2])] = 1 + am[which(V(j)$name == md2_edgelist[e,1]),which(V(j)$name == md2_edgelist[e,2])]
      am[which(V(j)$name == md2_edgelist[e,2]),which(V(j)$name == md2_edgelist[e,1])] = 1 + am[which(V(j)$name == md2_edgelist[e,2]),which(V(j)$name == md2_edgelist[e,1])]
    }
    
    # compute degrees
    md2_degree = apply(am, 1, sum)
    
    # make modularity matrix
    mm = matrix(0,length(V(j)),length(V(j)))
    for(ii in 1:num_vertex){
      for(jj in 1:num_vertex){
        mm[ii,jj] = am[ii,jj] - md2_degree[ii] * md2_degree[jj]/(2.0 * num_edges)
      }
    }
    
    
    divisible_list = list()
    divisible_list[[1]] = V(j)$name
    md2_res = list()
    # md2_res[[1]] = V(j)$name
    while(length(divisible_list)){
      remove_list = c()
      for(kkk in 1:length(divisible_list)){
        
        
        kk = sort(divisible_list[[kkk]], decreasing = TRUE)
        num_vertex_temp = length(kk)
        mm_temp = matrix(0,num_vertex_temp,num_vertex_temp)
        for(ii in 1:num_vertex_temp){
          for(jj in 1:num_vertex_temp){
            mm_temp[ii,jj] = mm[which(V(j)$name == kk[ii]),which(V(j)$name == kk[jj])]
          }
        }
        # 0 * 0 matrix
        if(!nrow(mm_temp)) {
          # divisible_list = divisible_list[-kkk]
          remove_list = c(remove_list,kkk)
          next
        }
        # eigen divide
        eg <- eigen(mm_temp)
        eigen_values = eg$values
        eigen_vectors = eg$vectors
        
        if(max(eigen_values) > 0) {
          lead_vector = eigen_vectors[,which(eigen_values == max(eigen_values))]
          group1 = c()
          group2 = c()
          for(ii in 1:num_vertex_temp){
            if(lead_vector[ii] < 0) group1 = c(group1, kk[ii])
            else group2 = c(group2, kk[ii])
          }
          #calculate s
          md2_s = replicate(num_vertex_temp, 0)
          md2_s[which(kk %in% group1)] = 1
          md2_s[which(kk %in% group2)] = -1
          
          # compute deltaQ
          dQ = 0
          for(ii in 1:num_vertex_temp){
            for(jj in 1:num_vertex_temp){
              dQ = dQ + (md2_s[ii] * md2_s[jj] - 1) * mm_temp[ii,jj]
            }
          }
          dQ = dQ / (4.0 * num_edges)
          dQs[[length(dQs) + 1]] = dQ
        }
        
        
        
        
        
        # remove from temp list
        # divisible_list = divisible_list[-kkk]
        remove_list = c(remove_list,kkk)
        if(dQ > 0 & length(group1) > 0) {
          # #remove from res list
          # md2_res = md2_res[unlist(md2_res) != unlist(kk)]
          # md2_res = discard(md2_res, ~ all(is.null(.x)))
          # add to temp list
          if(!is.null(group1)) divisible_list[[length(divisible_list) + 1]] = group1
          if(!is.null(group2)) divisible_list[[length(divisible_list) + 1]] = group2

          
          
        } else{
          # add to res list
          if(!is.null(group1))  md2_res[[length(md2_res) + 1]] = group1
          if(!is.null(group2))  md2_res[[length(md2_res) + 1]] = group2
        }
      }
      divisible_list = divisible_list[-remove_list]
    }
    md2_subgraph = list()
    for(i in 1:length(md2_res)){
      if (length(md2_res[[i]]) > 3) {
        md2_subgraph[[length(md2_subgraph) + 1]] = induced_subgraph(human_data_cpd,md2_res[[i]])
      }
    }
    
    for(i in md2_subgraph){
        
      # shave
      excessive = names(which(degree(i,setdiff(V(i)$name,seed_cpds)) == 1))
      while(length(excessive)) {
        i = delete_vertices(i, excessive)
        excessive = names(which(degree(i,setdiff(V(i)$name,seed_cpds)) == 1))
      }
        
        # calculate activity score
        # get empirical cpds number
        # just use EID
        # Ns
        emp_cpds_num = nrow(unique(trio %>% filter(compounds %in% intersect(V(i)$name,seed_cpds)) %>% select(EID)))
        #Nm
        nodes_num = as.double(length(V(i)))
        Q = 0
        if(nodes_num > 0){
          # m 
          ref_edge_num = length(E(human_data_cpd))
          expected = 0
          for (ii in V(i)$name){
            for (jj in V(i)$name){
              if (ii != jj) {
                expected = expected + as.numeric(degree(human_data_cpd,ii)) * as.numeric(degree(human_data_cpd,jj))
              } 
            }
          }
          expected = expected / (4.0 * ref_edge_num)
          Q = (length(E(i)) - expected) / ref_edge_num
          as = sqrt(length(seed_cpds) / nodes_num) * Q * (emp_cpds_num / nodes_num) * 100
        }else{
          as = 0
          Q = 0
        }
        # filter duplicates
        exist_in_module = FALSE
        if(nrow(modules) > 0){
          for (md_row in 1:nrow(modules)){
            if ((length(intersect(V(i)$name, V(modules$subgraph[[md_row]])$name)) == length(V(i)$name)) & (length(V(i)$name) == length(V(modules$subgraph[[md_row]])$name)))
              exist_in_module = TRUE
          }
        }
        if(!exist_in_module)
          modules <- modules %>% add_row(subgraph = list(i), activity_score = as)
    }
    
  }
}

# filter th modules
temp = c()
for(i in 1:nrow(modules)){
  if(length(V(modules$subgraph[[i]])) <= 3){
    temp <- c(temp,i)
  }
}
if (!is.null(temp))
  sig_modules <- modules[-temp,]






# do permutations
N = length(input_feature_list)
permuation_mscores = c()
temp_record = c()
for(iii in 1:100){
  ran_sample = sample(ref_rows,N)
  ref_trio <- ecompounds_v2 %>% filter(row %in% intersect(ecompounds_v2$row, ran_sample))
  
  
  
  
  # copy
  seeds = unique(ref_trio$compounds)
  seed_cpds = unique(ref_trio$compounds)
  
  # the graph is not the same as compounds
  while(length(intersect(V(human_data_cpd)$name,seeds)) != length(seeds)){
    ran_sample = sample(ref_rows,N)
    ref_trio <- ecompounds_v2 %>% filter(row %in% intersect(ecompounds_v2$row, ran_sample))
    

    seeds = unique(ref_trio$compounds)
    seed_cpds = unique(ref_trio$compounds)
  }
  
  
  modules <- tibble()
  modules <- modules %>% mutate(subgraph = NA, activity_score = double())
  modules$subgraph <- as.list(modules$subgraph)
  
  
  for (j in 1:SEARCH_STEPS){
    cpds_edges = data.frame(temp = as_ids(E(human_data_cpd)[from(seeds)]))
    cpds_edges = as_tibble(cpds_edges)
    cpds_edges = separate(data = cpds_edges, col = temp, into = c("V1", "V2"), sep = "\\|")
    cpds_edges <- cpds_edges %>% mutate(exist_or_not = NA)
    if (j == 1){
      for(i in 1:nrow(cpds_edges)){
        cpds_edges$exist_or_not[i] = (as.logical(length(intersect(cpds_edges$V1[i],seeds))) & as.logical(length(intersect(cpds_edges$V2[i],seeds))) ) 
        sig_edges <- cpds_edges %>% filter(exist_or_not == TRUE) %>% select(-exist_or_not)
        sig_graph <- as.undirected(graph_from_edgelist(as.matrix(sig_edges)))
      }
    } else {
      sig_graph <- as.undirected(graph_from_edgelist(as.matrix(cpds_edges %>% select(-exist_or_not))))
      seeds = V(sig_graph)$name
    }
    
    
    for(i in decompose.graph(sig_graph)){
      if (length(V(i)) > 3 & length(V(i)) < MODULE_SIZE_LIMIT) {
        # shave
        excessive = names(which(degree(i,setdiff(V(i)$name,seed_cpds)) == 1))
        while(length(excessive)) {
          i = delete_vertices(i, excessive)
          excessive = names(which(degree(i,setdiff(V(i)$name,seed_cpds)) == 1))
        }
        
        # calculate activity score
        # get empirical cpds number
        # just use EID
        # Ns
        emp_cpds_num = nrow(unique(ref_trio %>% filter(compounds %in% intersect(V(i)$name,seed_cpds)) %>% select(EID)))
        #Nm
        nodes_num = as.double(length(V(i)))
        Q = 0
        if(nodes_num > 0){
          # m 
          ref_edge_num = length(E(human_data_cpd))
          expected = 0
          for (ii in V(i)$name){
            for (jj in V(i)$name){
              if (ii != jj) {
                expected = expected + as.numeric(degree(human_data_cpd,ii)) * as.numeric(degree(human_data_cpd,jj))
              } 
            }
          }
          expected = expected / (4.0 * ref_edge_num)
          Q = (length(E(i)) - expected) / ref_edge_num
          as = sqrt(length(seed_cpds) / nodes_num) * Q * (emp_cpds_num / nodes_num) * 100
        }else{
          as = 0
          Q = 0
        }
        
        # filter duplicates
        exist_in_module = FALSE
        if(nrow(modules) > 0){
          for (md_row in 1:nrow(modules)){
            if ((length(intersect(V(i)$name, V(modules$subgraph[[md_row]])$name)) == length(V(i)$name)) & (length(V(i)$name) == length(V(modules$subgraph[[md_row]])$name)))
              exist_in_module = TRUE
          }
        }
        if(!exist_in_module)
          modules <- modules %>% add_row(subgraph = list(i), activity_score = as)
      }
    }
  }
  modules1<-modules
  
  # model2 use Newman's spectral split method
  dQs = list()
  for(i in 1:nrow(modules)){
    j = modules$subgraph[[i]]
    if(length(V(j)) > 5 & length(E(j) >= 1)){
      num_edges = length(E(j))
      num_vertex = length(V(j))
      
      # make adjacency matrix
      am = matrix(0,length(V(j)),length(V(j)))
      md2_edgelist = get.edgelist(j)
      for(e in 1:nrow(md2_edgelist)){
        am[which(V(j)$name == md2_edgelist[e,1]),which(V(j)$name == md2_edgelist[e,2])] = 1 + am[which(V(j)$name == md2_edgelist[e,1]),which(V(j)$name == md2_edgelist[e,2])]
        am[which(V(j)$name == md2_edgelist[e,2]),which(V(j)$name == md2_edgelist[e,1])] = 1 + am[which(V(j)$name == md2_edgelist[e,2]),which(V(j)$name == md2_edgelist[e,1])]
      }
      
      # compute degrees
      md2_degree = apply(am, 1, sum)
      
      # make modularity matrix
      mm = matrix(0,length(V(j)),length(V(j)))
      for(ii in 1:num_vertex){
        for(jj in 1:num_vertex){
          mm[ii,jj] = am[ii,jj] - md2_degree[ii] * md2_degree[jj]/(2.0 * num_edges)
        }
      }
      
      
      divisible_list = list()
      divisible_list[[1]] = V(j)$name
      md2_res = list()
      # md2_res[[1]] = V(j)$name
      while(length(divisible_list)){
        remove_list = c()
        for(kkk in 1:length(divisible_list)){
          
          
          kk = sort(divisible_list[[kkk]], decreasing = TRUE)
          num_vertex_temp = length(kk)
          mm_temp = matrix(0,num_vertex_temp,num_vertex_temp)
          for(ii in 1:num_vertex_temp){
            for(jj in 1:num_vertex_temp){
              mm_temp[ii,jj] = mm[which(V(j)$name == kk[ii]),which(V(j)$name == kk[jj])]
            }
          }
          # 0 * 0 matrix
          if(!nrow(mm_temp)) {
            # divisible_list = divisible_list[-kkk]
            remove_list = c(remove_list,kkk)
            next
          }
          # eigen divide
          eg <- eigen(mm_temp)
          eigen_values = eg$values
          eigen_vectors = eg$vectors
          
          if(max(eigen_values) > 0) {
            lead_vector = eigen_vectors[,which(eigen_values == max(eigen_values))]
            group1 = c()
            group2 = c()
            for(ii in 1:num_vertex_temp){
              if(lead_vector[ii] < 0) group1 = c(group1, kk[ii])
              else group2 = c(group2, kk[ii])
            }
            #calculate s
            md2_s = replicate(num_vertex_temp, 0)
            md2_s[which(kk %in% group1)] = 1
            md2_s[which(kk %in% group2)] = -1
            
            # compute deltaQ
            dQ = 0
            for(ii in 1:num_vertex_temp){
              for(jj in 1:num_vertex_temp){
                dQ = dQ + (md2_s[ii] * md2_s[jj] - 1) * mm_temp[ii,jj]
              }
            }
            dQ = dQ / (4.0 * num_edges)
            dQs[[length(dQs) + 1]] = dQ
          }
          
          
          
          
          
          # remove from temp list
          # divisible_list = divisible_list[-kkk]
          remove_list = c(remove_list,kkk)
          if(dQ > 0 & length(group1) > 0) {
            # #remove from res list
            # md2_res = md2_res[unlist(md2_res) != unlist(kk)]
            # md2_res = discard(md2_res, ~ all(is.null(.x)))
            # add to temp list
            if(!is.null(group1)) divisible_list[[length(divisible_list) + 1]] = group1
            if(!is.null(group2)) divisible_list[[length(divisible_list) + 1]] = group2
            
            
            
          } else{
            # add to res list
            if(!is.null(group1))  md2_res[[length(md2_res) + 1]] = group1
            if(!is.null(group2))  md2_res[[length(md2_res) + 1]] = group2
          }
        }
        divisible_list = divisible_list[-remove_list]
      }
      md2_subgraph = list()
      for(i in 1:length(md2_res)){
        if (length(md2_res[[i]]) > 3) {
          md2_subgraph[[length(md2_subgraph) + 1]] = induced_subgraph(human_data_cpd,md2_res[[i]])
        }
      }
      
      for(i in md2_subgraph){
        
        # shave
        excessive = names(which(degree(i,setdiff(V(i)$name,seed_cpds)) == 1))
        while(length(excessive)) {
          i = delete_vertices(i, excessive)
          excessive = names(which(degree(i,setdiff(V(i)$name,seed_cpds)) == 1))
        }
        
        # calculate activity score
        # get empirical cpds number
        # just use EID
        # Ns
        emp_cpds_num = nrow(unique(ref_trio %>% filter(compounds %in% intersect(V(i)$name,seed_cpds)) %>% select(EID)))
        #Nm
        nodes_num = as.double(length(V(i)))
        Q = 0
        if(nodes_num > 0){
          # m 
          ref_edge_num = length(E(human_data_cpd))
          expected = 0
          for (ii in V(i)$name){
            for (jj in V(i)$name){
              if (ii != jj) {
                expected = expected + as.numeric(degree(human_data_cpd,ii)) * as.numeric(degree(human_data_cpd,jj))
              } 
            }
          }
          expected = expected / (4.0 * ref_edge_num)
          Q = (length(E(i)) - expected) / ref_edge_num
          as = sqrt(length(seed_cpds) / nodes_num) * Q * (emp_cpds_num / nodes_num) * 100
        }else{
          as = 0
          Q = 0
        }
        
        # filter duplicates
        exist_in_module = FALSE
        if(nrow(modules) > 0){
          for (md_row in 1:nrow(modules)){
            if ((length(intersect(V(i)$name, V(modules$subgraph[[md_row]])$name)) == length(V(i)$name)) & (length(V(i)$name) == length(V(modules$subgraph[[md_row]])$name)))
              exist_in_module = TRUE
          }
        }
        if(!exist_in_module)
          modules <- modules %>% add_row(subgraph = list(i), activity_score = as)
        
      }
      
    }
  }
  
  # filter th modules
  if(nrow(modules) > 0) {
    temp = c()
    for(i in 1:nrow(modules)){
      if((length(V(modules$subgraph[[i]])) <= 3) |(modules$activity_score[i] == 0)){
        temp <- c(temp,i)
      }
    }
    if (!is.null(temp))
      modules <- modules[-temp,]
  }
  
  # add the activity scores
  if(nrow(modules) > 0) {
    permuation_mscores = c(permuation_mscores, modules$activity_score)
  }
  temp_record <- c(temp_record, nrow(modules))
  
  
}


# rank significance
sig_modules <- sig_modules %>% mutate(p_value = NA)
for(i in 1:nrow(sig_modules)){
  # calculate p value
  total_scores <- c(sig_modules$activity_score[[i]])
  total_scores <- c(total_scores,permuation_mscores)
  total_scores <- sort(total_scores, decreasing = TRUE)
  
  D = length(permuation_mscores) + 1
  
  sig_modules$p_value[[i]] = (which(total_scores == sig_modules$activity_score[[i]])) / D
  
  
}

top_modules <- sig_modules %>% filter(p_value < SIGNIFICANCE_CUTOFF) %>% arrange(desc(p_value))




# ==============================================================================================
# Pathway analysis
q_set  = unique(trio$EID)
query_set_size = length(q_set)
total_feature_num = length(unique(ecompounds_v2$EID))


meta_pathway = data.frame(path_id=character(),
                                 name=character())
meta_pathway <- tibble(meta_pathway)
pathway_data = human.metabolic$metabolic_pathways

for(i in pathway_data){
  meta_pathway <- meta_pathway %>% add_row(path_id = i$id, name = i$name)
}
meta_pathway <- meta_pathway %>% add_column(compounds = NA, reactions = NA, enzymes = NA)
for(i in 1:length(pathway_data)){
  meta_pathway$compounds[[i]] <- list(pathway_data[[i]]$cpds)
  meta_pathway$reactions[[i]] <- list(pathway_data[[i]]$rxns)
  meta_pathway$enzymes[[i]] <- list(pathway_data[[i]]$ecs)
}
meta_pathway <- meta_pathway %>% add_column(ecompounds = NA)
for(i in 1:length(pathway_data)){
  meta_pathway$ecompounds[[i]] <- list(unique(ecompounds_v2[which(ecompounds_v2$compounds %in% unlist(meta_pathway$compounds[i])),]$EID))
}
meta_pathway <- meta_pathway %>% add_column(p_FET = NA)
for(i in 1:nrow(meta_pathway)){
  overlap_eompounds = intersect(unlist(meta_pathway$ecompounds[[i]]), q_set)
  overlap_size = length(overlap_eompounds)
  emp_size = length(meta_pathway$ecompounds[[i]][[1]])
  if(overlap_size > 0){
    neg_neg = total_feature_num + overlap_size - emp_size - query_set_size
    fisher_input <- matrix(c(overlap_size, query_set_size - overlap_size, emp_size - overlap_size, neg_neg), nrow = 2)
    meta_pathway$p_FET[[i]] = (fisher.test(fisher_input,alternative = "greater"))$p.value
  } else {
    meta_pathway$p_FET[[i]] = 1
  }
}


# do permutations
N = length(which(ListOfMassFeatures$is_significant==TRUE))
# permutation_pathway = tibble(data.frame(matrix(ncol = 0, nrow = 119)))
permutation_FET = c()
for (i in 1:NUM_PERM){
  ran_sample = sample(ref_rows,N)
  ref_trio <- ecompounds_v2 %>% filter(row %in% intersect(ecompounds_v2$row, ran_sample))
  query_ecompounds = unique(ref_trio$EID)
  per_query_set_size = length(query_ecompounds)
  
  for (j in 1:nrow(meta_pathway)){
    per_overlap_eompounds = intersect(unlist(meta_pathway$ecompounds[[j]]), query_ecompounds)
    per_overlap_size = length(per_overlap_eompounds)
    per_emp_size = length(meta_pathway$ecompounds[[j]][[1]])
    
    if(per_overlap_size > 0){
      per_neg_neg = total_feature_num + per_overlap_size - per_emp_size - per_query_set_size
      per_fisher_input <- matrix(c(per_overlap_size, per_query_set_size - per_overlap_size, per_emp_size - per_overlap_size, per_neg_neg), nrow = 2)
      permutation_FET = c(permutation_FET,(fisher.test(per_fisher_input,alternative = "greater"))$p.value)
    } else {
      permutation_FET = c(permutation_FET, 1)
    }
    
  }
  # permutation_pathway <- permutation_pathway %>% add_column(temp_p_FET,.name_repair = "minimal")
  
}

# calculate adjusted P
meta_pathway <- meta_pathway %>% add_column(adjusted_p = NA)
for(i in 1:nrow(meta_pathway)){
  # calculate p value
  total_scores <- c(meta_pathway$p_FET[[i]])
  total_scores <- c(total_scores,permutation_FET)
  total_scores <- sort(total_scores)
  
  D = length(permutation_FET) + 1
  
  meta_pathway$adjusted_p[[i]] = (which(total_scores == meta_pathway$p_FET[[i]]))[1] / D
  
  
}

# pathway collect_hit_Trios
overlap_sig_ecompounds <- c()
sig_meta_pathway <- meta_pathway %>% filter(adjusted_p < SIGNIFICANCE_CUTOFF)
for(i in 1:nrow(sig_meta_pathway)){
  overlap_sig_ecompounds <- union(sig_meta_pathway$ecompounds[[i]][[1]], overlap_sig_ecompounds)
}
# for(i in 1:nrow(trio)){
#   if((trio$EID[[i]] %in% overlap_sig_ecompounds) & (trio$row[[i]] %in% sig_rows)){
#     
#   }
# }
pathway_trio = trio %>% filter((EID %in% overlap_sig_ecompounds) & (trio$row %in% sig_rows))

# module collect_hit_Trios
overlap_sig_ecompounds2 <- c()
for(i in 1:nrow(top_modules)){
  overlap_sig_ecompounds2 <- c(overlap_sig_ecompounds2, V(top_modules$subgraph[[i]])$name)
}
overlap_sig_ecompounds2 <- unique(overlap_sig_ecompounds2)
module_trio = trio %>% filter(compounds %in% overlap_sig_ecompounds2)


new_trio = unique(module_trio %>% add_row(pathway_trio))


if (!is.null(new_trio$compounds)) {
  # find the largest subgraph
  index = 1
  max_sub_nodes_num = 0
  for(i in decompose.graph(induced_subgraph(human_data_cpd,unique(new_trio$compounds)))){
    if(length(V(i)) > max_sub_nodes_num) {
      index = i
      max_sub_nodes_num = length(V(i))
    }
  }
  sub1 = index
}



# clusters(sig_graph)


plot(human_data_cpd, edge.arrow.size=.5, vertex.label.color="black", vertex.label.dist=1.5, vertex.size = 0.5,vertex.label=NA)

network_print = human_data_cpd
highlight.these=V(top_modules$subgraph[[2]])$name
vertex.attributes(network_print)$color=ifelse(vertex.attributes(network_print)$name%in%highlight.these,"red","black")
vertex.attributes(network_print)$size=ifelse(vertex.attributes(network_print)$name%in%highlight.these,5,0.5)
plot(network_print, edge.arrow.size=.5, vertex.label.dist=1.5,vertex.label=NA)


plot(sig_graph, edge.arrow.size=.5, vertex.label.color="black", vertex.label.dist=1.5, vertex.size = 5,vertex.label=NA)

sub1_label = c()
for(i in V(sub1)$name){
  print(i)
  sub_name = human.metabolic$Compounds[[i]]$name
  print(sub_name)
  sub1_label = c(sub1_label,sub_name)
}
sub1_statistics = c()
for(i in V(sub1)$name){
  sub_row = new_trio[new_trio$compounds == i,]$row
  sub1_statistic = ListOfMassFeatures[ListOfMassFeatures$raw_row_number == sub_row[1],]$t.score
  sub1_statistics = c(sub1_statistics,sub1_statistic)
}
plot(sub1, edge.arrow.size=.5, vertex.label.color="black", vertex.label.dist=1.5, vertex.size = 5,vertex.label=sub1_label)

sub2 = top_modules$subgraph[[2]]
sub2_label = c()
for(i in V(sub2)$name){
  print(i)
  sub_name = human.metabolic$Compounds[[i]]$name
  print(sub_name)
  sub_name = str_split(sub_name, ";")[[1]][1]
  sub2_label = c(sub2_label,sub_name)
}
sub2_statistics = c()
for(i in V(sub2)$name){
  sub_row = new_trio[new_trio$compounds == i,]$row
  sub2_statistic = ListOfMassFeatures[ListOfMassFeatures$raw_row_number == sub_row[1],]$t.score
  sub2_statistics = c(sub2_statistics,sub2_statistic)
}
V(sub2)$label = sub2_label
V(sub2)$statistics = ifelse(sub2_statistics > 0, log10(sub2_statistics), -log10(-sub2_statistics))
write_graph(sub2, "sub2.graphml", format = "graphml")


# export seeds
human_data_cpd_export = human_data_cpd
human_mark1 = c()
human_mark2 = c()
human_mark3 = c()
human_mark4 = c()
for(i in V(human_data_cpd_export)$name){
  if(i %in% unique(seeds1)){
    human_mark1 = c(human_mark1, 1)
  } else {
    human_mark1 = c(human_mark1, 0)
  }
  
}
V(human_data_cpd_export)$mark1 = human_mark1
for(i in V(human_data_cpd_export)$name){
  if(i %in% unique(seed2)){
    human_mark2 = c(human_mark2, 1)
  } else {
    human_mark2 = c(human_mark2, 0)
  }
  
}
V(human_data_cpd_export)$mark2 = human_mark2
for(i in V(human_data_cpd_export)$name){
  if(i %in% unique(seeds3)){
    human_mark3 = c(human_mark3, 1)
  } else {
    human_mark3 = c(human_mark3, 0)
  }
  
}
V(human_data_cpd_export)$mark3 = human_mark3
for(i in V(human_data_cpd_export)$name){
  if(i %in% unique(seeds4)){
    human_mark4 = c(human_mark4, 1)
  } else {
    human_mark4 = c(human_mark4, 0)
  }
  
}
V(human_data_cpd_export)$mark4 = human_mark4
write_graph(human_data_cpd_export, "human_seeds.graphml", format = "graphml")

# barplot of adducts of compounds
adducts_num = sapply(human.metabolic$Compounds, function(x){length(x$adducts)})
adducts_num = sort(adducts_num)
adducts_num_df = tibble(name = names(adducts_num),adducts_num = adducts_num)
adducts_num_sum = adducts_num_df %>% group_by(adducts_num) %>% summarise(adducts_num_counts = length(adducts_num))
ggplot(adducts_num_sum, aes(x = adducts_num, y = adducts_num_counts)) + geom_bar(stat = "identity")+xlab("# adducts of compounds") + ylab("# cmopounds")
# barplot of adducts of m/z
mz_num = sapply(ion_cpd_tree, length)
mz_num_df = tibble(mz_num)
mz_num_df = mz_num_df %>% mutate(mz = 1:2000)
ggplot(mz_num_df, aes(x = mz, y = mz_num)) + geom_bar(stat = "identity")+xlab("mz") + ylab("# adducts")


# print activity network 
an = induced_subgraph(human_data_cpd,unique(new_trio$compounds))
module_cpds = module_trio$compounds
pathway_cpds =pathway_trio$compounds
intersect_cpds = intersect(pathway_trio$compounds,module_trio$compounds)
vertex.attributes(an)$color=ifelse(vertex.attributes(an)$name%in%intersect_cpds,"red",ifelse(vertex.attributes(an)$name%in%module_cpds, "blue","yellow"))
plot(an, edge.arrow.size=.5, vertex.label.dist=1.5,vertex.size = 5,vertex.label=NA)








