# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# Main pipeline.

# Study title: "Characterising distinct subgroups of very preterm born children and 
#   exploring differences in neonatal structural and functional brain patterns" 
# Date: 19/11/2021
# Authors: L Hadaya, K Dimitrakopoulou, L Vanes, D Kanel, S Fenn-Moltu, D Pecheva, G Ball, AD Edwards, 
# SJ Counsell, M Saqi, D Batalle and C Nosarti
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =  

### ==== Install Packages ====
install.packages("SNFtool")
install.packages("readxl")
install.packages("openxlsx")
install.packages("BiocManager")
BiocManager::install("ConsensusClusterPlus")
install.packages("cluster")
BiocManager::install("CancerSubtypes")
install.packages("purrr")
install.packages("ggplot2")
install.packages("gdata")
install.packages("diceR")

### ==== Load libraries ====
library(SNFtool)
library(readxl)
library(openxlsx)
library(ConsensusClusterPlus)
library(cluster)
library(CancerSubtypes)
library(data.table)
library(purrr)
library(gdata)
library(ggplot2)
library(diceR)
library(alluvial)
library(dplyr)
library(ggplot2)

### ==== I. Load .RData for selected data types + sample ====

load("data_types.RData") # load your data
df_data_type_1 <-  df_ex_func # first data type - all numeric variables
df_data_type_2 <- df_socio_emo # second data type - all numeric variables
df_data_type_3 <- df_clin_social_risk # third data type - mixed data type

### ==== II. Cluster number estimation table creation ====
#For each combination of SNF hyperparameters K={10, 15, 20, 25, 30} & alpha={0.3,, 0.4, 0.5, 0.6, 0.7, 0.8} 
#Calculate best and second best C based on eigengap and rotation cost using, estimateNumberOfClustersGivenGraph()
#Choose the cluster numbers C that most combinations agree on and are clinically relevant.

# Set param combinations
N = dim(df_data_type_1)[1]
N #number of subjects
parameter_combos_K <- c(rep(10,6), rep(15,6), rep(round((N/10), digits=0),6),rep(25,6), rep(30,6))
parameter_combos_alpha <- c(rep(c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8),5))

param_combos <- data.frame(c(1:length(parameter_combos_alpha)), parameter_combos_K, parameter_combos_alpha)
colnames(param_combos) <- c("combo_Num", "K", "alpha")

# Get Cluster Estimations for all 30 combination runs:
estClustNumb_paramCombos <- data.frame()
for(i in 1:dim(param_combos)[1]){
  paste("K =", param_combos[i,2], ", alpha =", param_combos[i,3])
  #SNF pipeline:
  N_df_data_type_1 <- standardNormalization(df_data_type_1)
  N_df_data_type_2 <- standardNormalization(df_data_type_2)
  sq_dist_N_df_data_type_1 <- (dist2(as.matrix(N_df_data_type_1), as.matrix(N_df_data_type_1)))^(1/2)
  sq_dist_N_df_data_type_2 <- (dist2(as.matrix(N_df_data_type_2), as.matrix(N_df_data_type_2)))^(1/2)
  dist_gower_df_data_type_3 <- daisy(df_data_type_3, metric = c("gower")) 
  w.data_type_1 <- affinityMatrix(sq_dist_N_df_data_type_1, K=param_combos[i,2], sigma=param_combos[i,3]) # sigma is what we also refer to as alpha
  w.data_type_2 <- affinityMatrix(sq_dist_N_df_data_type_2, K=param_combos[i,2], sigma=param_combos[i,3]) # sigma is what we also refer to as alpha
  w.data_type_3 <- affinityMatrix(as.matrix(dist_gower_df_data_type_3), K=param_combos[i,2], sigma=param_combos[i,3]) # sigma is what we also refer to as alpha
  W_SNF <- SNF(list(w.data_type_1, w.data_type_2, w.data_type_3), K=param_combos[i,2], t=20)
  estClustNumb <- print(estimateNumberOfClustersGivenGraph(W_SNF, 2:10))
  estClustNumb_paramCombos[1,i] <- estClustNumb$`Eigen-gap best`
  estClustNumb_paramCombos[2,i] <- estClustNumb$`Eigen-gap 2nd best`
  estClustNumb_paramCombos[3,i] <- estClustNumb$`Rotation cost best`
  estClustNumb_paramCombos[4,i] <- estClustNumb$`Rotation cost 2nd best`
  rownames(estClustNumb_paramCombos) <- c("Eigengap-best", "Eigengap-2nd-best", "rotation-cost-best", "rotation-cost-2nd-best")
  print(estClustNumb_paramCombos)
}

for(i in 1:dim(param_combos)[1]){
  print(paste("K =", param_combos[i,2], ", alpha =", param_combos[i,3]))
}

### ==== III. Run SNF bootstrap code ====
# === 1. Create 1000 randomly resampled IDs with 80% of subjects 
df_data_type_1 <- df_ex_func
bootstrap_80 <- purrr::rerun(1000, resample(rownames(df_data_type_1), size=0.8*(dim(df_data_type_1)[1])))
save(bootstrap_80, file="bootstrap_80.RData") # copy this R object over to Rosalind

# === 2. Run bootstrap - code for this is in "Hadaya_bootstrap_code.R" - this was run on Rosalind.

#choose the cluster numbers based the cluster numbers that came up most in section (II.)

### ==== IV. MERGE the 1000 outputs that the bootstrap code creates ====
# copy the outputs from Rosalind to local folder

### ==== i. LOAD and MERGE the iterations ===
### Groupings: LOAD
list_ALL_bootstraps_c2 =list_ALL_bootstraps_c3 = list() # list with results from the 10 or 100 or X iterations
for (i in 1:1000) {
  list_ALL_bootstraps_c2[[i]] <- as.data.frame(read_excel(paste("run_", i, "param_combo_multi_full_BIG_bootstrap_c2.xlsx", sep ="")))
  list_ALL_bootstraps_c3[[i]] <- as.data.frame(read_excel(paste("run_", i, "param_combo_multi_full_BIG_bootstrap_c3.xlsx", sep ="")))
  list_ALL_bootstraps_c4[[i]] <- as.data.frame(read_excel(paste("run_", i, "param_combo_multi_full_BIG_bootstrap_c4.xlsx", sep ="")))
}

### Groupings: MERGE
#===== for c=2
multi_full_c2 <- Reduce(
  function(x, y, ...) merge(x, y, by="ID", all = TRUE, ...),
  list_ALL_bootstraps_c2  
) #this is creating a full/outer join between all elements in the list


#===== for c=3
multi_full_c3 <- Reduce(
  function(x, y, ...) merge(x, y, by="ID", all = TRUE, ...),
  list_ALL_bootstraps_c3 
) #this is creating a full/outer join between all elements in the list

#===== for c=4
multi_full_c4 <- Reduce(
  function(x, y, ...) merge(x, y, by="ID", all = TRUE, ...),
  list_ALL_bootstraps_c4
) #this is creating a full/outer join between all elements in the list


### Silh Scores: 
# LOAD 
list_ALL_sil_c2 = list_ALL_sil_c3 = list_ALL_sil_c4 = list() # list with results from the 1000 iterations
for (i in 1:1000) {
  list_ALL_sil_c2[[i]] <- as.data.frame(read_excel(paste("run_", i, "param_combo_mean_silh_scores_c2.xlsx", sep =""), col_names=FALSE))
  list_ALL_sil_c3[[i]] <- as.data.frame(read_excel(paste("run_", i, "param_combo_mean_silh_scores_c3.xlsx", sep =""), col_names=FALSE))
  list_ALL_sil_c4[[i]] <- as.data.frame(read_excel(paste("run_", i, "param_combo_mean_silh_scores_c4.xlsx", sep =""), col_names=FALSE))
  
  }

# MERGE: 
multi_mean_silh_c2 <- unlist(list_ALL_sil_c2)
multi_mean_silh_c3 <- unlist(list_ALL_sil_c3)
multi_mean_silh_c4 <- unlist(list_ALL_sil_c4)

### ==== V. SAVE the merged iterations file ====

write.xlsx(multi_full_c2, file = "multi_full_c2.xlsx")
write.xlsx(multi_full_c3, file = "multi_full_c3.xlsx")
write.xlsx(multi_full_c4, file = "multi_full_c4.xlsx")

write.xlsx(multi_mean_silh_c2, file = "multi_mean_silh_c2.xlsx")
write.xlsx(multi_mean_silh_c3, file = "multi_mean_silh_c3.xlsx")
write.xlsx(multi_mean_silh_c4, file = "multi_mean_silh_c4.xlsx")

### ==== VI. Filter based on silhouette scores: keep highest scoring grouping, from each iteration, only ====
#data.frame(1st column: number of run, then: 30 rows with the combos, repeated 1000 times, so end up with 30000 rows, 2nd column: index subsample)
full_silh_data_QC <- data.frame("run_number"=c(1:30000), "seednumber"=c(sort(rep(1:1000,30))), "combo_numb"=rep(1:30, 1000), "K"= rep(param_combos[,2], 1000), "alpha"=rep(param_combos[,3], 1000), multi_mean_silh_c2, multi_mean_silh_c3, multi_mean_silh_c4)

# 1. look at keeping one top silh per seednumber (end up with 1000)
keep_highest_only_full_silh_data_QC_c2 <- full_silh_data_QC %>%
  group_by(seednumber) %>%
  slice(which.max(multi_mean_silh_c2))

keep_highest_only_full_silh_data_QC_c3 <- full_silh_data_QC %>%
  group_by(seednumber) %>%
  slice(which.max(multi_mean_silh_c3))

keep_highest_only_full_silh_data_QC_c4 <- full_silh_data_QC %>%
  group_by(seednumber) %>%
  slice(which.max(multi_mean_silh_c4))

# 2. filtered columns to feed into the consensus_combine algorithms
c2_GROUPS_highest_silh_only <- multi_full_c2[,c(1, keep_highest_only_full_silh_data_QC_c2$run_number+1)]
c3_GROUPS_highest_silh_only <- multi_full_c3[,c(1, keep_highest_only_full_silh_data_QC_c3$run_number+1)]
c4_GROUPS_highest_silh_only <- multi_full_c4[,c(1, keep_highest_only_full_silh_data_QC_c4$run_number+1)]

### ==== VII. Apply Consensus Combine and use for stats =======
rownames(c2_GROUPS_highest_silh_only) <- c2_GROUPS_highest_silh_only$ID
rownames(c3_GROUPS_highest_silh_only) <- c3_GROUPS_highest_silh_only$ID
rownames(c4_GROUPS_highest_silh_only) <- c4_GROUPS_highest_silh_only$ID

#c2
sub1.c2<- as.data.frame(sapply(c2_GROUPS_highest_silh_only[,-1], as.integer))
sub2.c2<-array(as.matrix(as.data.frame(sub1.c2)),dim=c(198,1000,1,1),dimnames=list(rownames(sub1.c2),colnames(sub1.c2),"SNF","2")) # N = 198
cc_c2 <- as.numeric(unlist(consensus_combine(sub2.c2, element = "class")))
cc_matrix_c2 <- matrix(unlist(consensus_combine(sub2.c2, element = "matrix")), nrow=198, ncol=198) # N = 198

#c3
sub1.c3 <- as.data.frame(sapply(c3_GROUPS_highest_silh_only[,-1], as.integer))
sub2.c3<-array(as.matrix(as.data.frame(sub1.c3)),dim=c(198,1000,1,1),dimnames=list(rownames(sub1.c3),colnames(sub1.c3),"SNF","3")) # N = 198
cc_c3 <- as.numeric(unlist(consensus_combine(sub2.c3, element = "class")))
cc_matrix_c3 <- matrix(unlist(consensus_combine(sub2.c3, element = "matrix")), nrow=198, ncol=198) # N = 198

#c4
sub1.c4 <- as.data.frame(sapply(c4_GROUPS_highest_silh_only[,-1], as.integer))
sub2.c4 <-array(as.matrix(as.data.frame(sub1.c4)),dim=c(198,1000,1,1),dimnames=list(rownames(sub1.c4),colnames(sub1.c4),"SNF","4")) # N = 198
cc_c4 <- as.numeric(unlist(consensus_combine(sub2.c4, element = "class")))
cc_matrix_c4 <- matrix(unlist(consensus_combine(sub2.c4, element = "matrix")), nrow=198, ncol=198) # N = 198

### ==== VIII. Get Silhouette scores: =======

#c2
distmatrix=cc_matrix_c2
sil_CC_c2 <- silhouette_SimilarityMatrix(group = cc_c2 ,similarity_matrix=distmatrix)
plot(silhouette_SimilarityMatrix(group = cc_c2 ,similarity_matrix=distmatrix))
plot(sil_CC_c2, do.n.k=F, do.clus.stat=F, col=c("darkolivegreen2","blanchedalmond"))       

#c3
### --- Here, relabelled cc_c3: 2 as 3 and 3 as 2, because in the manuscript: 
### --- 1 = resilient group, 2=risk subgroup and 3 - is the NEW emerging (intermediate subgroup)
cc_c3_relabelled 
cc_c3_relabelled <- as.numeric(cc_c3_relabelled)
distmatrix <- cc_matrix_c3
sil_CC_c3_relabelled <- silhouette_SimilarityMatrix(group = cc_c3_relabelled ,similarity_matrix=distmatrix)
plot(sil_CC_c3_relabelled, do.n.k=F, do.clus.stat=F, col=c("darkolivegreen2","blanchedalmond", "darksalmon"))       

#c4
cc_c4_relabelled <- as.numeric(cc_c4_relabelled)
distmatrix <- cc_matrix_c4
sil_CC_c4_relabelled <- silhouette_SimilarityMatrix(group = cc_c4_relabelled ,similarity_matrix=distmatrix)
plot(sil_CC_c4_relabelled, do.n.k=F, do.clus.stat=F, col=c("darkolivegreen2","blanchedalmond", "darksalmon", "cyan4"))       
