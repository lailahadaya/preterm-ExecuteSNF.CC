# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# Bootstrap code - to run on an HPC cluster (we used the King's College London Rosalind computing interface: https://rosalind.kcl.ac.uk/) 

# Study title: "Characterising distinct subgroups of very preterm born children and 
#   exploring differences in neonatal structural and functional brain patterns" 
# Date: 19/11/2021
# Authors: L Hadaya, K Dimitrakopoulou, L Vanes, D Kanel, S Fenn-Moltu, D Pecheva, G Ball, AD Edwards, 
# SJ Counsell, M Saqi, D Batalle and C Nosarti
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =  

### ==== Install Packages ====
#install.packages("SNFtool")
#install.packages("readxl")
#install.packages("openxlsx")
#install.packages("BiocManager")
#BiocManager::install("ConsensusClusterPlus")
#install.packages("cluster")
#BiocManager::install("CancerSubtypes")
#install.packages("purrr")
#install.packages("ggplot2")
#install.packages("gdata")
#install.packages("diceR")

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

### ==== 1. Set things up for Rosalind parallisation =====

args = commandArgs(trailingOnly=TRUE) # This reads input arguments 
paste0('Seed #', as.numeric(args[1])) # This converts the input to numeric and set as seed
seednumber <- as.numeric(args[1])
set.seed(seednumber)

#df_data_type_1 <- df_ex_func 
#bootstrap_80 <- purrr::rerun(1, resample(rownames(df_data_type_1), size=0.8*(dim(df_data_type_1)[1])))
#bootstrap_80 <- unlist(bootstrap_80)

load("bootstrap_80.RData")
bootstrap_80 <- unlist(bootstrap_80[[seednumber]])

### ==== 2. Load your data types from your RData.env for the specific run of interest =====
## !!!! NOTE !!!!!! 
## make sure data type 3 is the one with categorical data - this data type will have gowers distance in the code below

load("data_types.RData")

df_data_type_1 <-  df_ex_func
df_data_type_2 <- df_socio_emo
df_data_type_3 <- as.data.frame(Filter(is.data.frame, mget(ls(pattern="risk$"))))
print(ls(Filter(is.data.frame, mget(ls(pattern="risk$")))))

### ==== 3. Set param combos =====

N= dim(df_data_type_1)[1]
N #number of subjects
parameter_combos_K <- c(rep(10,6), rep(15,6), rep(round((N/10), digits=0),6),rep(25,6), rep(30,6))
parameter_combos_alpha <- c(rep(c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8),5))

param_combos <- data.frame(c(1:length(parameter_combos_alpha)), parameter_combos_K, parameter_combos_alpha)
colnames(param_combos) <- c("combo_Num", "K", "alpha")

### ==== 4. Load functions with plot = NULL =====
# set plot = "png" instead of NULL, if you get an error message saying:
# "Error in .External.graphics(C_layout, num.rows, num.cols, mat, as.integer(num.figures),  : 
#  invalid graphics state

ExecuteCC2_noPlot <- function (clusterNum, d, maxK = 10, clusterAlg = "hc", distance = "pearson", 
                               title = "ConsensusClusterResult", reps = 1000, pItem = 0.8, # set plot = "png" instead of NULL, if you get an error message
                               pFeature = 1, plot = NULL, innerLinkage = "average", finalLinkage = "average", 
                               writeTable = FALSE, weightsItem = NULL, weightsFeature = NULL, 
                               verbose = FALSE, corUse = "everything") 
{
  if (is.list(d)) {
    temp = NULL
    for (i in 1:length(d)) {
      temp = rbind(temp, d[[i]])
    }
    temp = t(scale(t(temp)))
  }
  else temp = d
  originalResult = ConsensusClusterPlus(temp, maxK = maxK, 
                                        clusterAlg = clusterAlg, distance = distance, title = title, 
                                        reps = reps, pItem = pItem, pFeature = pFeature, plot = NULL, # set plot = "png" instead of NULL, if you get an error message
                                        innerLinkage = innerLinkage, finalLinkage = finalLinkage, 
                                        writeTable = writeTable, weightsItem = weightsItem, weightsFeature = weightsFeature, 
                                        verbose = verbose, corUse = corUse, seed=123456)
  group = originalResult[[clusterNum]][["consensusClass"]]
  distanceMatrix = originalResult[[clusterNum]][["consensusMatrix"]]
  attr(distanceMatrix, "class") = "Similarity"
  result = list(group = group, distanceMatrix = distanceMatrix, 
                originalResult = originalResult)
  result
}

ExecuteSNF.CC_new_noPlot <- function (W_fused, clusterNum, 
                                      maxK = 10, pItem = 0.8, reps = 1000, title = "ConsensusClusterResult", 
                                      plot = NULL, finalLinkage = "average") 
{
  
  W_fused = as.dist(W_fused)
  result = ExecuteCC2_noPlot(clusterNum = clusterNum, d = W_fused, maxK = maxK, 
                             clusterAlg = "spectralAlg", title = title, reps = reps, 
                             pItem = pItem, plot = NULL, finalLinkage = finalLinkage) # set plot = "png" instead of NULL, if you get an error message
  result
}


### ==== 5. Bootstrap function ====

bootstrap_SNFcc_silh_hyperparam <- function(df_data_type_1, df_data_type_2, df_data_type_3, K, alpha, t=20) {
  group_SNF_bootstrap_c2 <- list()
  group_SNF_bootstrap_c3 <- list()
  group_SNF_bootstrap_c4 <- list()
  mean_silh_score_c2 <- list()
  mean_silh_score_c3 <- list()
  mean_silh_score_c4 <- list()
  for (j in 1:dim(param_combos)[1]) { 
    # 1. normalise data frames
    N_df_data_type_1 <- standardNormalization(df_data_type_1[bootstrap_80,]) # the bootstrap_80[[i]] has a random resampled 80% of the data 
    N_df_data_type_2 <- standardNormalization(df_data_type_2[bootstrap_80,]) # the bootstrap_80[[i]] has a random resampled 80% of the data 
    # 2. get distance matrices
    sq_dist_N_df_data_type_1 <- (dist2(as.matrix(N_df_data_type_1), as.matrix(N_df_data_type_1)))^(1/2)
    sq_dist_N_df_data_type_2 <- (dist2(as.matrix(N_df_data_type_2), as.matrix(N_df_data_type_2)))^(1/2)
    dist_gower_df_data_type_3 <- daisy(df_data_type_3[bootstrap_80,], metric = c("gower"))  
    # 3. get similarity matrices (w)
    w.data_type_1 <- affinityMatrix(sq_dist_N_df_data_type_1, K=param_combos[j,2], sigma=param_combos[j,3])
    w.data_type_2 <- affinityMatrix(sq_dist_N_df_data_type_2, K=param_combos[j,2], sigma=param_combos[j,3])
    w.data_type_3 <- affinityMatrix(as.matrix(dist_gower_df_data_type_3), K=param_combos[j,2], sigma=param_combos[j,3])
    # 4. get fused W
    W_SNF <- SNF(list(w.data_type_1, w.data_type_2, w.data_type_3), K=param_combos[j,2], t=20)
    # 5. get SNF_CC clusters
    SNF_CC_c2 <- ExecuteSNF.CC_new_noPlot(W_SNF, clusterNum=2, maxK=5, pItem=0.8, reps=1000) 
    SNF_CC_c3 <- ExecuteSNF.CC_new_noPlot(W_SNF, clusterNum=3, maxK=5, pItem=0.8, reps=1000)
    SNF_CC_c4 <- ExecuteSNF.CC_new_noPlot(W_SNF, clusterNum=4, maxK=5, pItem=0.8, reps=1000)
    # 6. get SNF_CC groupings from each iteration  
    # c2
    group_SNF_bootstrap_c2[[j]] <- cbind(bootstrap_80, SNF_CC_c2$group) 
    colnames(group_SNF_bootstrap_c2[[j]]) <- c("ID", "c2")
    # c3
    group_SNF_bootstrap_c3[[j]] <- cbind(bootstrap_80, SNF_CC_c3$group) 
    colnames(group_SNF_bootstrap_c3[[j]]) <- c("ID", "c3")
    # c4
    group_SNF_bootstrap_c4[[j]] <- cbind(bootstrap_80, SNF_CC_c4$group) 
    colnames(group_SNF_bootstrap_c4[[j]]) <- c("ID", "c4")
    # 7. get silhouette scores from each iteration
    # 7a. sil score for c2
    dist_matrix_c2 <- SNF_CC_c2$distanceMatrix
    silh_scores_from_SNF_bootstrap_c2 <- silhouette_SimilarityMatrix(group = SNF_CC_c2$group, similarity_matrix=dist_matrix_c2)
    mean_silh_score_c2[[j]] <- mean(silh_scores_from_SNF_bootstrap_c2[,3])
    # 7b. sil score for c3
    dist_matrix_c3 <- SNF_CC_c3$distanceMatrix
    silh_scores_from_SNF_bootstrap_c3 <- silhouette_SimilarityMatrix(group = SNF_CC_c3$group, similarity_matrix=dist_matrix_c3)
    mean_silh_score_c3[[j]] <- mean(silh_scores_from_SNF_bootstrap_c3[,3])
    # 7b. sil score for c4
    dist_matrix_c4 <- SNF_CC_c4$distanceMatrix
    silh_scores_from_SNF_bootstrap_c4 <- silhouette_SimilarityMatrix(group = SNF_CC_c4$group, similarity_matrix=dist_matrix_c4)
    mean_silh_score_c4[[j]] <- mean(silh_scores_from_SNF_bootstrap_c4[,3])
  }
  # }
  # 8. Append all outputs to one list 
  print(c(group_SNF_bootstrap_c2, group_SNF_bootstrap_c3, group_SNF_bootstrap_c4, mean_silh_score_c2, mean_silh_score_c3, mean_silh_score_c4)) # to print out all lists, appended
  output_list <- list("group_c2"=group_SNF_bootstrap_c2, "group_c3"=group_SNF_bootstrap_c3, "group_c4"=group_SNF_bootstrap_c4, "silh_c2" = mean_silh_score_c2, "silh_c3"=mean_silh_score_c3, "silh_c4"=mean_silh_score_c4)
  return(output_list)
}
### ==== 6. Run the bootstrap ====
# make sure data type 3 is the one with categorical data - this data type will have gowers dist in the code below

bootstrap_groups_and_silhs <- bootstrap_SNFcc_silh_hyperparam(df_data_type_1, df_data_type_2, df_data_type_3, K, alpha, t=20)


### ==== 7. Extract OUTPUT: bootstrap_groups and mean_silh scores ====
# Extract from the main output (bootstrap_groups_and_silhs), SNF groupings and sil scores: 

list_bootstrap_groups_c2 <- bootstrap_groups_and_silhs$group_c2
list_bootstrap_groups_c3 <- bootstrap_groups_and_silhs$group_c3
list_bootstrap_groups_c4 <- bootstrap_groups_and_silhs$group_c4

# turn list into a data.frame - and we are merging by ID number
#===== for c=2
multi_full_c2 <- Reduce(
  function(x, y, ...) merge(x, y, by="ID", all = TRUE, ...),
  list_bootstrap_groups_c2
) 

#===== for c=3
multi_full_c3 <- Reduce(
  function(x, y, ...) merge(x, y, by="ID", all = TRUE, ...),
  list_bootstrap_groups_c3
)  

#===== for c=4
multi_full_c4 <- Reduce(
  function(x, y, ...) merge(x, y, by="ID", all = TRUE, ...),
  list_bootstrap_groups_c4
)  

write.xlsx(multi_full_c2, file=paste("run_", seednumber, "param_combo_multi_full_BIG_bootstrap_c2.xlsx", sep =""))
write.xlsx(multi_full_c3, file=paste("run_", seednumber, "param_combo_multi_full_BIG_bootstrap_c3.xlsx", sep =""))
write.xlsx(multi_full_c4, file=paste("run_", seednumber, "param_combo_multi_full_BIG_bootstrap_c4.xlsx", sep =""))

### ==== 8. Extract OUTPUT: mean_silh scores ====
# Extract from the main output (bootstrap_groups_and_silhs): sil scores

#===== for c=2
list_mean_silhs_c2 <- bootstrap_groups_and_silhs$silh_c2
final_mean_silh_c2 <- unlist(list_mean_silhs_c2)
write.xlsx(final_mean_silh_c2, file=paste("run_", seednumber, "param_combo_mean_silh_scores_c2.xlsx", sep =""))

#===== for c=3
list_mean_silhs_c3 <- bootstrap_groups_and_silhs$silh_c3
final_mean_silh_c3 <- unlist(list_mean_silhs_c3)
write.xlsx(final_mean_silh_c3, file=paste("run_", seednumber, "param_combo_mean_silh_scores_c3.xlsx", sep =""))

#===== for c=4
list_mean_silhs_c4 <- bootstrap_groups_and_silhs$silh_c4
final_mean_silh_c4 <- unlist(list_mean_silhs_c4)
write.xlsx(final_mean_silh_c4, file=paste("run_", seednumber, "param_combo_mean_silh_scores_c4.xlsx", sep =""))
