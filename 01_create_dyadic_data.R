library(amen)
library(rstan)
data(IR90s)
str(IR90s$dyadvars)


# A function to convert a matrix of pairs to an edgelist
# Assumes that rows == columns, but will check
sociomatrix_to_convert <- IR90s$dyadvars[, , 2]

matrix_to_edgelist <- function(sociomatrix_to_convert){
  
  if(nrow(sociomatrix_to_convert) != ncol(sociomatrix_to_convert)){
    stop("nrows != ncols")
  }
  
  all_nodes <- expand.grid(
    list(rownames(sociomatrix_to_convert), colnames(sociomatrix_to_convert)))
  all_nodes[, 1] <- as.numeric(all_nodes[, 1])
  all_nodes[, 2] <- as.numeric(all_nodes[, 2])
  all_nodes <- all_nodes[all_nodes[, 1] != all_nodes[, 2], ]

  lookup_table <- data.frame(
    node_names = rownames(sociomatrix_to_convert), 
    idx = as.numeric(factor(rownames(sociomatrix_to_convert))))
  obs_values <- which(sociomatrix_to_convert != 0, arr.ind = TRUE)
  obs_values <- cbind(obs_values, sociomatrix_to_convert[obs_values])
 
   edgelist <- merge(
    all_nodes, obs_values, by.x = c("Var1", "Var2"), 
    by.y = c("row", "col"), all.x = TRUE)
  
  edgelist$V3 <- ifelse(is.na(edgelist$V3), 0, edgelist$V3)
  edgelist <- cbind(
    edgelist, 
    rep(seq.int(from = 1, to = sum(edgelist[, 1] < edgelist[, 2])), 2))
  
  return(list(edgelist = edgelist, lookup_table = lookup_table, 
       n_nodes = max(edgelist[, 1]), n_dyads = max(edgelist[, 4]), N = nrow(edgelist)))
  
}

