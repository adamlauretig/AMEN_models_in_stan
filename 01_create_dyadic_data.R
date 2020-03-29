library(amen)
library(rstan)
data(IR90s)
str(IR90s$dyadvars)


# A function to convert a matrix of pairs to an edgelist
# Assumes that rows == columns, but will check
sociomatrix_to_convert <- IR90s$dyadvars[, , 2] # trade data

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
  edgelist <- cbind(edgelist, ifelse(edgelist[, 1] < edgelist[, 2], 1, 2))
  colnames(edgelist) <- c("sender", "receiver", "y", "dyad_id", "sr_indicator")
  # for dyad list include a (1, 2) indicator for send/receive
  # if s < r -> 1, else 2
  
  return(list(edgelist = edgelist, lookup_table = lookup_table, 
       n_nodes = max(edgelist[, 1]), n_dyads = max(edgelist[, 4]), 
       N = nrow(edgelist)))
  
}

data_for_stan <- matrix_to_edgelist(sociomatrix_to_convert)

m1 <- stan(file = "01_srm_stan.stan", 
     data = list(
       N = data_for_stan$N,
       n_nodes = data_for_stan$n_nodes,
       n_dyads = data_for_stan$n_dyads,
       sender_id = data_for_stan$edgelist[, 1],
       receiver_id = data_for_stan$edgelist[, 2],
       Y = data_for_stan$edgelist[, 3]), chains = 4, iter = 2000, cores = 4)

m1_params <- extract(m1)
preds <- apply(m1_params$Y_sim, 2, mean)
plot(data_for_stan$edgelist[, 3], preds)
