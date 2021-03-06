library(amen)
library(rstan)
data(IR90s)


# A function to convert a matrix of pairs to an edgelist
# Assumes that rows == columns, but will check

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

data_for_stan <- matrix_to_edgelist(IR90s$dyadvars[, , 2]) # trade data




m0_code <- stan_model(file =  "00_fixed_effect.stan"
  )

m0 <- vb(m0_code, 
     data = list(
       N = data_for_stan$N,
       n_nodes = data_for_stan$n_nodes,
       n_dyads = data_for_stan$n_dyads,
       sender_id = data_for_stan$edgelist[, 1],
       receiver_id = data_for_stan$edgelist[, 2],
       Y = data_for_stan$edgelist[, 3]), 
     seed = 123,
     # chains = 4, cores = 4,
     iter = 10000 
     )

m0_params <- extract(m0)
preds0 <- apply(m0_params$Y_sim, 2, mean)
plot(data_for_stan$edgelist[, 3], preds0)

m1_code <- stan_model(file =  "01_srm_stan.stan"
  )

m1 <- vb(m1_code, 
     data = list(
       N = data_for_stan$N,
       n_nodes = data_for_stan$n_nodes,
       n_dyads = data_for_stan$n_dyads,
       sender_id = data_for_stan$edgelist[, 1],
       receiver_id = data_for_stan$edgelist[, 2],
       Y = data_for_stan$edgelist[, 3]), 
     seed = 123,
     # chains = 4, cores = 4,
     iter = 10000 
     )

m1_params <- extract(m1)
preds <- apply(m1_params$Y_sim, 2, mean)
plot(data_for_stan$edgelist[, 3], preds)

m2_code <- stan_model(file =  "02_srm_stan_dyad.stan"
  )

m2 <- vb(m2_code, 
     data = list(
       N = data_for_stan$N,
       n_nodes = data_for_stan$n_nodes,
       n_dyads = data_for_stan$n_dyads,
       sender_id = data_for_stan$edgelist[, 1],
       receiver_id = data_for_stan$edgelist[, 2],
       dyad_id = data_for_stan$edgelist[, 4],
       send_receive = data_for_stan$edgelist[, 5],
       Y = data_for_stan$edgelist[, 3]), 
     seed = 123,
     # chains = 4, cores = 4,
     iter = 10000 
     )


m2_params <- extract(m2)
preds2 <- apply(m2_params$Y_sim, 2, mean)
plot(data_for_stan$edgelist[, 3], preds2)


m3_code <- stan_model(file =  "03_amen_stan.stan"
  )

m3 <- vb(m3_code, 
     data = list(
       N = data_for_stan$N,
       n_nodes = data_for_stan$n_nodes,
       n_dyads = data_for_stan$n_dyads,
       sender_id = data_for_stan$edgelist[, 1],
       receiver_id = data_for_stan$edgelist[, 2],
       K = 10,
       Y = data_for_stan$edgelist[, 3]), 
     # chains = 4, cores = 4, 
     iter = 10000, 
     seed = 123)

m3_params <- extract(m3)
preds3 <- apply(m3_params$Y_sim, 2, mean)
plot(data_for_stan$edgelist[, 3], preds3)

mean((data_for_stan$edgelist[, 3] - preds0)^2)
mean((data_for_stan$edgelist[, 3] - preds)^2)
mean((data_for_stan$edgelist[, 3] - preds2)^2)

mean((data_for_stan$edgelist[, 3] - preds3)^2)
save(m0, m1, m2, m3, file = "srm_amen_stan.rdata")
