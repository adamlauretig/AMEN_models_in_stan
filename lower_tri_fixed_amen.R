# fitting a stan model with an AMEN model with an identified lower triangle

library(amen)
library(rstan)
library(data.table)
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

m4_code <- stan_model(file =  "04_stan_fixed_lower_tri.stan"
  )

m4 <- vb(m4_code, 
     data = list(
       N = data_for_stan$N,
       n_nodes = data_for_stan$n_nodes,
       n_dyads = data_for_stan$n_dyads,
       sender_id = data_for_stan$edgelist[, 1],
       receiver_id = data_for_stan$edgelist[, 2],
       K = 5,
       Y = data_for_stan$edgelist[, 3]), 
     # chains = 4, cores = 4, 
     iter = 10000, 
     seed = 123)

m4_params <- extract(m4)

preds4 <- apply(m4_params$Y_sim, 2, mean)
plot(data_for_stan$edgelist[, 3], preds4)
mean((data_for_stan$edgelist[, 3] - preds4)^2)


latent_params <- apply(m4_params$mean_multi_effects, c(2:3), mean)
latent_params_dt <- data.table(latent_params)
setnames(latent_params_dt, c(paste0("sender_", 1:5), paste0("receiver_", 1:5)))
latent_params_dt$country <- data_for_stan$lookup_table[, 1]
latent_params_dt <- melt(latent_params_dt, id.vars = "country")
latent_params_dt[, c("sr", "dim_num") := tstrsplit(variable, "_")]
latent_params_dt <- dcast(latent_params_dt, country + dim_num ~ sr)
ggplot(data = latent_params_dt, aes(x = sender, y = receiver)) + 
  geom_text(aes(label = country)) + 
  facet_wrap(~dim_num)
ggplot(data = latent_params_dt, aes(x = sender, y = receiver)) + 
  geom_text(aes(label = country)) + 
  facet_wrap(~dim_num) + 
  scale_x_continuous(limits = c(-4, 4)) + scale_y_continuous(limits = c(-4, 4))
ggplot(data = latent_params_dt, aes(x = sender, y = receiver)) + 
  geom_text(aes(label = country)) + 
  facet_wrap(~dim_num) + 
  scale_x_continuous(limits = c(-1, 1)) + scale_y_continuous(limits = c(-1, 1))
