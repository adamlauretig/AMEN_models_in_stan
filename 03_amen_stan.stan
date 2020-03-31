
data{
  int n_nodes ;
  int n_dyads ;
  int N; //total obs. should be n_dyads * 2
  int sender_id[n_dyads * 2] ;
  int receiver_id[n_dyads * 2] ;
  int K ; // number of latent dimensions
  real Y[n_dyads * 2] ;
  
}
parameters{
  real intercept ;
  cholesky_factor_corr[2] corr_nodes; // correlation matrix w/in hh
  vector<lower=0>[2] sigma_nodes; // sd w/in nodes
  matrix[2, n_nodes] z_nodes ; // for node non-centered parameterization
  cholesky_factor_corr[K * 2] corr_multi_effects ; // correlation matrix for multiplicative effect
  vector<lower=0>[K * 2] sigma_multi_effects ; // sd
  matrix[K * 2, n_nodes] z_multi_effects; // Multi-effect non-centered term 
  
}
transformed parameters{
    matrix[n_nodes, 2] mean_nodes; // node parameter mean
    matrix[n_nodes, K * 2] mean_multi_effects ; // multi-effect mean
    
    mean_nodes = (diag_pre_multiply(
      sigma_nodes, corr_nodes) * z_nodes)'; // sd *correlation
    mean_nodes = (diag_pre_multiply(
      sigma_multi_effects, corr_multi_effects) * z_multi_effects)'; // sd *correlation  
    
}
model{
  intercept ~ normal(0, 5) ;
  to_vector(z_nodes) ~ normal(0, 1) ;
  corr_nodes ~ lkj_corr_cholesky(5) ;
  sigma_nodes ~ gamma(1, 1) ;

  for(n in 1:N){
  Y[n] ~ normal(intercept + 
        mean_nodes[sender_id[n], 1] + mean_nodes[receiver_id[n], 2], 1 );
  }
}
generated quantities{
  real Y_sim[N] ;
  for(n in 1:N){
    Y_sim[n] = normal_rng(intercept + 
        mean_nodes[sender_id[n], 1] + mean_nodes[receiver_id[n], 2], 1) ;
  }
  
}
