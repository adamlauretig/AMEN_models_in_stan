
data{
  int n_nodes ;
  int n_dyads ;
  int N; //total obs. should be n_dyads * 2
  int sender_id[n_dyads * 2] ;
  int receiver_id[n_dyads * 2] ;
  real Y[n_dyads * 2] ;
  
}
parameters{
  real intercept ;
  cholesky_factor_corr[2] corr_nodes; // correlation matrix 
  vector<lower=0>[2] sigma_nodes; // sd 
  matrix[2, n_nodes] z_nodes ; // for node non-centered parameterization
}
transformed parameters{
    matrix[n_nodes, 2] mean_nodes; 
    mean_nodes = (diag_pre_multiply(
      sigma_nodes, corr_nodes) * z_nodes)'; // non-centered parameterization
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
