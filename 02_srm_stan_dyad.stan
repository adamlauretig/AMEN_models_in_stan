
data{
  int n_nodes ;
  int n_dyads ;
  int N; //total obs. should be n_dyads * 2
  int sender_id[n_dyads * 2] ;
  int receiver_id[n_dyads * 2] ;
  int dyad_id[n_dyads * 2] ;
  int send_receive[n_dyads * 2] ;
  real Y[n_dyads * 2] ;
  
}
parameters{
  real intercept ;
  cholesky_factor_corr[2] corr_nodes; // correlation matrix w/in hh
  vector<lower=0>[2] sigma_nodes; // sd w/in household
  matrix[2, n_nodes] z_nodes ; // for node non-centered parameterization
  cholesky_factor_corr[2] corr_dyads; // correlation matrix w/in hh
  real<lower=0> sigma_dyads; // sd w/in household
  matrix[2, n_dyads] z_dyads ; // for node non-centered parameterization
}
transformed parameters{
    matrix[n_dyads,2] mean_dyads;
    matrix[n_nodes, 2] mean_nodes; 
    // real linear_predictor[N] ;
    
    mean_dyads = (diag_pre_multiply(
      rep_vector(sigma_dyads, 2), corr_dyads) * z_dyads)'; // sd *correlation
    mean_nodes = (diag_pre_multiply(
      sigma_nodes, corr_nodes) * z_nodes)'; // sd *correlation
    
    
    // for(n in 1:N){
    //   linear_predictor[n] = 
    // }
}
model{
  intercept ~ normal(0, 5) ;
  to_vector(z_nodes) ~ normal(0, 1) ;
  // to_vector(z_dyads) ~ normal(0, 1) ;
  corr_nodes ~ lkj_corr_cholesky(5) ;
  // corr_dyads ~ lkj_corr_cholesky(5) ;
  sigma_nodes ~ gamma(1, 1) ;
  // sigma_dyads ~ gamma(1, 1) ;
  // Y ~ normal(linear_predictor, 1) ;
  
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
