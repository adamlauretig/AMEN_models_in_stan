
data{
  int n_nodes ;
  int n_dyads ;
  int N; //total obs. should be n_dyads * 2
  int sender_id[n_dyads * 2] ;
  int receiver_id[n_dyads * 2] ;
  int sender_idx[n_dyads * 2] ;
  int receiver_idx[n_dyads * 2] ;
//  int dyad_idx[n_dyads * 2] ;
  real Y[n_dyads * 2] ;
  
}
parameters{
  real intercept ;
  cholesky_factor_corr[2] corr_nodes; // correlation matrix w/in hh
  vector<lower=0>[2] sigma_nodes; // sd w/in household
  matrix[2, n_nodes] z_nodes ; // for node non-centered parameterization
  // cholesky_factor_corr[2] corr_dyads; // correlation matrix w/in hh
  // real<lower=0> sigma_dyads; // sd w/in household
  // matrix[2, n_dyads] z_dyads ; // for node non-centered parameterization
}
transformed parameters{
    // matrix[n_dyads,2] mean_dyads; 
    matrix[n_nodes, 2] mean_nodes; 
    real linear_predictor[N] ;
    
  for(n in 1:N){
    linear_predictor[n] = mean_nodes[sender_idx[n], 1] + mean_nodes[receiver_idx[n], 2] ;
  }
    
    // mean_dyads = (diag_pre_multiply(
    //   rep_vector(sigma_dyads, 2), corr_dyads) * z_dyads)'; // sd *corredation
    mean_nodes = (diag_pre_multiply(
      sigma_nodes, corr_nodes) * z_nodes)'; // sd *corredation
    
}
model{
  intercept ~ normal(0, 5) ;
  to_vector(z_nodes) ~ normal(0, 1) ;
  // to_vector(z_dyads) ~ normal(0, 1) ;
  corr_nodes ~ lkj_corr_cholesky(5) ;
  // corr_dyads ~ lkj_corr_cholesky(5) ;
  sigma_nodes ~ gamma(1, 1) ;
  // sigma_dyads ~ gamma(1, 1) ;
  Y ~ normal(linear_predictor, 1) ;
  
}

