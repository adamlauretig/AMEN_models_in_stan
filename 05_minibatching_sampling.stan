// TODO: go through, figure out how to make batch indexing work for sending/receiving nodes and dyads.
data{
  int n_nodes ;
  int n_dyads ;
  int N; //total obs. should be n_dyads * 2
  int sender_id[n_dyads * 2] ;
  int receiver_id[n_dyads * 2] ;
  int K ; // number of latent dimensions
  real Y[n_dyads * 2] ;
  int batch_size ;
  int dyad_idx[n_dyads] ;
  int n_batches ;
}
transformed data{
  int lower_tri_size = (K * (K + 1))/2 ;
  simplex[n_dyads] uniform = rep_vector(1.0 / n_dyads, n_dyads);
  int<lower = 1, upper = n_dyads> batch_idxs[batch_size, n_batches];
  for (n in 1:n_batches)
    for (i in 1:batch_size)
      batch_idxs[i , n] = categorical_rng(uniform) ;
}

parameters{
  real intercept ;
  cholesky_factor_cov[2] cov_nodes; // correlation matrix w/in hh
  matrix[2, n_nodes] z_nodes ; // for node non-centered parameterization
  // vector<lower=-pi()/2, upper=pi()/2>[lower_tri_size] cov_raw_values;
  cholesky_factor_cov[K * 2] cov_multi_effects ; // correlation matrix for multiplicative effect
  matrix[K * 2, n_nodes] z_multi_effects; // Multi-effect non-centered term 
  
}
transformed parameters{
    matrix[n_nodes, 2] mean_nodes; // node parameter mean
    matrix[n_nodes, K * 2] mean_multi_effects ; // multi-effect mean
    // vector [lower_tri_size] cov_values;
   
    // cov_values = 3 * tan(cov_raw_values);
    
    mean_nodes = (cov_nodes * z_nodes)' ;
    mean_multi_effects = (cov_multi_effects * z_multi_effects)'; // sd *correlation  
    
}
model{
  intercept ~ normal(0, 5) ;
  
  //node terms
  to_vector(z_nodes) ~ normal(0, 1) ;
  to_vector(cov_nodes) ~ student_t(3, 0, 2) ;

  // multi-effect terms
  to_vector(z_multi_effects) ~ normal(0, 1) ;
  to_vector(cov_multi_effects) ~ student_t(3, 0, 2) ;

  for(batch in 1:n_batches){
    for(obs in 1:batch_size){
      
      
    Y[batch_idxs[obs, batch]] ~ normal(intercept + 
          mean_nodes[sender_id[n], 1] + mean_nodes[receiver_id[n], 2] + 
          mean_multi_effects[sender_id[n], 1:K] * 
          (mean_multi_effects[receiver_id[n], (K+1):(K*2)])', 
          1 );
    }
  }
}
generated quantities{
  real Y_sim[N] ;
  for(n in 1:N){
    Y_sim[n] = normal_rng(intercept + 
        mean_nodes[sender_id[n], 1] + mean_nodes[receiver_id[n], 2] + 
        mean_multi_effects[sender_id[n], 1:K] * 
        (mean_multi_effects[receiver_id[n], (K+1):(K*2)])', 1) ;
  }
  
}
