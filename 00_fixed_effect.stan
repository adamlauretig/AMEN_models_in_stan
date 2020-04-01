
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
  vector[n_nodes] sender_beta; 
  vector[n_nodes] receiver_beta; 
}

model{
  intercept ~ normal(0, 5) ;
  sender_beta ~ normal(0, 5) ;
  receiver_beta ~ normal(0, 5) ;

  for(n in 1:N){
  Y[n] ~ normal(intercept + 
        sender_beta[sender_id[n]] + receiver_beta[receiver_id[n]], 1 );
  }
}
generated quantities{
  real Y_sim[N] ;
  for(n in 1:N){
    Y_sim[n] = normal_rng(intercept + 
        sender_beta[sender_id[n]] + receiver_beta[receiver_id[n]], 1) ;
  }
  
}
