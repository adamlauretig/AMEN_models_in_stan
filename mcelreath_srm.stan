//mcelreath srm

data{
    int N_households;
    int N; // number of dyads
    int giftsAB[300]; // Directed: a -> b
    int giftsBA[300]; // Directed: b -> a
    int did[300]; // dyad indicator
    int hidA[300];
    int hidB[300];
}
parameters{
    real a;
    vector[2] gr[N_households]; // giving & receiving effects per household
    corr_matrix[2] Rho_gr; // correlation matrix w/in hh
    vector<lower=0>[2] sigma_gr; // sd w/in household
    matrix[2,N] z; // for non-centered parameterization
    cholesky_factor_corr[2] L_Rho_d; // correlation matrix for dyads
    real<lower=0> sigma_d; // dyad sd
}
transformed parameters{
    matrix[N,2] d; 
    d = (diag_pre_multiply(rep_vector(sigma_d, 2), L_Rho_d) * z)'; // sd *corredation
}
model{
    vector[300] lambdaAB; // rate param for poisson
    vector[300] lambdaBA; // rate param for poisson
    sigma_d ~ exponential( 1 ); //dyad sd
    L_Rho_d ~ lkj_corr_cholesky( 8 ); // correlation matrix for dyads
    to_vector( z ) ~ normal( 0 , 1 ); //none-centered parameterization
    sigma_gr ~ exponential( 1 ); 
    Rho_gr ~ lkj_corr( 4 );
    gr ~ multi_normal( rep_vector(0,2) , quad_form_diag(Rho_gr , sigma_gr) ); // sd * corr
    a ~ normal( 0 , 1 ); // intercept
    for ( i in 1:300 ) {
        lambdaBA[i] = a + gr[hidB[i], 1] + gr[hidA[i], 2] + d[did[i], 2]; // random effect
        lambdaBA[i] = exp(lambdaBA[i]); // link fn
    }
    for ( i in 1:300 ) {
        lambdaAB[i] = a + gr[hidA[i], 1] + gr[hidB[i], 2] + d[did[i], 1];
        lambdaAB[i] = exp(lambdaAB[i]);
    }
    giftsBA ~ poisson( lambdaBA );
    giftsAB ~ poisson( lambdaAB );
}
generated quantities{
    matrix[2,2] Rho_d;
    Rho_d = multiply_lower_tri_self_transpose(L_Rho_d);
}

