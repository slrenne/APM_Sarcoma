
data{
    int DM[126];
    vector[126] months_to_event;
    int A[126];
    int D[126];
    int G[126];
}
parameters{
    matrix[8,2] z_a;
    vector<lower=0>[8] sigma_a;
    cholesky_factor_corr[8] L_Rho_a;
    matrix[8,3] z_b;
    vector<lower=0>[8] sigma_b;
    cholesky_factor_corr[8] L_Rho_b;
}
transformed parameters{
    matrix[2,8] alpha;
    alpha = (diag_pre_multiply(sigma_a, L_Rho_a) * z_a)';

    matrix[3,8] beta;
    beta = (diag_pre_multiply(sigma_b, L_Rho_b) * z_b)';
}
model{
    vector[126] lambda;
    vector[126] mu;
    
    L_Rho_a ~ lkj_corr_cholesky( 2 );
    sigma_a ~ exponential( 1 );
    to_vector( z_a ) ~ normal( 0 , 1 );
    
    L_Rho_b ~ lkj_corr_cholesky( 2 );
    sigma_b ~ exponential( 1 );
    to_vector( z_b ) ~ normal( 0 , 1 );
    
    for ( i in 1:126 ) {
        mu[i] = alpha[A[i], D[i]] + beta[G[i], D[i]] ;
        mu[i] = exp(mu[i]);
    }
    for ( i in 1:126 ) {
        lambda[i] = 1/mu[i];
    }
    for ( i in 1:126 ) 
        if ( DM[i] == 0 ) target += exponential_lccdf(months_to_event[i] | lambda[i]);
    for ( i in 1:126 ) 
        if ( DM[i] == 1 ) months_to_event[i] ~ exponential( lambda[i] );
}
generated quantities{
    matrix[8,8] Rho_a;
    Rho_a = multiply_lower_tri_self_transpose(L_Rho_a);   
    matrix[8,8] Rho_b;
    Rho_b = multiply_lower_tri_self_transpose(L_Rho_b);
}
