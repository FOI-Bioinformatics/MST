data {

  int L;	// number of raw water treatment plants (rwtp)
  int N;	// number of samples	
  int K;	// number of measured outcomes

  real yn[N];	// proporions obtained without background water
  real yc[N];	// proporions obtained using a DM-background
  real yl[N];	// proporions obtained using a local background

  int<lower=0> yi[N]; // FIB 1 counts
  int<lower=0> yj[N]; // FIB 2 counts
  int<lower=0> yk[N]; // FIB 3 counts	 
  int x[N]; 	// design matrix relating at which location each sample is taken 
 
}
parameters {
     
  vector[K] mu_a;	// overall intercept
  corr_matrix[K] Omega; 
  vector<lower=0>[K] sigma; 
  vector[K] aj[N];    	// outcome specific means, K for each of N samples
  vector[K] a[L]; 	// group coefficient for intercept of each rwtp
  real<lower=0.001> sigma_g[K]; // std dev of random intercept
  real<lower=0.001> M_sd[3]; // std dev for normal distribution of contaminations

}
transformed parameters {
  cov_matrix[K] Sigma; 
  Sigma = quad_form_diag(Omega, sigma); 
}
model {

  // priors
  for (j in 1:L)
    a[j] ~ normal(mu_a,sigma_g); // group coefficient for rwtp
  sigma ~ cauchy(0, 5); // prior on the standard deviations
  M_sd ~ cauchy(0,5);
  Omega ~ lkj_corr(1); // LKJ prior on the correlation matrix   
  for ( j in 1:N ) aj[j] ~ multi_normal( rep_vector(0,K) , Sigma );

  // likelihoods
  for ( j in 1:N ) {
    yn[j] ~ normal(a[x[j],1] + aj[j,1], M_sd[1]);
    yc[j] ~ normal(a[x[j],2] + aj[j,2], M_sd[2]);
    yl[j] ~ normal(a[x[j],3] + aj[j,3], M_sd[3]);
    yi[j] ~ poisson_log(a[x[j],4] + aj[j,4] );
    yj[j] ~ poisson_log(a[x[j],5] + aj[j,5] );
    yk[j] ~ poisson_log(a[x[j],6] + aj[j,6] );
  }

}

