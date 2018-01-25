data {

  int L;	// number of raw water treatment plants (Locations)
  int N;	// number of observations
  int M; 	// number of samples	
  int K;	// number of sources

  // covariates
  int<lower=0, upper=1> DM[N]; 	// DM selected background model 
  int<lower=0, upper=1> Local[N]; // Local selected background model

  // indicator variables
  int<lower=1, upper=L> lev3ForLev2[N];	// Location factor, level 3 in hierarchy
  int<lower=1, upper=M> level2[N];	// Sample per location, level 2 in hierarchy
  int<lower=1, upper=K> level2S[N];	// Sources indicator 

  // offset
  real offset[N];	// log total number of OTUs behind the signal

  // Count outcome
  int<lower=0> y[N];	// the number of OTUs with posterior probability > 0.7

}
parameters {
  // Define parameters to estimate

  real a;	// mean across sites
  real beta[2];
  real<lower=0.001> sigma_g; // std dev of random intercept for site
  real<lower=0.001> sigma_t; // std dev of random intercept for site x location
  real<lower=0.001> sigma_s; // std dev of random intercept for source
  real<lower=0.001> sigma_u; // std dev of random intercept for sample
  // residuals
  real u_j[L]; 	// site specific residual
  real u_jk[M];	// site x location residual
  real u_s[K]; 	// source specific residual

}
model {

  vector[L] a_j; 	// group coefficient for intercept of each site
  vector[M] a_jk; // group coefficient for intercept of each sample x site
  vector[K] b_s; 	// group coefficient for intercept of each site
  real mu[N];

  // priors
  sigma_g ~ cauchy(0, 5); // prior on the standard deviations
  sigma_t ~ cauchy(0, 5); 
  sigma_u ~ cauchy(0, 5); 
  sigma_s ~ cauchy(0, 5); 
  // Level-3 site level (L level-3 random intercepts)
  for (j in 1:L) {
    a_j[j] = a + u_j[j];
  }
  // Level-2 (M level-2 random intercepts)
  for (k in 1:M) {
    a_jk[k] = a_j[lev3ForLev2[k]] + u_jk[k];
  }
  // Level-2 source level 
  for (j in 1:K) {
    b_s[j] = u_s[j]; 
  }

  for (j in 1:L)
    u_j[j] ~ normal(0,sigma_g); // residual for location
  for (k in 1:M)
    u_jk[k] ~ normal(0,sigma_t); // residual for sample x location
  for (j in 1:K)
    u_s[j] ~ normal(0,sigma_s); // residual for source
  // Individual mean at the replicate level
  for (i in 1:N) {
     // Linear predictor
     mu[i] = a_jk[level2[i]] + b_s[level2S[i]] + beta[1] * DM[i] + beta[2] * Local[i] + offset[i]; 
  }
  // likelihood
  y ~ poisson_log(mu);

}

