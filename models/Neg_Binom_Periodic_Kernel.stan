data {
  // Observed Datapoints and the Unobserved Points At Which We Want to Predict
  int<lower=1> N_train;             // Number of observations (i.e. datapoints)
  int<lower=1> N_pred;              // Number of points to be predicted (i.e. not associated observation)
  int y[N_train];                   // Actual observation, output value
  real x_train[N_train];            // Actual observation, index position
  real x_pred[N_pred];              // Point to be predicted, index position
  
  // Model Parameters
	real<lower=0> length_scale_mean;  // Length scale prior mean and sd
	real<lower=0> length_scale_sd;    
	real<lower=0> period_mean;        // Period prior mean and sd
	real<lower=0> period_sd;          
	real<lower=0> alpha_mean;         // Alpha prior mean and sd (controlling variability of GP) 
  real<lower=0> alpha_sd;  
  real<lower=0> overdispersion_mean;  // Negative binomial overdispersion prior mean and sd 
  real<lower=0> overdispersion_sd;  
}

transformed data {
  
  int<lower=1> N_tot;           // Length of vector containing both actual datapoints and those to be predicted
  int<lower=1> k;               // A Counter
  real x_tot[N_pred + N_train]; // Vector storing the index positions for both actual datapoints and those to be predicted
  
  N_tot = N_pred + N_train;     
  k = 1;
  
  for (n in 1:N_train) {        // Filling the vector x_tot with index positions for 1) the observed datapoints and 2) the points we're predicting at
    x_tot[k] = x_train[n];      // Filling with the index of observed datapoints first
    k = k + 1;
  }
  for (n in 1:N_pred) {   
    x_tot[k] = x_pred[n];       // Then filling with the index of points we want to predict at 
    k = k + 1;
  }
}

parameters {
  real<lower=0.5,upper=12> length_scale;
  real<lower=4,upper=18> period;
  real<lower=0> alpha;
  vector[N_tot] eta;
  real<lower=0> overdispersion;
}

transformed parameters {
  vector[N_tot] f;        
  {                                 // Think this stops this stuff from being stored by STAN? Unsure. 
    matrix[N_tot, N_tot] L;
    matrix[N_tot, N_tot] Sigma;     // Initialising and filling the covariance matrix, with values determined by the periodic kernel below
	  for (i in 1:N_tot) {   
  		for (j in 1:i) {
  			if (i == j) 
  			  Sigma[i,j] = pow(alpha, 2) * exp(-(2/pow(length_scale, 2)) * pow((sin((pi() * (x_tot[i] - x_tot[j]))/period)), 2));
  			else 
  			  Sigma[i,j] = pow(alpha, 2) * exp(-(2/pow(length_scale, 2)) * pow((sin((pi() * (x_tot[i] - x_tot[j]))/period)), 2));
    		  Sigma[j,i] = Sigma[i,j];
  		}
		  Sigma[i,i] = Sigma[i,i] + 1e-12;
	  }
    L = cholesky_decompose(Sigma);  // cholesky decomposition of the covariance matrix
    f = L * eta;                    // choleksy decompose * normal random variate is equivalent to sampling from the required MVN
  }
}

// Model Block 
model{
	  // Priors 
	  length_scale ~ normal(length_scale_mean, length_scale_sd);
	  period ~ normal(period_mean, period_sd);
	  alpha ~ normal(alpha_mean, alpha_sd);
    overdispersion ~ normal(overdispersion_mean, overdispersion_sd);
    
    // Sample from Normal for MVN Sampling
	  eta ~ normal(0, 1);
	  
	  // Calculating the Likelihood
	  y ~ neg_binomial_2(exp(f[1:N_train]), overdispersion); 
}
