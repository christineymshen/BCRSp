// adapted from crsm_mgp6.stan, updated everything to be consistent with CRS working files including dimensions of parameters also to allow for different bandwidths
functions {
	// competing risks cox proportional hazards model, log likelihood
  vector crcox_loglik(int n, int m, matrix Xbeta, matrix delta, matrix lambda_cal, matrix H) {
    matrix[m,n] loglik_mat = delta .* (log(lambda_cal) + Xbeta) - exp(Xbeta) .* H;
    vector[n] loglik;
    for (i in 1:n){
      loglik[i] = sum(loglik_mat[,i]);
    }
    return loglik;
  }
}
data {
  int<lower=1> n; // number of data points
  int<lower=1> p; // number of covariates
  int<lower=1> k; // number of bins for piecewise linear baseline hazard rates
  int<lower=1> m; // number of risk types
  vector[n] y; // time to event data
  vector[k+1] s; // knots for the piecewise linear baseline hazard rates
  matrix[n,m] delta; // indicator of risk type
  matrix[n,p] X; // design matrix
  real<lower=0> a0; // hyperparameter for psi0
  real<lower=0> b0; // hyperparameter for psi0
  real<lower=0> a1; // hyperparameter for kappa
  real<lower=0> b1; // hyperparameter for kappa
}
transformed data {
  array[n] int<lower=1,upper=k> ki; // bin for each event
  matrix[m,k] s_diff; // bin widths, matrix form for calculation
  matrix[m,n] T; // time since last knot, matrix form for calculation
	vector[m] zeros_m = zeros_vector(m);
  for (i in 1:n) {
    for (r in 1:k){
      if (y[i]>s[r] && y[i] <= s[r+1]){
        ki[i] = r;
        T[,i] = rep_vector(y[i] - s[r],m);
      }
    }
  }
  for (r in 1:k){
    s_diff[,r] = rep_vector(s[r+1]-s[r],m);
  }
  vector[m*k] s_diff_inv_1d = to_vector(inv(s_diff));
}
parameters {
  matrix[p,m] beta;
  vector<lower=0>[m] psi0;
  matrix<lower=0>[m,k] psi;
  real<lower=0> s2;
  real<lower=0> kappa; // for multiplicative gamma process
}
transformed parameters {
  matrix[m,k] lambda; // baseline hazard rates
  lambda[,1] = psi0 .* psi[,1]; // psi0 corresponds to (imaginary) baseline hazard rate at time 0
  for (r in 2:k){
    lambda[,r] = psi[,r] .* lambda[,r-1]; // setup in this way because matrix is column major ordered in stan
  }
  vector[n] log_lik;
  {
    matrix[m,n] H; // baseline cumulative hazard rates for each ti
    matrix[n,m] Xbeta = X * beta;
    matrix[m,n] lambda_cal = lambda[,ki]; // lambda for each ti
    matrix[m,k] cum_lambda_tmp = lambda .* s_diff;
    matrix[m,k] cum_lambda;
    cum_lambda[,1] = zeros_m;
    // do it in this way because matrix in stan follows column major order, should be more efficient than using cumulative_sum for each row vector?
    for (r in 2:k){
      cum_lambda[,r] = cum_lambda[,r-1] + cum_lambda_tmp[,r-1];
    }
    H = cum_lambda[,ki] + lambda_cal .* T;
    log_lik = crcox_loglik(n,m,Xbeta', delta', lambda_cal, H);
  }
}
model {
  vector[k*m] psi_1d = to_vector(psi);
  vector[p*m] beta_1d = to_vector(beta);
  target += normal_lupdf(beta_1d | 0, sqrt(s2));
  target += inv_gamma_lupdf(s2 |1,1);
  target += gamma_lupdf(psi0 | a0, b0);
  target += gamma_lupdf(psi_1d | kappa*s_diff_inv_1d, kappa*s_diff_inv_1d);
  target += inv_gamma_lupdf(kappa | a1, b1);
  target += sum(log_lik);
}

