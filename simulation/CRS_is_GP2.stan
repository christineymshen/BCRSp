functions {
	matrix L_matern32(array[] vector d, real alpha, real l, real epsilon) { 
		matrix[size(d), size(d)] cov = add_diag(gp_matern32_cov(d,alpha,l), epsilon);
		return cholesky_decompose(cov);
	}
  vector theta_pred_rng(array[] vector d, array[] vector d_pred, real alpha, real l, 
                        matrix LK, vector theta_obs, real epsilon) {
		int n = size(d);
		int n_pred = size(d_pred);
    vector[n_pred] theta_pred;
		{
      matrix[n,n_pred] K12 = gp_matern32_cov(d,d_pred,alpha,l);
      matrix[n,n_pred] LKinv_K12 = mdivide_left_tri_low(LK,K12);
      vector[n] LKinv_theta_obs = mdivide_left_tri_low(LK,theta_obs);
      vector[n_pred] theta_pred_mu = LKinv_K12' * LKinv_theta_obs;
      matrix[n_pred,n_pred] K22 = add_diag(gp_matern32_cov(d_pred,alpha,l), epsilon);
      matrix[n_pred,n_pred] theta_pred_cov = K22 - LKinv_K12'*LKinv_K12;
      theta_pred = multi_normal_rng(theta_pred_mu, theta_pred_cov);
		}
		return theta_pred;
	}
	// competing risks cox proportional hazards model, log likelihood
  real crcox_lpdf(vector y, matrix Xbeta, matrix delta, matrix lambda_cal, matrix H) {
    return sum(delta .* (log(lambda_cal) + Xbeta) - exp(Xbeta) .* H);
  }
}
data {
  int<lower=1> n; // number of observed data points
  int<lower=1> n_pred; // total number of points for prediction
  int<lower=1> p; // number of fixed covariates
  int<lower=1> k; // number of bins for piecewise linear baseline hazard rates
  int<lower=1> m; // number of risk types
  vector[n] y; // time to event data
  vector[k+1] s; // knots for the piecewise linear baseline hazard rates
  matrix[n,m] delta; // indicator of risk type
  matrix[n,p] X; // design matrix - fixed effects excluding W
  vector[n] W; // design vector - for spatial coefficients
  array[n] vector[2] d; // locations for observed data
  array[n_pred] vector[2] d_pred; // locations for prediction
  real<lower=0> kappa0; // for multiplicative gamma process
  real<lower=0> a; // hyperparameter for kappa1
  real<lower=0> b; // hyperparameter for kappa1
}
transformed data {
  real epsilon = 0.000001;
	vector[m] zeros_m = zeros_vector(m);
  array[n] int<lower=1,upper=k> ki; // bin for each event
  matrix[m,k] s_diff; // bin widths, matrix form for calculation
  matrix[m,n] T; // time since last knot, matrix form for calculation
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
  vector[m] beta_w;
  vector<lower=0>[m] psi0;
  matrix<lower=0>[m,k] psi;
 	vector<lower=0>[m] alpha0; // GP scale for intercept
 	vector<lower=0>[m] alpha1; // GP scale for slope
	vector<lower=0>[m] l0; // GP length scale for intercept
	vector<lower=0>[m] l1; // GP length scale for slope
  matrix[n,m] z0; // standard normal
  matrix[n,m] z1; // standard normal
  real<lower=0> s2;
  real<lower=0> kappa1; // for multiplicative gamma process
}
transformed parameters {
  array[m] matrix[n,n] LK0;
  array[m] matrix[n,n] LK1;
  array[m] vector[n] theta0;
  array[m] vector[n] theta1;
  {
    for (j in 1:m){
      LK0[j] = L_matern32(d,alpha0[j],l0[j],epsilon);
      theta0[j] = LK0[j] * col(z0,j); // spatial intercept
      LK1[j] = L_matern32(d,alpha1[j],l1[j],epsilon);
      theta1[j] = LK1[j] * col(z1,j); // spatial intercept
    }    
  }
  matrix[m,k] lambda; // baseline hazard rates
  lambda[,1] = psi0 .* psi[,1];
  for (r in 2:k){
    lambda[,r] = psi[,r] .* lambda[,r-1];
  }
}
model {
  matrix[m,n] H; // baseline cumulative hazard rates for each ti
  matrix[n,m] Xbeta = X * beta;
  for (j in 1:m){
    Xbeta[,j] = Xbeta[,j] + theta0[j] + (theta1[j]+beta_w[j]).*W;
  }
  vector[k*m] psi_1d = to_vector(psi);
  vector[p*m] beta_1d = to_vector(beta);
  vector[m*n] z0_1d =to_vector(z0);
  vector[m*n] z1_1d =to_vector(z1);
  matrix[m,n] lambda_cal = lambda[,ki]; // lambda for each ti
  {
    matrix[m,k] cum_lambda1 = lambda .* s_diff;
    matrix[m,k] cum_lambda2;
    cum_lambda2[,1] = zeros_m;
    for (r in 2:k){
      cum_lambda2[,r] = col(cum_lambda2,r-1) + col(cum_lambda1,r-1);
    }
    H = cum_lambda2[,ki] + lambda_cal .* T;
  }
  target += normal_lupdf(beta_1d | 0, sqrt(s2));
  target += normal_lupdf(beta_w | 0, sqrt(s2));
  target += inv_gamma_lupdf(s2 |1,1);
  target += gamma_lupdf(psi0 | kappa0, kappa0);
  target += gamma_lupdf(psi_1d | kappa1*s_diff_inv_1d, kappa1*s_diff_inv_1d);
  target += inv_gamma_lupdf(kappa1 | a, b);
  target += inv_gamma_lupdf(l0 |2,1); // same prior as in demo
  target += normal_lupdf(alpha0 |0,4); // same prior as in demo
  target += inv_gamma_lupdf(l1 |2,1); // same prior as in demo
  target += normal_lupdf(alpha1 |0,4); // same prior as in demo
  target += normal_lupdf(z0_1d |0,1);
  target += normal_lupdf(z1_1d |0,1);
  target += crcox_lupdf(y | Xbeta, delta, lambda_cal', H');
}
generated quantities {
  array[m] matrix[n_pred,2] theta_pred;
  for (j in 1:m){
    theta_pred[j][,1] = theta_pred_rng(d,d_pred,alpha0[j],l0[j],LK0[j],theta0[j],epsilon);
    theta_pred[j][,2] = theta_pred_rng(d,d_pred,alpha1[j],l1[j],LK1[j],theta1[j],epsilon);      
  }
}

