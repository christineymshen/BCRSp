functions{
	// matern 32 kernel
	vector sqrt_diagSPD_2D(real alpha, real l, vector L, matrix indices) {
		return alpha * sqrt(6*pi())* 3^(3/4) * l * (3+(indices^2 * (l*pi() / (2*L))^2))^(-5/4);
  	}
	matrix PHI_2D(int N, int M1, int M2, real L1, real L2, vector x1, vector x2) {
		matrix[N,M1*M2] PHI;
		matrix[N,M1] PHI_1 = sin(diag_post_multiply(rep_matrix(pi()/(2*L1) * (x1+L1), M1), linspaced_vector(M1, 1, M1)))/sqrt(L1);
		matrix[N,M2] PHI_2 = sin(diag_post_multiply(rep_matrix(pi()/(2*L2) * (x2+L2), M2), linspaced_vector(M2, 1, M2)))/sqrt(L2);
		PHI[,1:M2] = rep_matrix(PHI_1[,1], M2) .* PHI_2;
		for(i in 2:M1)
		  PHI[,1:(M2*i)] = append_col(PHI[,1:(M2*(i-1))], rep_matrix(PHI_1[,i], M2) .* PHI_2);
		return PHI;
	}
	// competing risks cox proportional hazards model, log likelihood
  vector crcox_loglik(int n, int m, matrix Xbeta, matrix delta, matrix lambda_cal, matrix H) {
    matrix[m,n] loglik_mat = delta .* (log(lambda_cal) + Xbeta) - exp(Xbeta) .* H;
    vector[n] loglik;
    for (i in 1:n){
      loglik[i] = sum(loglik_mat[,i]);
    }
    return loglik;
  }
  vector theta_pred_rng(matrix PHI_sqrt_SPD, matrix PHI_pred_sqrt_SPD, vector z_obs){
    int n = rows(PHI_sqrt_SPD);
    int m = cols(PHI_sqrt_SPD);
		int n_pred = rows(PHI_pred_sqrt_SPD);
    matrix[m,n] Q = qr_thin_Q(PHI_sqrt_SPD');
    matrix[m,m] QQT = tcrossprod(Q);
    vector[n_pred] theta_pred_mu = PHI_pred_sqrt_SPD * QQT * z_obs;
    matrix[n_pred,m] theta_pred_cov_H = PHI_pred_sqrt_SPD * add_diag(-QQT,1);
    vector[m] b = to_vector(normal_rng(rep_vector(0, m), rep_vector(1, m)));
		vector[n_pred] theta_pred = theta_pred_cov_H * b + theta_pred_mu;
    return theta_pred;
  }
}
data {
	int<lower=1> n; // number of observed data points
	int<lower=1> n_pred; // total number of points for prediction
	int<lower=1> p; // number of fixed covariates
  int<lower=1> k; // number of bins for piecewise linear baseline hazard rates
  vector[n] y; // time to event data
  vector[k+1] s; // knots for the piecewise linear baseline hazard rates
  matrix[n,2] delta; // indicator of risk type
	matrix[n,p] X; // design matrix - fixed effects
  vector[n] W; // design vector - for spatial coefficients
	matrix[n,2] d; // locations for observed data
	matrix[n_pred,2] d_pred; // locations for prediction
  real<lower=0> a0; // hyperparameter for psi0
  real<lower=0> b0; // hyperparameter for psi0
  real<lower=0> a1; // hyperparameter for kappa
  real<lower=0> b1; // hyperparameter for kappa
  real<lower=0> al; // hyperparameter for GP length scale parameter(s)
  real<lower=0> bl; // hyperparameter for GP length scale parameter(s)
  real<lower=0> tau; // hyperparameter for GP scale parameter(s)
	matrix[2,2] L0; // boundary condition
	array[2,2] int<lower=1> M0; // number of basis function of each dim
	matrix[2,2] L1; // boundary condition
	array[2,2] int<lower=1> M1; // number of basis function of each dim
}
transformed data {
  int m = 2;
  real epsilon = 0.000001;
	vector[m] zeros_m = zeros_vector(m);
  array[n] int<lower=1,upper=k> ki; // bin idx for each observation
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
	array[m] int<lower=1> M_nD0; // total number of basis functions for spatial intercept
	array[m] int<lower=1> M_nD1; // total number of basis functions for spatial slope
	for (j in 1:m){
	  M_nD0[j] = M0[1,j] * M0[2,j];
	  M_nD1[j] = M1[1,j] * M1[2,j];
	}
  matrix[M_nD0[1],2] idx_basis01; // grid for the basis functions for risk type 1
  matrix[M_nD0[2],2] idx_basis02; // grid for the basis functions for risk type 2
  matrix[M_nD1[1],2] idx_basis11; // grid for the basis functions for risk type 1
  matrix[M_nD1[2],2] idx_basis12; // grid for the basis functions for risk type 2
	int<lower=0> mm = 0;
	for (i in 1:M0[1,1]){
	  for (r in 1:M0[2,1]){
		  mm = mm + 1;
		  idx_basis01[mm,1] = i;
		  idx_basis01[mm,2] = r; 
		}
	}
	mm = 0;
	for (i in 1:M0[1,2]){
	  for (r in 1:M0[2,2]){
		  mm = mm + 1;
			idx_basis02[mm,1] = i;
			idx_basis02[mm,2] = r; 
		}
	}
	mm = 0;
	for (i in 1:M1[1,1]){
	  for (r in 1:M1[2,1]){
		  mm = mm + 1;
		  idx_basis11[mm,1] = i;
		  idx_basis11[mm,2] = r; 
		}
	}
	mm = 0;
	for (i in 1:M1[1,2]){
	  for (r in 1:M1[2,2]){
		  mm = mm + 1;
			idx_basis12[mm,1] = i;
			idx_basis12[mm,2] = r; 
		}
	}
	matrix[n,M_nD0[1]] PHI01 = PHI_2D(n, M0[1,1], M0[2,1], L0[1,1], L0[2,1], d[,1], d[,2]);
	matrix[n_pred,M_nD0[1]] PHI01_pred = PHI_2D(n_pred, M0[1,1], M0[2,1], L0[1,1], L0[2,1], d_pred[,1], d_pred[,2]);
	matrix[n,M_nD0[2]] PHI02 = PHI_2D(n, M0[1,2], M0[2,2], L0[1,2], L0[2,2], d[,1], d[,2]);
	matrix[n_pred,M_nD0[2]] PHI02_pred = PHI_2D(n_pred, M0[1,2], M0[2,2], L0[1,2], L0[2,2], d_pred[,1], d_pred[,2]);
	matrix[n,M_nD1[1]] PHI11 = PHI_2D(n, M1[1,1], M1[2,1], L1[1,1], L1[2,1], d[,1], d[,2]);
	matrix[n_pred,M_nD1[1]] PHI11_pred = PHI_2D(n_pred, M1[1,1], M1[2,1], L1[1,1], L1[2,1], d_pred[,1], d_pred[,2]);
	matrix[n,M_nD1[2]] PHI12 = PHI_2D(n, M1[1,2], M1[2,2], L1[1,2], L1[2,2], d[,1], d[,2]);
	matrix[n_pred,M_nD1[2]] PHI12_pred = PHI_2D(n_pred, M1[1,2], M1[2,2], L1[1,2], L1[2,2], d_pred[,1], d_pred[,2]);
}
parameters {
  matrix[p,m] beta;
  vector[m] beta_w;
  vector<lower=0>[m] psi0;
  matrix<lower=0>[m,k] psi;
 	vector<lower=0>[m] alpha0; // GP scale for intercept
	vector<lower=0,upper=10>[m] l0; // GP length scale for intercept
 	vector<lower=0>[m] alpha1; // GP scale for slope
	vector<lower=0,upper=10>[m] l1; // GP length scale for slope
  vector[M_nD0[1]] z01;
  vector[M_nD0[2]] z02;
  vector[M_nD1[1]] z11;
  vector[M_nD1[2]] z12;
  real<lower=0> s2;
  real<lower=0> kappa; // for multiplicative gamma process
}
transformed parameters {
  matrix[m,k] lambda; // baseline hazard rates
  lambda[,1] = psi0 .* psi[,1]; // psi0 corresponds to (imaginary) baseline hazard rate at time 0
  for (r in 2:k){
    lambda[,r] = psi[,r] .* lambda[,r-1]; // setup in this way because matrix is column major ordered in stan
  }
  vector[M_nD0[1]] sqrt_SPD01 = sqrt_diagSPD_2D(alpha0[1], l0[1], L0[,1], idx_basis01);
  vector[M_nD0[2]] sqrt_SPD02 = sqrt_diagSPD_2D(alpha0[2], l0[2], L0[,2], idx_basis02);
  vector[M_nD1[1]] sqrt_SPD11 = sqrt_diagSPD_2D(alpha1[1], l1[1], L1[,1], idx_basis11);
  vector[M_nD1[2]] sqrt_SPD12 = sqrt_diagSPD_2D(alpha1[2], l1[2], L1[,2], idx_basis12);
  vector[n] loglik;
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
    vector[n] theta01 = PHI01 * (sqrt_SPD01 .* z01);
    vector[n] theta02 = PHI02 * (sqrt_SPD02 .* z02);
    vector[n] theta11 = PHI11 * (sqrt_SPD11 .* z11);
    vector[n] theta12 = PHI12 * (sqrt_SPD12 .* z12);
    Xbeta[,1] = Xbeta[,1] + theta01 + (theta11 + beta_w[1]).*W;
    Xbeta[,2] = Xbeta[,2] + theta02 + (theta12 + beta_w[2]).*W;
    loglik = crcox_loglik(n,m,Xbeta', delta', lambda_cal, H);
  }
}
model {
  vector[k*m] psi_1d = to_vector(psi);
  vector[p*m] beta_1d = to_vector(beta);
  target += normal_lupdf(beta_1d | 0, sqrt(s2));
  target += normal_lupdf(beta_w | 0, sqrt(s2));
  target += inv_gamma_lupdf(s2 |1,1);
  target += gamma_lupdf(psi0 | a0, b0);
  target += gamma_lupdf(psi_1d | kappa*s_diff_inv_1d, kappa*s_diff_inv_1d);
  target += inv_gamma_lupdf(kappa | a1, b1);
  target += inv_gamma_lupdf(l0 |al, bl); // same prior as in demo
  target += inv_gamma_lupdf(l1 |al, bl); // same prior as in demo
  target += normal_lupdf(alpha0 |0, tau);
	target += normal_lupdf(alpha1 |0, tau);
  target += normal_lupdf(z01 |0,1);
  target += normal_lupdf(z02 |0,1);
  target += normal_lupdf(z11 |0,1);
  target += normal_lupdf(z12 |0,1);
  target += sum(loglik);
}
generated quantities {
  array[m] matrix[n_pred,2] theta_pred;
  vector<lower=0>[m] alpha0_adj;
  vector<lower=0>[m] alpha1_adj;
  {
    matrix[n,M_nD0[1]] PHI01_sqrt_SPD01 = diag_post_multiply(PHI01, sqrt_SPD01);
    matrix[n,M_nD0[2]] PHI02_sqrt_SPD02 = diag_post_multiply(PHI02, sqrt_SPD02);
    matrix[n,M_nD1[1]] PHI11_sqrt_SPD11 = diag_post_multiply(PHI11, sqrt_SPD11);
    matrix[n,M_nD1[2]] PHI12_sqrt_SPD12 = diag_post_multiply(PHI12, sqrt_SPD12);
    alpha0_adj[1] = sqrt(sum(square(PHI01_sqrt_SPD01))/n);
    alpha0_adj[2] = sqrt(sum(square(PHI02_sqrt_SPD02))/n);
    alpha1_adj[1] = sqrt(sum(square(PHI11_sqrt_SPD11))/n);
    alpha1_adj[2] = sqrt(sum(square(PHI12_sqrt_SPD12))/n);
    if (M_nD0[1]>n){
      matrix[n_pred,M_nD0[1]] PHI01_pred_sqrt_SPD01 = diag_post_multiply(PHI01_pred, sqrt_SPD01);
      theta_pred[1][,1] = theta_pred_rng(PHI01_sqrt_SPD01, PHI01_pred_sqrt_SPD01, z01);
    } else {
      theta_pred[1][,1] = PHI01_pred * (sqrt_SPD01 .* z01);
    }
    if (M_nD0[2]>n){
      matrix[n_pred,M_nD0[2]] PHI02_pred_sqrt_SPD02 = diag_post_multiply(PHI02_pred, sqrt_SPD02);
      theta_pred[2][,1] = theta_pred_rng(PHI02_sqrt_SPD02, PHI02_pred_sqrt_SPD02, z02); 
    } else {
      theta_pred[2][,1] = PHI02_pred * (sqrt_SPD02 .* z02);
    }
    if (M_nD1[1]>n){
      matrix[n_pred,M_nD1[1]] PHI11_pred_sqrt_SPD11 = diag_post_multiply(PHI11_pred, sqrt_SPD11);
      theta_pred[1][,2] = theta_pred_rng(PHI11_sqrt_SPD11, PHI11_pred_sqrt_SPD11, z11);
    } else {
      theta_pred[1][,2] = PHI11_pred * (sqrt_SPD11 .* z11);
    }
    if (M_nD1[2]>n){
      matrix[n_pred,M_nD1[2]] PHI12_pred_sqrt_SPD12 = diag_post_multiply(PHI12_pred, sqrt_SPD12);
      theta_pred[2][,2] = theta_pred_rng(PHI12_sqrt_SPD12, PHI12_pred_sqrt_SPD12, z12);      
    } else {
      theta_pred[2][,2] = PHI12_pred * (sqrt_SPD12 .* z12);
    }
  }
}
