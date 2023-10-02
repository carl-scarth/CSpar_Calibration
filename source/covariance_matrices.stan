// Contains expressions for the ARD covariance function which are used by 
// multiple other Stan scripts.

functions {
  
  // Function for defining a symmetric covariance matrix (covariance of set of
  // inputs, x, with itself)
  matrix ARD_SE_cov(matrix x, real lambda, row_vector beta, real delta) {
    // lambda = precision, beta = correlation length, delta = nugget
    int N = rows(x);    // number of training data points
    int d = cols(x);    // number of inputs
    matrix[N, N] K;
    {
    row_vector[d] dx;   // variable used to determine difference between inputs
    K = diag_matrix(rep_vector((1/lambda), N));
    for (i in 1:N) {
      for (j in (i+1):N) {
        dx = x[i] - x[j];
        K[i, j] = (beta .* dx) * dx';
        K[i, j] = exp(-K[i, j])/lambda;
        K[j, i] = K[i, j];
      }
    }
    K = K + diag_matrix(rep_vector(delta,N));
    }
    return K;
  }
  
  // Non-symmetric covariance matrix between two different sets of input values,
  // x1 and x2.
  matrix ARD_SE_cov_non_sym(matrix x1, matrix x2, real lambda, row_vector beta) {
    // lambda = precision, beta = correlation length, delta = nugget
    int N1 = rows(x1);  // number of training data points
    int N2 = rows(x2);  // number of points at which predictions are required
    int d = cols(x1);   // number of inputs
    matrix[N1, N2] K;
    {
    row_vector[d] dx;   // variable used to determine difference between inputs
    K = rep_matrix(0,N1,N2);
    for (i in 1:N1) {
      for (j in 1:N2) {
        dx = x1[i] - x2[j];
        K[i, j] = (beta .* dx) * dx';
        K[i, j] = exp(-K[i, j])/lambda;
      }
    }
    }
    return K;
  }
  
  // Covariance matrix between two sets of input values, x1 and x2, where x1 is 
  // a row vector with a single sample, and x2 is matrix of multiple samples
  matrix ARD_SE_cov_vect_mat(row_vector x1, matrix x2, real lambda, row_vector beta) {
    // lambda = precision, beta = correlation length, delta = nugget
    int N1 = rows(x1);  // number of training data points
    int N2 = rows(x2);  // number of points at which predictions are required
    int d = cols(x1);   // number of inputs
    matrix[N1, N2] K;
    {
    row_vector[d] dx;   // variable used to determine difference between inputs
    K = rep_matrix(0,N1,N2);
    for (i in 1:N1) {
      for (j in 1:N2) {
        dx = x1 - x2[j];
        K[i, j] = (beta .* dx) * dx';
        K[i, j] = exp(-K[i, j])/lambda;
      }
    }
    }
    return K;
  }
 
}
