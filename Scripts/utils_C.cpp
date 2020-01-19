#include <RcppArmadillo.h>
#include <armadillo>
#include <assert.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]



// [[Rcpp::export]]
bool any_naC(NumericVector x) {
  return is_true(any(is_na(x)));
}


// Model_Dist Function - returns probability of data 'y' under the marginal t-distribution from a Bayesian DLM
// [[Rcpp::export]]
double model_distC(const arma::vec& X, const arma::vec& m, const arma::mat& C, double n, double s, const NumericVector& y, bool mean) {
  //ft = t(X) %*% dlm$m
  //qt = as.numeric(t(X) %*% dlm$C %*% as.matrix(X)) + dlm$s
  double answer = 0.0;
  double ft = arma::dot(X, m);
  arma::vec qt = X.t() * C * X + s;
  
  if(mean){
    return ft;
  }
  
  if(any_naC(y)){
    arma::vec tmp = rt(1, n)*sqrt(qt(0)) + ft;
    answer = (double) tmp(0);
    return answer;
  }else{
    arma::vec tmp = dt((y - ft)/sqrt(qt(0)), n)/sqrt(qt(0));
    answer = (double) tmp(0);
    return answer;
  }
}

// Helper function to extract elements from a matrix
// [[Rcpp::export]]
arma::vec matrix_locs(arma::mat M, arma::umat locs) {
  
  arma::uvec eids = sub2ind( size(M), locs ); // Obtain Element IDs
  arma::vec v  = M.elem( eids );              // Values of the Elements
  
  return v;
}

// Alternative helper function to extract elements of a matrix
// [[Rcpp::export]]
arma::mat submat(arma::mat M, arma::uvec rows, arma::uvec cols){
  arma::mat subm = arma::mat(rows.size(), cols.size());
  subm = M.submat(rows, cols);
  return subm;
  
}

// Model_Dist wrapper function that evaluates probability of data 'y'
// [[Rcpp::export]]
double model_dist_dlmC(int time, arma::mat F, arma::vec model, arma::vec m, arma::mat C, double n, double s, NumericVector y, bool mean){
  arma::umat locs = arma::umat(2, model.size());
  for(int i = 0; i < model.size(); ++i){
    locs(0, i) = time - 1; // Subtracting 1 to convert to C++ indices
    locs(1, i) = model(i) - 1; // Subtracting 1 to convert to C++ indices
  }

  arma::vec X = matrix_locs(F, locs);
  return model_distC(X, m, C, n, s, y, mean);
}


// Alternative wrapper function that evaluates probability of data 'y'
// [[Rcpp::export]]
double model_dist_dlmC2(arma::uvec time, arma::mat F, arma::uvec model, arma::vec m, arma::mat C, double s, double n, NumericVector y, bool mean){
  arma::vec X = vectorise(submat(F, time - 1, model - 1));
  return model_distC(X, m, C, n, s, y, mean);
}

// [[Rcpp::export]]
LogicalVector logical_index(arma::uvec idx, unsigned int n) {
  LogicalVector result(n, false);

  for (int i = 0; i < idx.size(); i++) {
    result[idx[i]] = true;
  }

  return result;
}


// Initialize a DLM
// [[Rcpp::export]]
List initialize_dynamic_linear_modelC(const arma::mat& F_prior, const arma::vec& y_prior, arma::uvec model, double delta, double beta){
  model = model - 1; // Subtracting 1 to convert to C++ indices
  int p = model.size();
  int initial_points = y_prior.size();
  
  // Initialize the dlm parameters
  double n = 1;
  double s = .1;
  vec m(p); m.zeros();
  mat C(p, p); C.eye();
  C = C*10;
  
  vec X;
  umat locs(2, model.size());
  double ft, qt;
  NumericVector obs(1), error(1);
  
  for(int t = 0; t < initial_points; ++t){
    n = delta * n;
    C = C / beta;
    
    
    // Selecting out my X data
    for(int i = 0; i < model.size(); ++i){
      locs(0, i) = t;
      locs(1, i) = model(i);
    }
    
    X = matrix_locs(F_prior, locs);
    
    ft = arma::dot(X, m);
    arma::vec tmp1 = X.t() * C * X + s;
    qt = tmp1(0);
    
    // Evaluating observation likelihood
    obs = y_prior(t);
    error = obs - ft;

    // Kalman filter updating
    arma::mat At = C * X / qt;
    double rt = (n + pow(error(0), 2)/qt)/(n+1);
    ++n;
    s = s * rt;
    m = m + At * as<arma::vec>(error);
    C = rt * (C - qt * At * At.t());
    
    
  }
  

  
  
  
  return List::create(_["model"] = model, _["m"] = m, _["C"] = C, _["n"] = n, _["s"] = s, _["time"] = 0);
}

// Dynamic_Linear_Model Function - Kalman filtering, stores model likelihood aka 1-step ahead forecast densities
// [[Rcpp::export]]
List dynamic_linear_model_updateC(double time, const arma::mat& F, const arma::vec& y, arma::uvec model, arma::vec m, arma::mat C, double n, double s, NumericVector likelihood, int start_time, double delta, double beta, bool intervention, vec interventions){
  
  model = model - 1; // Subtracting 1 to convert to C++ indices

  arma::vec X;
  arma::umat locs(2, model.size());
  double ft, qt;
  NumericVector obs(1), error(1);
  //NumericVector obs;
  if(time > 0){
    for(int t = start_time; t < time; ++t){
 
      if(intervention){
        for(int i = 0; i < interventions.size(); i++){
          if(t == (int) (interventions(i) - 1)){ // Subtracting 1 to convert to C++ indices
            n = 0.5 * n / delta;
            C = C / 0.5 * beta;
          }
        }
      }
      
      n = delta * n;
      C = C / beta;
      
      // Selecting out my X data
      for(int i = 0; i < model.size(); ++i){
        locs(0, i) = t;
        locs(1, i) = model(i);
      }

      X = matrix_locs(F, locs);

      // Rcout << X << std::endl;
      ft = arma::dot(X, m);
      arma::vec tmp1 = X.t() * C * X + s;
      qt = tmp1(0);
      
      // Evaluating observation likelihood
      obs = y(t);
      error = obs - ft;
      NumericVector tmp = dt(error/sqrt(qt), n)/sqrt(qt);
      likelihood(t) = (double) tmp(0);
      
      // Kalman filter updating
      arma::mat At = C * X / qt;
      double rt = (n + pow(error(0), 2)/qt)/(n+1);
      ++n;
      s = s * rt;
      m = m + At * as<arma::vec>(error);
      C = rt * (C - qt * At * At.t());
    }
  }
  
  return List::create(_["model"] = model, _["m"] = m, _["C"] = C, _["n"] = n, _["s"] = s, _["time"] = time);
}

// Input DLM components and evaluate a forecast path of length k
// Evaluate the log-probability of of sequence of 'true' values of length k
// This is a univariate setting, can be used in a DDNM model to evaluate the LPFDS score for a single series' model
// Use true values for the other series, rather than simulating from a full DDNM (combination of many univariate DLMs, plus contemporaneous predictors)
// [[Rcpp::export]]
double LPFDSC(vec m, mat C, double n, double s, uvec model, int k, int series_to_eval, int time_to_eval, mat data, vec observations, double delta, double beta){

  model = model - 1;
  vec true_log_likelihoods(k);
  NumericVector likelihood(time_to_eval);
  arma::umat locs(2, model.size());
  arma::vec X;
  double ft;
  double g1; double g2;
  
  // Initialize samples for the precision and the state evolution variance
  double v = 1/R::rgamma(n/2, 2/(n*s));
  mat W = ((1 - delta)/delta) * C;
  mat L = chol(W).t();
  int p = m.size();
  vec theta = m + sqrt(v / s) * chol(C).t() * randn(p);

  for(int t = 0; t < k; ++t){
    
    // Simulate the precisions and the states
    
    // Simulating a beta RV as fxn of 2 gammas:
    // https://en.wikipedia.org/wiki/Beta_distribution#Generating_beta-distributed_random_variates
    g1 = R::rgamma(beta*n/2, 1);
    g2 = R::rgamma((1-beta)*n/2, 1);
    
    v = v * (beta / (g1 / (g1 + g2)));

    theta = theta + sqrt(v / s) * L * randn(p);

    // # Log 1-step ahead predictive distribution of observed values
    // Selecting out X data
    for(int i = 0; i < model.size(); ++i){
      locs(0, i) = time_to_eval - k + t;
      locs(1, i) = model(i);
    }
    X = matrix_locs(data, locs);
    ft = dot(X, theta);
    true_log_likelihoods(t) = log(R::dnorm(observations(time_to_eval - k + t), ft, sqrt(v), false));  
  }

  return(sum(true_log_likelihoods));
}

// Wrapper function for simulate_path_LPFDSC_univar that calculates a model score through all times
// Up to time t
// Being efficient by only initializing the dlm once, and updating only as needed
// [[Rcpp::export]]
double LPFDSC_MC(int nsamps, uvec model, int k, int series_to_eval, int start_time, int end_time, mat data, vec observations, double delta, double beta, double alpha, bool intervention, vec interventions, const arma::mat& F_prior, const arma::vec& y_prior){

  // Initialize the DLM at time 0
  List dlm_init = initialize_dynamic_linear_modelC(F_prior, y_prior, model, delta, beta);
  
  vec m = as<arma::vec>(dlm_init["m"]);
  mat C = as<mat>(dlm_init["C"]);
  double n = as<double>(dlm_init["n"]);
  double s = as<double>(dlm_init["s"]);
  List dlm = dlm_init;
  int dlm_time = as<int>(dlm_init["time"]);

  NumericVector likelihood(end_time);
  vec scores(end_time - start_time + 1); scores.zeros();

  // Setting the weights - using exponential discounting
  vec weights(end_time - start_time + 1);
  weights.fill(1);
  for(int i = 0; i < (end_time - start_time); i++){
    weights.head(end_time - start_time - i) = weights.head(end_time - start_time - i) * alpha;
  }

  // Scale the number of samples used by the weights - don't need many samples for old, discounted times
  vec samps = floor(weights * nsamps);
  for(int time_to_eval = start_time; time_to_eval < (end_time + 1); time_to_eval++){
    // Update the dynamic linear model to the evaluation point
    dlm_time = as<int>(dlm["time"]);
    // Rcout << dlm_time << std::endl;
    
    dlm = dynamic_linear_model_updateC((double) time_to_eval - k,
                                            data,
                                            observations,
                                            model,
                                            m,
                                            C,
                                            n,
                                            s,
                                            likelihood,
                                            dlm_time,
                                            delta,
                                            beta,
                                            intervention,
                                            interventions);

    m = as<arma::vec>(dlm["m"]);
    C = as<mat>(dlm["C"]);
    n = as<double>(dlm["n"]);
    s = as<double>(dlm["s"]);

    if(samps(time_to_eval - start_time) > 0){ // Don't evaluate if we take 0 samples
      //  Not being parallel in here, because I'm parallel at the MODEL level, not at the nsamps level or the time level
      vec scores_MC(samps(time_to_eval - start_time)); scores_MC.zeros();
      for(int i = 0; i < samps(time_to_eval - start_time); i++){
        scores_MC(i) = LPFDSC(m, C, n, s, model, k, series_to_eval, time_to_eval, data, observations, delta, beta);
      }
      scores(time_to_eval - start_time) = log(mean(exp(scores_MC)));
    }

  }

  // Adding in the intervention - within this function, these are correct, do not need to subtract 1 compared to the interventions in R
  if(intervention){
    for(int i = 0; i < interventions.size(); i++){
      if(start_time < interventions(i) & end_time > interventions(i)){
        scores.head(interventions(i)-start_time) = 0.5 * scores.head(interventions(i)-start_time);
      }
    }
  }

  // Rcout << scores << std::endl;

  return(dot(scores, weights));
}


// # Function to construct the VAR representation of a DDNM
// # Decompose that representation to test for stationarity
// Need to pass in a theta_list, model_list, and p_list with the intercept removed
// [[Rcpp::export]]
bool is_stationaryC(std::vector<vec> theta_list, std::vector<uvec> model_list, std::vector<int> p_list, int max_lag, int num_series){
// # Initialize the dense G matrix
  mat G = mat(max_lag*num_series, max_lag*num_series, fill::zeros);
  mat coefs = mat(num_series, (1+max_lag)*num_series, fill::zeros);
// # Initialize the coefficients & contemporaneous predictors from each model filled in with zeroes
// # Ignore the intercept
  for(int i = 0; i < num_series; i++){
    for(int j = 0; j < p_list[i]; j++){
      coefs(i, model_list[i](j) + i - 1) = theta_list[i](j);
    }
  }
  

// # Set up the contemporaneous predictor matrix and get I-Gam
    mat IGam(num_series, num_series); IGam.eye();

    IGam -= coefs.head_cols(num_series);
    
// # Go from DDNM to VAR representation for coefficient matrices by multiplying coef =  (I-Gam)^-1 * coef for each coef matrix
// # Record the locations & values for a sparse matrix representation as well
    G.head_rows(num_series) = solve(IGam, coefs.tail_cols(max_lag*num_series));
    
// # Add in the identity matrices to G
    G(span(num_series, num_series*max_lag-1), span(0, num_series*max_lag - num_series - 1)).diag().ones();

    
// # I only need the LARGEST eigenvalue - which luckily I can extract just that
    cx_vec eigval = eig_gen(G);

// # Test for stationarity
  if(max(abs(eigval)) < 1){
    return(true);
  } else{
    return(false);
  }
}

// Wrapper function that calculates the KL between a model's forecast dist. and sampled values
// This uses the LPFDS function, but instead of using just a vector of true, future values
// It accepts a data cube coming from some other forecast dist (in this case, the model averaged forecast)
// [[Rcpp::export]]
double KL_indep(int nsamps, const arma::vec& m, const arma::mat& C, double n, double s, uvec model, int k, int series_to_eval, int time, cube data, vec observations, mat future_values, int num_series, int ncol, double delta, double beta){
  
  mat dat(data.n_rows, data.n_cols, fill::zeros);
  vec full_obs(time, fill::zeros);
  full_obs.head(time - k) = observations;
  
  NumericVector likelihood(time);
  
  //  Not being parallel in here, because I'm parallel at the MODEL level, not at the nsamps level or the time level
  vec scores_MC(nsamps); scores_MC.zeros();
  int replicates = 10;
  vec scores(replicates);
  for(int i = 0; i < nsamps; i++){
    dat = data.slice(i);
    full_obs.tail(k) = future_values.col(i);
    
    // For each synthetic future generated by the mixture of models
    // Doing 10 MC draws from each individual model to get a sense of how closely they match
    for(int j = 0; j < replicates; j++){
      scores(j) = LPFDSC(m, C, n, s, model, k, series_to_eval, time, dat, full_obs, delta, beta);
    }
    
    scores_MC(i) = log(mean(exp(scores)));
  }

  return(log(mean(exp(scores_MC))));
}
