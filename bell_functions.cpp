#include <Rcpp.h>
#include <vector>
#include <cmath>

// Helper function to calculate Stirling numbers of the second kind
double Stirling2(int n, int k) {
  if (k == 0 && n == 0) return 1;
  if (k == 0 || n == 0) return 0;
  return k * Stirling2(n - 1, k) + Stirling2(n - 1, k - 1);
}

// C++ version of the Bell function
// [[Rcpp::export]]
Rcpp::NumericVector Bell(Rcpp::IntegerVector n) {
  int nN = n.size();
  Rcpp::NumericVector result(nN);
  
  for (int r = 0; r < nN; ++r) {
    if (n[r] == 0) {
      result[r] = 1;  // Define Bell(0) = 1
    } else if (n[r] < 0) {
      Rcpp::stop("n must be greater than or equal to 0");
    } else {
      double sum = 0;
      for (int k = 1; k <= n[r]; ++k) {
        sum += Stirling2(n[r], k);
      }
      result[r] = sum;
    }
  }
  
  return result;
}

// C++ version of the dbell function
// [[Rcpp::export]]
Rcpp::NumericVector dbell(Rcpp::NumericVector x, double shape = 1, bool log_prob = false) {
  int n = x.size();
  Rcpp::NumericVector logdensity(n, R_NegInf);  // Initialize with log(0)
  
  for (int i = 0; i < n; ++i) {
    if (x[i] >= 0 && std::isfinite(x[i]) && x[i] == std::round(x[i])) {
      if (shape > 0) {
        Rcpp::NumericVector bell_val = Bell(Rcpp::IntegerVector::create(x[i]));
        double ld = x[i] * std::log(shape) - std::expm1(shape) +
          std::log(bell_val[0]) - std::lgamma(x[i] + 1);
        logdensity[i] = ld;
      } else {
        logdensity[i] = NAN;  // Set to NaN if shape is not positive
      }
    }
  }
  
  // Return either the log-density or the density depending on log_prob
  if (log_prob) {
    return logdensity;
  } else {
    return Rcpp::exp(logdensity);
  }
}
// Helper function to generate random positive Poisson variates
int rpospois(double lambda) {
  if (lambda <= 0) return 1;
  int count = 0;
  while (true) {
    count = R::rpois(lambda);
    if (count > 0) return count;
  }
}

// C++ version of the rbell function
// [[Rcpp::export]]
Rcpp::IntegerVector rbell(int n, double shape = 1) {
  if (n <= 0) {
    Rcpp::stop("Argument 'n' must be a positive integer.");
  }
  
  Rcpp::IntegerVector ans(n, 0);
  Rcpp::IntegerVector N(n);
  Rcpp::LogicalVector Nok(n);
  
  for (int i = 0; i < n; ++i) {
    N[i] = R::rpois(std::expm1(shape));  // expm1(shape) returns exp(shape) - 1
    Nok[i] = !Rcpp::NumericVector::is_na(shape) && N[i] > 0;
  }
  
  int maxN = *std::max_element(N.begin(), N.end());
  
  for (int kk = 1; kk <= maxN; ++kk) {
    for (int i = 0; i < n; ++i) {
      if (Nok[i] && kk <= N[i]) {
        ans[i] += rpospois(shape);
      }
    }
  }
  
  for (int i = 0; i < n; ++i) {
    if (shape <= 0) {
      ans[i] = NAN;
    }
    if (Rcpp::NumericVector::is_na(shape)) {
      ans[i] = NA_INTEGER;
    }
  }
  
  return ans;
}
// Improved pbell function in C++
// [[Rcpp::export]]
Rcpp::NumericVector pbell(Rcpp::IntegerVector q, double theta, bool lower_tail = true, bool log_p = false) {
  int n = q.size();
  Rcpp::NumericVector result(n);
  
  for (int i = 0; i < n; ++i) {
    if (q[i] < 0) {
      result[i] = log_p ? R_NegInf : 0.0; // Log(0) or 0 probability for negative q
      continue;
    }
    
    double cumulative_prob = 0.0;
    for (int k = 0; k <= q[i]; ++k) {
      cumulative_prob += dbell(Rcpp::NumericVector::create(k), theta, false)[0];
    }
    
    if (!lower_tail) {
      cumulative_prob = 1.0 - cumulative_prob;
    }
    
    result[i] = log_p ? std::log(cumulative_prob) : cumulative_prob;
  }
  
  return result;
}
// Function to compute the quantile function of the Bell distribution
// [[Rcpp::export]]
Rcpp::IntegerVector qbell(Rcpp::NumericVector p, double theta, bool log_p = false) {
  int n = p.size();
  Rcpp::IntegerVector result(n);
  
  for (int i = 0; i < n; ++i) {
    double target_p = log_p ? std::exp(p[i]) : p[i];
    double cumulative_p = 0;
    int k = 0;
    
    // Ensure the target probability is within the valid range [0, 1]
    if (target_p < 0 || target_p > 1) {
      result[i] = -1; // Sentinel value for invalid probability
      continue; // Skip further computation for this element
    }
    
    // Increment k until the cumulative probability exceeds the target
    while (cumulative_p < target_p) {
      cumulative_p = pbell(Rcpp::IntegerVector::create(k), theta, true, false)[0];
      ++k;
      
      // Safety check: Avoid infinite loop in case of unexpected behavior
      if (k > 10000) { // Arbitrary large number to prevent infinite loop
        Rcpp::stop("Exceeded iteration limit in qbell.");
      }
    }
    
    result[i] = k - 1; // Adjust k to get the correct quantile
  }
  
  return result;
}
// Lambert W function for the principal branch (W_0)
// Using Halley's method for numerical approximation
double lambertW0(double z) {
  if (z == 0.0) {
    return 0.0;
  }
  if (z == -1.0 / std::exp(1.0)) {
    return -1.0;
  }
  if (z < -1.0 / std::exp(1.0)) {
    return NAN; // No real solution exists
  }
  
  double w = z;
  double tolerance = 1e-10;
  for (int i = 0; i < 100; ++i) {
    double ew = std::exp(w);
    double wew = w * ew;
    double wewz = wew - z;
    double w1 = w + 1.0;
    double delta_w = wewz / (ew * w1 - (w + 2.0) * wewz / (2.0 * w1));
    w -= delta_w;
    if (std::abs(delta_w) < tolerance) {
      return w;
    }
  }
  return w; // Return the last approximation if convergence not achieved
}

// Expose the function to R
// [[Rcpp::export]]
double lambertW(double z) {
  return lambertW0(z);
}
// Function to compute quantile residuals for a Bell regression model
// [[Rcpp::export]]
Rcpp::NumericVector qresiduals(Rcpp::List object) {
  Rcpp::NumericVector y = object["y"];
  Rcpp::NumericVector mu = object["mu"];
  int n = y.size();
  Rcpp::NumericVector residuals(n);
  
  for (int i = 0; i < n; ++i) {
    double theta = lambertW(mu[i]);  // Calculate theta using Lambert W function
    
    // Calculate the probabilities for the CDF at y[i] - 1 and y[i]
    double a = pbell(Rcpp::IntegerVector::create(y[i] - 1), theta, true, false)[0];
    double b = pbell(Rcpp::IntegerVector::create(y[i]), theta, true, false)[0];
    
    // Ensure that a is non-negative and less than b
    if (a < 0) a = 0;
    if (a > b) std::swap(a, b);
    
    // Generate a uniform random number between a and b
    double u = R::runif(a, b);
    
    // Calculate the quantile residuals using the normal quantile function
    residuals[i] = R::qnorm(u, 0.0, 1.0, true, false);
  }
  
  return residuals;
}

// Improved log-likelihood function for the Bell distribution
// [[Rcpp::export]]
double log_likelihood_bell(Rcpp::IntegerVector x, double theta) {
  int n = x.size();
  double log_likelihood = 0.0;
  
  Rcpp::NumericVector bell_numbers = Bell(x);  // Calculate Bell numbers for all x[i]
  
  for (int i = 0; i < n; ++i) {
    if (x[i] < 0) {
      Rcpp::stop("x must be non-negative integers.");
    }
    
    // Calculate the log-likelihood component
    double lf = x[i] * std::log(theta) - std::expm1(theta) + 
      std::log(bell_numbers[i]) - std::lgamma(x[i] + 1);
    
    // Ensure lf is finite before adding it
    if (std::isfinite(lf)) {
      log_likelihood += lf;
    } else {
      Rcpp::Rcout << "Warning: Non-finite log-likelihood component detected for x[" 
                  << i << "] = " << x[i] << std::endl;
    }
  }
  
  return log_likelihood;
}
// Function to compute the variance function of mu for the Bell distribution
Rcpp::NumericVector vmu(Rcpp::NumericVector mu) {
  int n = mu.size();
  Rcpp::NumericVector variance(n);
  
  for (int i = 0; i < n; ++i) {
    if (mu[i] < 0) {
      Rcpp::stop("mu must be non-negative.");
    }
    
    double lambert_w0 = lambertW0(mu[i]);
    variance[i] = mu[i] * (1.0 + lambert_w0);
  }
  
  return variance;
}

// Rcpp export for vectorized vmu function
// [[Rcpp::export]]
Rcpp::NumericVector cpp_vmu(Rcpp::NumericVector mu) {
  return vmu(mu);
}
