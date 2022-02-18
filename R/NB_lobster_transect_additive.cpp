// Negative binomial lobster transect model:

#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
   // Data variable declarations:
   DATA_VECTOR(y);          // Count observations.
   DATA_VECTOR(area);       // Area coverage.
   DATA_FACTOR(year);       // Year indicators.
   DATA_FACTOR(region);     // Region indicators.
   DATA_FACTOR(cohort);     // Cohort indicators.
   DATA_FACTOR(transect);   // Transect indicators.

   // Parameter and random effect declarations:
   PARAMETER(alpha);                            // Global intercept parameter.
   PARAMETER_VECTOR(year_effect);               // Year random effect.
   PARAMETER_VECTOR(region_effect);             // Region random effect.
   PARAMETER_VECTOR(cohort_effect);             // Cohort random effect.
   PARAMETER_VECTOR(transect_effect);           // Transect random effect.
   PARAMETER(log_sigma_year);                   // Log-scale error for year effect.
   PARAMETER(log_sigma_region);                 // Log-scale error for region effect.
   PARAMETER(log_sigma_cohort);                 // Log-scale error for cohort effect.
   PARAMETER(log_sigma_transect);               // Log-scale error for transect effect.
   PARAMETER(log_r);                            // Precision parameter for negative binomial (log-scale).

   // Transform random effect error parameters:
   Type sigma_year = exp(log_sigma_year);                             // Standard error for year effect.
   Type sigma_region = exp(log_sigma_region);                         // Standard error for region effect.
   Type sigma_cohort = exp(log_sigma_cohort);                         // Standard error for cohort effect.
   Type sigma_transect = exp(log_sigma_transect);                     // Standard error for transect effect.
   
   // Transformed parameters:
   int n_obs = y.size();                        // Number of observations.
   int n_year = year_effect.size();             // Number of year factor levels.
   int n_region = region_effect.size();         // Number of region factor levels.
   int n_cohort = cohort_effect.size();         // Number of cohort factor levels.  
   int n_transect = transect_effect.size();     // Number of transect factor levels.  
   Type r = exp(log_r);                         // Negative binomial precision parameter.
   Type muhat;                                  // Regression mean temporary variable.
   Type res = 0;                                // Log-likelihood accumulator.
   Type pi = 3.14159265359;                     // Define pi.
   
   // Priors over random effect standard errors (sigma ~ Half-Cauchy(0,5)):
   res -= log(2 * 5) - log(pi) - log(5*5 + (sigma_year * sigma_year));
   res -= log(2 * 5) - log(pi) - log(5*5 + (sigma_region * sigma_region));
   res -= log(2 * 5) - log(pi) - log(5*5 + (sigma_cohort * sigma_cohort));
   res -= log(2 * 5) - log(pi) - log(5*5 + (sigma_transect * sigma_transect));
   
   // Hierarchical prior contributions to log-likelihood:
   res -= sum(dnorm(year_effect, Type(0), sigma_year, true));
   res -= sum(dnorm(region_effect, Type(0), sigma_region, true));
   res -= sum(dnorm(cohort_effect, Type(0), sigma_cohort, true));
   res -= sum(dnorm(transect_effect, Type(0), sigma_transect, true));
               
   // Compute categorical mean values:
   array<Type> log_mu(n_year, n_region, n_cohort);
   for(int i = 0; i < n_year; i++){
      for(int j = 0; j < n_region; j++){
         for(int k = 0; k < n_cohort; k++){
            log_mu(i,j,k) = alpha +
                            year_effect[i] +
                            region_effect[j] +
                            cohort_effect[k];             
         }
      }
   }
   
   // Observation model:
   for(int i = 0; i < n_obs; i++){
      // Define mean:
      muhat = exp(log_mu(year[i]-1, region[i]-1, cohort[i]-1) + 
                  transect_effect[transect[i]-1] +
                  log(area[i]) - log(100));
      
      // Negative binomial distribution:
      res -= lgamma(y[i]+r) - lgamma(r) - lgamma(y[i]+1) + r*log(r) + y[i]*log(muhat) - (r+y[i])*log(r+muhat);
   }

   return res;
}
