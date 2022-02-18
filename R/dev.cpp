// Negative binomial lobster transect model:

#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() (){
   
   // Data variable declarations:
   DATA_VECTOR(y);                              // Count observations.
   DATA_IVECTOR(year);                          // Year indicators.
   
   // Fixed and random effect parameters:
   PARAMETER(alpha);                            // Global intercept parameter.
   PARAMETER_VECTOR(year_effect);               // Year random effect.
   
   // Log-scale error parameters:
   PARAMETER(log_sigma_year);                   // Log-scale error for year effect.
   
   // Log-scale error parameters:
   PARAMETER(log_r);                            // Precision parameter for negative binomial (log-scale).
   
   // Transformed parameters:
   int n_obs = y.size();                        // Number of observations.
   int n_year = year_effect.size();             // Number of year factor levels.
   Type r = exp(log_r);                         // Negative binomial precision parameter.
   Type sigma_year = exp(log_sigma_year);       // Standard error for year effect.
   Type res = 0;                                // Log-likelihood accumulator.
   
   // Random effects contributions:
   res -= sum(dnorm(year_effect, Type(0), sigma_year, true));
   
   // Observation log-likelihood:
   for(int i = 0; i < n_obs; i++){
      // Define log-linear mean:
      Type mu = exp(alpha + year_effect[year[i]]);
      
      // Negative binomial log-density:
      res -= lgamma(y[i]+r) - lgamma(r) - lgamma(y[i]+1) + r*log(r) + y[i]*log(mu) - (r+y[i])*log(r+mu);
   }
   
   return res;
}


