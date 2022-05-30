// Negative binomial lobster transect model:

#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() (){
   
   // Data variable declarations:
   DATA_VECTOR(y);          // Count observations.
   
   // Fixed and random effect parameters:
   PARAMETER(alpha);                            // Global intercept parameter..
   
   // Log-scale error parameters:
   PARAMETER(log_r);                            // Precision parameter for negative binomial (log-scale).
   
   // Transformed parameters:
   int n_obs = y.size();                        // Number of observations.
   Type r = exp(log_r);                         // Negative binomial precision parameter.
   Type res = 0;                                // Log-likelihood accumulator.
   
   // Observation log-likelihood:
   for(int i = 0; i < n_obs; i++){
      // Define log-linear mean:
      Type mu = exp(alpha + log(area[i]) - log(100));
      
      // Negative binomial log-density:
      res -= lgamma(y[i]+r) - lgamma(r) - lgamma(y[i]+1) + r*log(r) + y[i]*log(mu) - (r+y[i])*log(r+mu);
   }
   
   return res;
}


