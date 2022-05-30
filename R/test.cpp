// Negative binomial lobster transect model:

#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() (){
   // Data variable declarations:
   DATA_VECTOR(y);          // Count observations.

   // Parameter and random effect declarations:
   PARAMETER(alpha);                            // Global intercept parameter.
   PARAMETER(log_r);                            // Precision parameter for negative binomial (log-scale).

   // Transformed parameters:
   int n_obs = y.size();                        // Number of observations.

   Type r = exp(log_r);                         // Negative binomial precision parameter.
   Type mu;                                     // Regression mean variable.
   Type res = 0;                                // Log-likelihood accumulator.


   for(int i = 0; i < n_obs; i++){
      // Define log-linear mean:
      mu = exp(alpha);

      // Negative binomial distribution:
      res -= lgamma(y[i]+r) - lgamma(r) - lgamma(y[i]+1) + r*log(r) + y[i]*log(mu) - (r+y[i])*log(r+mu);
   }

   return res;
}
