#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() (){
   // Data variable declarations:
   DATA_VECTOR(z);                              // Count observations.
   DATA_IVECTOR(len);                           // Length indicator.
   DATA_IVECTOR(year);                          // Year indicator.
   DATA_IVECTOR(region);                        // Region indicator.
   DATA_IVECTOR(transect);                      // Transect indicator.
   DATA_IVECTOR(diver);                         // Diver indicator
   DATA_VECTOR(p_sampled);                      // Proportion sampled.
   DATA_VECTOR(swept_area);                     // Swept area.
   
   // Parameter and random effects:
   PARAMETER(alpha);                            // Intercept parameter.
   PARAMETER(beta_sampled);                     // Sampling proportion coefficient.
   PARAMETER_VECTOR(length_effect);             // Length  effect.
   PARAMETER_VECTOR(year_effect);               // Year effect.
   PARAMETER_VECTOR(region_effect);             // Region effect.
   PARAMETER_VECTOR(transect_effect);           // Transect effect.
   PARAMETER_VECTOR(diver_effect);              // Diver effect.
   
   // Error parameters:
   PARAMETER(log_sigma_length);                 // Error for length effect.
   PARAMETER(log_sigma_year);                   // Error for year effect.
   PARAMETER(log_sigma_region);                 // Error for region effect.
   PARAMETER(log_sigma_transect);               // Error for transect effect.
   PARAMETER(log_sigma_diver);                  // Error for diver effect.
   PARAMETER(log_r);                            // Negative binomial dispersion parameter.
   
   // Auto-correlation parameter for random effects:
   PARAMETER(logit_phi_year);                                            // Correlation parameter for year effect.
   
   // Data and parameter dimensions:
   int n = z.size();                            // Number of observations.
  // int n_length = length_effect.size();         // Number of length bins.
//   int n_year = year_effect.size();             // Number of years.
 //  int n_region = region_effect.size();         // Number of regions.
//   int n_transect = transect_effect.size();     // Number of transects.
 //  int n_diver = diver_effect.size();           // Number of divers.
   
   // Transformations:
   Type res = 0;                                       // Log-likelihood accumulator.
   Type r = exp(log_r);                                // Negative binomial precision parameter.
   Type phi_year = 1.0 / (1.0 + exp(-logit_phi_year)); // Year effect correlation parameter.
   
   // Primary random effects:
   using namespace density;
   res -= sum(dnorm(length_effect, 0.0, exp(log_sigma_length), true));
   res += SCALE(AR1(phi_year), exp(log_sigma_year))(year_effect);        // Year effect.
   // res -= sum(dnorm(year_effect, 0.0, exp(log_sigma_year), true));
   res -= sum(dnorm(region_effect, 0.0, exp(log_sigma_region), true));
   res -= sum(dnorm(transect_effect, 0.0, exp(log_sigma_transect), true));
   res -= sum(dnorm(diver_effect, 0.0, exp(log_sigma_diver), true));
   
   // Observation model:
   for (int i = 0; i < n; i++){
      // Define mean:
      Type mu = exp(alpha + 
                    length_effect[len[i]] + 
                    year_effect[year[i]] + 
                    region_effect[region[i]] +
                    transect_effect[transect[i]] +
                    diver_effect[diver[i]] + 
                    beta_sampled * (1-p_sampled[i]) + 
                    log(swept_area[i]));
      
      // Negative binomial distribution:
      res -= lgamma(z[i]+r) - lgamma(r) - lgamma(z[i]+1) + r*log(r) + z[i]*log(mu) - (r+z[i])*log(r+mu);;
   }
   
   return res;
}



