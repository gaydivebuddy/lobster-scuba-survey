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
   DATA_IMATRIX(mls);                           // Minimum legal sizes (n_year x n_region).
   
   // Parameter and random effects:
   PARAMETER(alpha);                            // Intercept parameter.
   PARAMETER(beta_sampled);                     // Sampling proportion coefficient.
   PARAMETER_VECTOR(length_effect);             // Length  effect.
   PARAMETER_VECTOR(year_effect);               // Year effect.
   PARAMETER_VECTOR(region_effect);             // Region effect.
   PARAMETER_VECTOR(transect_effect);           // Transect effect.
   PARAMETER_VECTOR(diver_effect);              // Diver effect.
   
   // Interaction random effect terms:
   PARAMETER_VECTOR(length_year_effect);        // Length x year interaction effect.
   PARAMETER_VECTOR(region_year_effect);        // Region x year interaction effect.
   PARAMETER_VECTOR(transect_year_effect);      // Transect x year interaction effect.
   PARAMETER_VECTOR(diver_year_effect);         // Diver x year interaction effect.
   PARAMETER_VECTOR(length_region_effect);      // Length x region interaction effect.
   PARAMETER_VECTOR(length_diver_effect);       // Length x diver interaction effect.
   PARAMETER_VECTOR(length_year_region_effect); // Length x year x region interaction effect.
   
   // Error parameters for random effects:
   PARAMETER(log_sigma_length);                 // Error for length effect.
   PARAMETER(log_sigma_year);                   // Error for year effect.
   PARAMETER(log_sigma_region);                 // Error for region effect.
   PARAMETER(log_sigma_transect);               // Error for transect effect.
   PARAMETER(log_sigma_diver);                  // Error for diver effect.
   
   // Error parameters for interaction random effects:
   PARAMETER(log_sigma_length_year);            // Error for length x year interaction effect.
   PARAMETER(log_sigma_region_year);            // Error for region x year interaction effect.
   PARAMETER(log_sigma_length_region);          // Error for length x region interaction effect.
   PARAMETER(log_sigma_transect_year);          // Error for transect x year interaction effect.
   PARAMETER(log_sigma_diver_year);             // Error for diver x year interaction effect.
   PARAMETER(log_sigma_length_diver);           // Error for length x diver interaction effect.
   PARAMETER(log_sigma_length_year_region);     // Error for length x year x region interaction effect.  
      
   // Auto-correlation parameter for random effects:
   PARAMETER(logit_phi_year);                   // Correlation parameter for year effect.
   PARAMETER(logit_phi_length_year);            // Correlation parameter for length x year effect. 
   PARAMETER(logit_phi_region_year);            // Correlation parameter for region x year effect. 
   PARAMETER(logit_phi_transect_year);          // Correlation parameter for transect x year effect. 

   // Negative binomial dispersion parameter.
   PARAMETER(log_r);                            
   
   // Data and parameter dimensions:
   int n = z.size();                            // Number of observations.
   int n_length = length_effect.size();         // Number of length bins.
   int n_year = year_effect.size();             // Number of years.
   int n_region = region_effect.size();         // Number of regions.
   int n_transect = transect_effect.size();     // Number of transects.
   int n_diver = diver_effect.size();           // Number of divers.
   
   // Transformations:
   Type res = 0;                                       // Log-likelihood accumulator.
   Type r = exp(log_r);                                // Negative binomial precision parameter.
   Type phi_year = 1.0 / (1.0 + exp(-logit_phi_year));                   // Year effect correlation parameter.
   Type phi_length_year = 1.0 / (1.0 + exp(-logit_phi_length_year));     // Length x year effect correlation parameter. 
   Type phi_transect_year = 1.0 / (1.0 + exp(-logit_phi_transect_year)); // Transect x year effect correlation parameter.  
   Type phi_region_year = 1.0 / (1.0 + exp(-logit_phi_region_year));     // Region x year effect correlation parameter.  
   
   using namespace density;
   
   // Primary random effects:
   res -= sum(dnorm(length_effect, 0.0, exp(log_sigma_length), true));
   res += SCALE(AR1(phi_year), exp(log_sigma_year))(year_effect);
   res -= sum(dnorm(region_effect, 0.0, exp(log_sigma_region), true));
   res -= sum(dnorm(transect_effect, 0.0, exp(log_sigma_transect), true));
   res -= sum(dnorm(diver_effect, 0.0, exp(log_sigma_diver), true));
   
   // Interaction random effects:
   // res -= sum(dnorm(length_year_effect, 0.0, exp(log_sigma_length_year), true));
   for (int i = 0; i < n_length; i++){
      res += SCALE(AR1(phi_length_year), exp(log_sigma_length_year))(length_year_effect.segment(i * n_year, n_year));
   }                                                                                         
   res -= sum(dnorm(length_region_effect, 0.0, exp(log_sigma_length_region), true));
   // res -= sum(dnorm(year_region_effect, 0.0, exp(log_sigma_year_region), true));
   for (int i = 0; i < n_region; i++){
      res += SCALE(AR1(phi_region_year), exp(log_sigma_region_year))(region_year_effect.segment(i * n_year, n_year));
   }
   res -= sum(dnorm(length_diver_effect, 0.0, exp(log_sigma_length_diver), true));
   res -= sum(dnorm(diver_year_effect, 0.0, exp(log_sigma_diver_year), true));
   //res -= sum(dnorm(year_transect_effect, 0.0, exp(log_sigma_year_transect), true));
   for (int i = 0; i < n_transect; i++){
      res += SCALE(AR1(phi_transect_year), exp(log_sigma_transect_year))(transect_year_effect.segment(i * n_year, n_year));
   }
   res -= sum(dnorm(length_year_region_effect, 0.0, exp(log_sigma_length_year_region), true));
   
   // Observation model:
   for (int i = 0; i < n; i++){
      // Define mean:
      Type mu = exp(alpha + 
                    length_effect[len[i]] + 
                    year_effect[year[i]] + 
                    region_effect[region[i]] +
                    transect_effect[transect[i]] +
                    diver_effect[diver[i]] + 
                    length_year_effect[len[i] * n_year + year[i]] +
                    length_region_effect[len[i] * n_region + region[i]] +
                    region_year_effect[region[i] * n_year + year[i]] +
                    length_diver_effect[len[i] * n_diver + diver[i]] + 
                    diver_year_effect[diver[i] * n_year + year[i]] +
                    transect_year_effect[transect[i] * n_year + year[i]] +
                    length_year_region_effect[(len[i] * n_year + year[i]) * n_region + region[i]] + 
                    beta_sampled * (1-p_sampled[i]) + 
                    log(swept_area[i]));
      
      // Negative binomial distribution:
      res -= lgamma(z[i]+r) - lgamma(r) - lgamma(z[i]+1) + r*log(r) + z[i]*log(mu) - (r+z[i])*log(r+mu);;
   }
   
   // Calculate predicted length frequencies by year and region:
   array<Type> mu(n_length, n_year, n_region);
   for (int i = 0; i < n_length; i++){
      for (int j = 0; j < n_year; j++){
         for (int k = 0; k < n_region; k++){
            mu(i,j,k) = exp(alpha + 
                            length_effect[i] +
                            year_effect[j] + 
                            region_effect[k] + 
                            length_year_effect[i * n_year + j] +   
                            length_region_effect[i * n_region + k] + 
                            region_year_effect[k * n_year + j] +
                            length_year_region_effect[(i * n_year + j) * n_region + k]);
         }
      }
   }
   
   // Calculate stock indices:
   matrix<Type> abundance_L40(n_year, n_region); 
   matrix<Type> abundance_G20L40(n_year, n_region);
   matrix<Type> com(n_year, n_region);
   abundance_L40.fill(0);
   abundance_G20L40.fill(0);
   com.fill(0);
   for (int j = 0; j < n_year; j++){
      for (int k = 0; k < n_region; k++){
         for (int i = 0; i < 40; i++)     abundance_L40(j,k) += mu(i,j,k);
         for (int i = 20; i < 40; i++) abundance_G20L40(j,k) += mu(i,j,k);
         for (int i = mls(j,k); i < n_length; i++)  com(j,k) += mu(i,j,k);
      }
   }
         
   // Output:
   REPORT(mu);
   REPORT(abundance_L40);
   REPORT(abundance_G20L40);
   REPORT(com);
   
   return res;
}
