// Negative binomial lobster transect model:

#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
   // Data variable declarations:
   DATA_VECTOR(y);        // Count observations.
   DATA_VECTOR(area);     // Area coverage.
   DATA_FACTOR(year);     // Year indicators.
   DATA_FACTOR(region);   // Region indicators.
   DATA_FACTOR(cohort);   // Cohort indicators.
   DATA_FACTOR(transect); // Transect indicators.
   
   // Parameter and random effect declarations:
   PARAMETER(alpha);                            // Global intercept parameter.
   PARAMETER_VECTOR(year_effect);               // Year random effect.
   PARAMETER_VECTOR(region_effect);             // Region random effect.
   PARAMETER_VECTOR(cohort_effect);             // Cohort random effect.
   PARAMETER_VECTOR(year_region_effect);        // Year-region random effect.
   PARAMETER_VECTOR(year_cohort_effect);        // Year-cohort random effect.
   PARAMETER_VECTOR(region_cohort_effect);      // Region-cohort random effect.
   PARAMETER_VECTOR(year_region_cohort_effect); // Year-region-cohort random effect.
   PARAMETER_VECTOR(transect_effect);           // Transect random effect.

   // Random effect error parameters:
   PARAMETER(log_sigma_year);                   // Log-scale error for year effect.
   PARAMETER(log_sigma_region);                 // Log-scale error for region effect.
   PARAMETER(log_sigma_cohort);                 // Log-scale error for cohort effect.
   PARAMETER(log_sigma_year_region);            // Log-scale error for year-region interaction effect.
   PARAMETER(log_sigma_year_cohort);            // Log-scale error for year-cohort interaction effect.
   PARAMETER(log_sigma_region_cohort);          // Log-scale error for region-cohort interaction effect.
   PARAMETER(log_sigma_year_region_cohort);     // Log-scale error for year-region-cohort interaction effect.
   PARAMETER(log_sigma_transect);               // Log-scale error for transect effect.
   PARAMETER(log_r);                            // Precision parameter for negative binomial (log-scale).

   // Transform random effect error parameters:
   Type sigma_year               = exp(log_sigma_year);               // Standard error for year effect.
   Type sigma_region             = exp(log_sigma_region);             // Standard error for region effect.
   Type sigma_cohort             = exp(log_sigma_cohort);             // Standard error for cohort effect.
   Type sigma_year_region        = exp(log_sigma_year_region);        // Standard error for year-region interaction effect.
   Type sigma_year_cohort        = exp(log_sigma_year_cohort);        // Standard error for year-cohort interaction effect.
   Type sigma_region_cohort      = exp(log_sigma_region_cohort);      // Standard error for region-cohort interaction effect.
   Type sigma_year_region_cohort = exp(log_sigma_year_region_cohort); // Standard error for year-region-cohort interaction effect.
   Type sigma_transect           = exp(log_sigma_transect);           // Standard error for transect effect.

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
   res -= log(2 * 5) - log(pi) - log(5*5 + (sigma_year_region * sigma_year_region));
   res -= log(2 * 5) - log(pi) - log(5*5 + (sigma_year_cohort * sigma_year_cohort));
   res -= log(2 * 5) - log(pi) - log(5*5 + (sigma_region_cohort * sigma_region_cohort));
   res -= log(2 * 5) - log(pi) - log(5*5 + (sigma_year_region_cohort * sigma_year_region_cohort));
   res -= log(2 * 5) - log(pi) - log(5*5 + (sigma_transect * sigma_transect));

   // Hierarchical prior contributions to log-likelihood:
   res -= sum(dnorm(year_effect, Type(0), sigma_year, true));
   res -= sum(dnorm(region_effect, Type(0), sigma_region, true));
   res -= sum(dnorm(cohort_effect, Type(0), sigma_cohort, true));
   res -= sum(dnorm(year_region_effect, Type(0), sigma_year_region, true));
   res -= sum(dnorm(year_cohort_effect, Type(0), sigma_year_cohort, true));
   res -= sum(dnorm(region_cohort_effect, Type(0), sigma_region_cohort, true));
   res -= sum(dnorm(year_region_cohort_effect, Type(0), sigma_year_region_cohort, true));
   res -= sum(dnorm(transect_effect, Type(0), sigma_transect, true));
                  
   // Compute categorical mean values:
   array<Type> log_mu(n_year, n_region, n_cohort);
   for(int i = 0; i < n_year; i++){
      for(int j = 0; j < n_region; j++){
         for(int k = 0; k < n_cohort; k++){
            log_mu(i,j,k) = alpha +
                            year_effect[i] +
                            region_effect[j] +
                            cohort_effect[k] +
                            year_region_effect[i * n_region + j] +
                            year_cohort_effect[i * n_cohort + k] +
                            region_cohort_effect[j * n_cohort + k] +
                            year_region_cohort_effect[(i * n_region + j) * n_cohort + k];                     
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
   
   // Compute year and region marginal means:
   array<Type> mu_year_region_cohort(n_year, n_region, n_cohort);
   for (int i = 0; i < n_year; i++){
      for (int j = 0; j < n_region; j++){
         for (int k = 0; k < n_cohort; k++){
            mu_year_region_cohort(i,j,k) = exp(log_mu(i,j,k));                  
         }
      }
   }
      
   // Compute year and region marginal means:
   array<Type> mu_year_region(n_year, n_region);
   for (int i = 0; i < n_year; i++){
      for (int j = 0; j < n_region; j++){
         mu_year_region(i,j) = 0;
         for (int k = 0; k < n_cohort; k++){
            mu_year_region(i,j) += exp(log_mu(i,j,k));                  
         }
         mu_year_region(i,j) = mu_year_region(i,j) / n_cohort;
      }
   }
   
   // Compute year and cohort marginal means:
   array<Type> mu_year_cohort(n_year, n_cohort);
   for (int i = 0; i < n_year; i++){
      for (int k = 0; k < n_cohort; k++){
         mu_year_cohort(i,k) = 0;
         for (int j = 0; j < n_region; j++){
            mu_year_cohort(i,k) += exp(log_mu(i,j,k));                  
         }
         mu_year_cohort(i,k) = mu_year_cohort(i,k) / n_region;
      }
   }
   
   // Compute region and cohort marginal means:
   array<Type> mu_region_cohort(n_region, n_cohort);
   for (int j = 0; j < n_region; j++){
      for (int k = 0; k < n_cohort; k++){
         mu_region_cohort(j,k) = 0;
         for (int i = 0; i < n_year; i++){
            mu_region_cohort(j,k) += exp(log_mu(i,j,k));                  
         }
         mu_region_cohort(j,k) = mu_region_cohort(j,k) / n_year;
      }
   }
      
   // Compute year marginal means:
   array<Type> mu_year(n_year);
   for (int i = 0; i < n_year; i++){
      mu_year(i) = 0;
      for (int j = 0; j < n_region; j++){
         mu_year(i) += mu_year_region(i,j);                  
      }
      mu_year(i) = mu_year(i) / n_region;
   }
   
   // Compute region marginal means:
   array<Type> mu_region(n_region);
   for (int j = 0; j < n_region; j++){
      mu_region(j) = 0;
      for (int i = 0; i < n_year; i++){
         mu_region(j) += mu_year_region(i,j);                  
      }
      mu_region(j) = mu_region(j) / n_year;
   }
   
   // Compute region marginal means (2000-2012) and (2013-2016):
   array<Type> mu_region_early(n_region);
   for (int j = 0; j < n_region; j++){
      mu_region_early(j) = 0;
      for (int i = 0; i < (n_year-4); i++){
         mu_region_early(j) += mu_year_region(i,j);                  
      }
      mu_region_early(j) = mu_region_early(j) / (n_year-Type(4));
   }
   array<Type> mu_region_late(n_region);
   for (int j = 0; j < n_region; j++){
      mu_region_late(j) = 0;
      for (int i = (n_year-4); i < n_year; i++){
         mu_region_late(j) += mu_year_region(i,j);                  
      }
      mu_region_late(j) = mu_region_late(j) / (Type(4));
   }
    
   // Compute cohort marginal means:
   array<Type> mu_cohort(n_cohort);
   for(int k = 0; k < n_cohort; k++){
      mu_cohort(k) = 0;
      for(int i = 0; i < n_year; i++){
         mu_cohort(k) += mu_year_cohort(i,k);                  
      }
      mu_cohort(k) = mu_cohort(k) / n_year;
   }
   
   // Compute cohort ratio values:
   //array<Type> R(n_year-1, n_region, n_cohort-1);
   //for(int i = 1; i < n_year; i++){
   //   for(int j = 0; j < n_region; j++){
   //      for(int k = 1; k < n_cohort; k++){
   //         R(i-1,j,k-1) = exp(log_mu(i,j,k) - log_mu(i-1,j,k-1));
   //      }
   //   }
   //}
   
   // Compute marginal cohort ratio values:
   // array<Type> R_cohort(n_cohort-1);
   // for(int k = 1; k < n_cohort; k++){
   //   R_cohort(k-1) = mu_cohort(k) / mu_cohort(k-1);
   // }
   
   // R Output:
   ADREPORT(mu_year);
   ADREPORT(mu_cohort);
   ADREPORT(mu_region);
   ADREPORT(mu_region_early);
   ADREPORT(mu_region_late);
   ADREPORT(mu_year_region);
   ADREPORT(mu_year_region_cohort);
   
   //ADREPORT(R_year);   
   //ADREPORT(R_region); 
   //ADREPORT(R_cohort); 
   //ADREPORT(R_year_region);
   //ADREPORT(R_year_cohort);
   //ADREPORT(R_region_cohort);
   
   return res;
}
