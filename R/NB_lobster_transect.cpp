// Negative binomial lobster transect model:

#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() (){
   
   // Data variable declarations:
   DATA_VECTOR(y);          // Count observations.
   DATA_VECTOR(area);       // Area coverage.
   DATA_FACTOR(year);       // Year indicators.
   DATA_FACTOR(region);     // Region indicators.
   DATA_FACTOR(cohort);     // Cohort indicators.

   // Fixed and random effect parameters:
   PARAMETER(alpha);                            // Global intercept parameter.
   PARAMETER_VECTOR(year_effect);               // Year random effect.
   PARAMETER_VECTOR(region_effect);             // Region random effect.
   PARAMETER_VECTOR(cohort_effect);             // Cohort random effect.
   PARAMETER_VECTOR(year_region_effect);        // Year-region random effect.
   PARAMETER_VECTOR(year_cohort_effect);        // Year-cohort random effect.
   PARAMETER_VECTOR(region_cohort_effect);      // Region-cohort random effect.
   PARAMETER_VECTOR(year_region_cohort_effect); // Year-region-cohort random effect.
   
   // Log-scale error parameters:
   PARAMETER(log_sigma_year);                   // Log-scale error for year effect.
   PARAMETER(log_sigma_region);                 // Log-scale error for region effect.
   PARAMETER(log_sigma_cohort);                 // Log-scale error for cohort effect.
   PARAMETER(log_sigma_year_region);            // Log-scale error for year-region interaction effect.
   PARAMETER(log_sigma_year_cohort);            // Log-scale error for year-cohort interaction effect.
   PARAMETER(log_sigma_region_cohort);          // Log-scale error for region-cohort interaction effect.
   PARAMETER(log_sigma_year_region_cohort);     // Log-scale error for year-region-cohort interaction effect.
   PARAMETER(log_r);                            // Precision parameter for negative binomial (log-scale).

   // Transform random effect error parameters:
   Type sigma_year = exp(log_sigma_year);                             // Standard error for year effect.
   Type sigma_region = exp(log_sigma_region);                         // Standard error for region effect.
   Type sigma_cohort = exp(log_sigma_cohort);                         // Standard error for cohort effect.
   Type sigma_year_region = exp(log_sigma_year_region);               // Standard error for year-region interaction effect.
   Type sigma_year_cohort = exp(log_sigma_year_cohort);               // Standard error for year-cohort interaction effect.
   Type sigma_region_cohort = exp(log_sigma_region_cohort);           // Standard error for region-cohort interaction effect.
   Type sigma_year_region_cohort = exp(log_sigma_year_region_cohort); // Standard error for year-region-cohort interaction effect.

   // Transformed parameters:
   int n_obs = y.size();                        // Number of observations.
   int n_year = year_effect.size();             // Number of year factor levels.
   int n_region = region_effect.size();         // Number of region factor levels.
   int n_cohort = cohort_effect.size();         // Number of cohort factor levels.  
   Type r = exp(log_r);                         // Negative binomial precision parameter.
   Type mu;                                     // Regression mean variable.
   Type res = 0;                                // Log-likelihood accumulator.

   // Random effects contributions:
   res -= sum(dnorm(year_effect, Type(0), sigma_year, true));
   res -= sum(dnorm(region_effect, Type(0), sigma_region, true));
   res -= sum(dnorm(cohort_effect, Type(0), sigma_cohort, true));
   res -= sum(dnorm(year_region_effect, Type(0), sigma_year_region, true));
   res -= sum(dnorm(year_cohort_effect, Type(0), sigma_year_cohort, true));
   res -= sum(dnorm(region_cohort_effect, Type(0), sigma_region_cohort, true));
   res -= sum(dnorm(year_region_cohort_effect, Type(0), sigma_year_region_cohort, true));
               
   // Observation log-likelihood:
   for(int i = 0; i < n_obs; i++){
      // Define log-linear mean:
      mu = exp(alpha +
               year_effect[year[i]-1] +
               region_effect[region[i]-1] +
               cohort_effect[cohort[i]-1] +
               year_region_effect[(year[i]-1) * n_region + (region[i]-1)] +
               year_cohort_effect[(year[i]-1) * n_cohort + (cohort[i]-1)] +
               region_cohort_effect[(region[i]-1) * n_cohort + (cohort[i]-1)] +
               year_region_cohort_effect[((year[i]-1) * n_region + (region[i]-1)) * n_cohort + (cohort[i]-1)] +
               log(area[i]) - log(100));

      // Negative binomial log-density:
      res -= lgamma(y[i]+r) - lgamma(r) - lgamma(y[i]+1) + r*log(r) + y[i]*log(mu) - (r+y[i])*log(r+mu);
   }

   return res;
}
