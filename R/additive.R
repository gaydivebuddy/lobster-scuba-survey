# Run 'Prepare data.R'

library(gulf.data)
library(gulf.graphics)
library(TMB)
source("U:/TMB/TMB utilities.R")

clc(); compile("R/additive.cpp")
dyn.load(dynlib("R/additive"))

# Categorical variables:
divers <- sort(unique(mm$diver))
transects <- unique(mm$transect.name)
ix <- order(unlist(lapply(strsplit(transects, "[(]"), function(x) x[2])),
            unlist(lapply(strsplit(transects, "[(]"), function(x) x[1])))
transects <- transects[ix]
regions <- sort(unique(mm$region))
years <- sort(unique(year(unique(x$sections$date))))
years <- min(years):max(years)
lens <- sort(unique(mm$carapace.length))

# Read minimum legal size table:
mls <- read.csv("data/minimum_legal_size.txt")
mls <- mls[mls$year >= 2003, ]
mls$sub.LFA[is.na(mls$sub.LFA)] <- ""
mls$sub.LFA <- substr(mls$sub.LFA,1,1)
mls <- mls[mls$LFA %in% c("23", "25", "26A"), ]
mls <- unique(mls) 
mls$size <- round(mls$size)

# For every unique sub.LFA, buffer early years with MLS from corresponding LFA:
u <- unique(mls[mls$sub.LFA != "", c("LFA", "sub.LFA")])
r <- NULL
for (i in 1:nrow(u)){
   ix <- which( (mls$LFA == u$LFA[i]) & (mls$sub.LFA == u$sub.LFA[i]))
   tmp <- data.frame(LFA = u$LFA[i], sub.LFA = u$sub.LFA[i], year = setdiff(2003:2021, mls$year[ix]))
   tmp2 <- tmp
   tmp2$sub.LFA <- ""
   iy <- match(tmp2, mls[c("LFA", "sub.LFA", "year")])
   tmp$size <- mls$size[iy]
   r <- rbind(r, tmp)
}
mls <- rbind(mls, r)
mls <- sort(mls, by = c("LFA", "sub.LFA", "year"))

# aggregate(mls$year, by = mls[c("LFA", "sub.LFA")], function(x) setdiff(2003:2021, unique(x)))

# Region - Minimum legal size index table:
tmp <- data.frame(region  = c("Caraquet", "Cocagne", "Fox Harbor", "Murray Corner", "Neguac", "Pointe-Verte", "Richibucto", "Shediac", "Toney River"),
                  LFA     = c("23", "25", "26A", "25", "23", "23", "25", "25", "26A"),
                  sub.LFA = c("C", "", "1", "", "C", "A", "", "", "1"))
tab <- expand.grid(region = regions, year = years)
ix <- match(tab$region, tmp$region)
tab$LFA <- tmp$LFA[ix]
tab$sub.LFA <- tmp$sub.LFA[ix]
ix <- match(tab[c("year", "LFA", "sub.LFA")], mls[c("year", "LFA", "sub.LFA")])
tab$size <- mls$size[ix]

# Convert to indices:
tab$region.id <- match(tab$region, regions)-1
tab$year.id <- tab$year - min(years)
tab$mls.id <- match(tab$size, lens)-1

# Square matrix of MLS values:
mls <- matrix(NA, nrow = length(regions), ncol = length(years))
dimnames(mls) <- list(region = regions, year = years)
for (i in 1:nrow(mls)){
   for (j in 1:ncol(mls)){
      ix <- (tab$region == regions[i]) & (tab$year == years[j])
      mls[i,j] <- tab$mls.id[ix]
   } 
}
mls <- t(mls)

# Prepare data:
data <- list(z = mm$n,
             len = match(mm$carapace.length, lens)-1,
             year = year(mm$date) - min(years),
             region = match(mm$region, regions)-1,
             transect = match(mm$transect.name, transects)-1, 
             diver = match(mm$diver, divers)-1, 
             p_sampled = mm$proportion.sampled,
             swept_area = mm$swept.area, 
             mls = mls)

# Initialize parameters:
parameters <- list(alpha = log(mean(data$z)),                      # Intercept parameter.
                   beta_sampled = 0,                               # Sampling proportion coefficient.
                   length_effect = rep(0, max(data$len)+1),        # Length  effect.
                   year_effect = rep(0, max(data$year)+1),         # Year effect.
                   region_effect = rep(0, max(data$region)+1),     # Region effect.
                   transect_effect = rep(0, max(data$transect)+1), # Transect effect.
                   diver_effect = rep(0, max(data$diver)+1),       # Diver effect.
                   log_sigma_length = -1,                          # Error for length effect.
                   log_sigma_year = -1,                            # Error for year effect.
                   log_sigma_region = -1,                          # Error for region effect.
                   log_sigma_transect = -1,                        # Error for transect effect.
                   log_sigma_diver = -1,                           # Error for diver effect.
                   logit_phi_year = 0,                             # Correlation para for year effect.
                   log_r = 0)                                      # Negative binomial dispersion parameter.

# Build TMB model:
random <- c("length_effect", "year_effect", "region_effect", "transect_effect", "diver_effect")
obj <- MakeADFun(data = data, parameters = parameters, random = random, DLL = "additive")

theta <- optim(obj$par, obj$fn, obj$gr, control = list(maxit = 5000, trace = 3))$par
obj$par <- theta

rep  <- sdreport(obj)
fixed <- summary(rep, "fixed")
random <- summary(rep, "random")

# Plot results:
gbarplot(exp(random[grep("length_effect", rownames(random)), 1]))
gbarplot(exp(random[grep("year_effect", rownames(random)), 1]))
gbarplot(exp(random[grep("diver_effect", rownames(random)), 1]))
hline(1, col = "red", lwd = 2)
text(1:length(divers), 0.5 * exp(random[grep("diver_effect", rownames(random)), 1]), divers, srt = 90)
gbarplot(exp(random[grep("transect_effect", rownames(random)), 1]))
gbarplot(exp(random[grep("region_effect", rownames(random)), 1]))
hline(1, col = "red", lwd = 2)
text(1:length(regions), 0.5 * exp(random[grep("region_effect", rownames(random)), 1]),regions, srt = 90)

# Compile and load interaction model:
clc(); compile("R/interaction.cpp")
dyn.load(dynlib("R/interaction"))

# Initialize parameters:
parameters <- update.parameters(parameters, fixed, random)
parameters$length_year_effect = rep(0, (max(data$len)+1) * (max(data$year)+1))         # Length x year interaction effect.
parameters$region_year_effect = rep(0, (max(data$year)+1) * (max(data$region)+1))      # Year x region interaction effect.
parameters$transect_year_effect = rep(0, (max(data$year)+1) * (max(data$transect)+1))  # Year x region interaction effect.
parameters$diver_year_effect = rep(0, (max(data$year)+1) * (max(data$diver)+1))        # Year x diver interaction effect.
parameters$length_region_effect = rep(0, (max(data$len)+1) * (max(data$region)+1))     # Length x region interaction effect.
parameters$length_diver_effect = rep(0, (max(data$len)+1) * (max(data$diver)+1))       # Length x diver interaction effect.
parameters$length_year_region_effect = rep(0, (max(data$len)+1) * (max(data$year)+1) * (max(data$region)+1))   # Length x year x region interaction effect.

parameters$log_sigma_length_year = -1         # Error for length x year interaction effect.
parameters$log_sigma_region_year = -1         # Error for year x region interaction effect.
parameters$log_sigma_transect_year = -1       # Error for year x transect interaction effect. 
parameters$log_sigma_diver_year = -1          # Error for diver x year interaction effect.
parameters$log_sigma_length_region = -1       # Error for length x region interaction effect.
parameters$log_sigma_length_diver = -1        # Error for length x diver interaction effect.
parameters$log_sigma_length_year_region = -1  # Error for length x year x region interaction effect.

parameters$logit_phi_year <- 0
parameters$logit_phi_transect_year <- 0
parameters$logit_phi_region_year <- 0
parameters$logit_phi_length_year <- 0

# Build TMB model:
map <- list(beta_sampled = factor(NA))
random <- c("length_effect", "year_effect", "region_effect", "transect_effect", "diver_effect", 
            "length_year_effect", "length_region_effect", "region_year_effect", "length_diver_effect",
            "diver_year_effect", "transect_year_effect", "length_year_region_effect")
obj <- MakeADFun(data = data, parameters = parameters, random = random, map = map, DLL = "interaction")

theta <- optim(obj$par, obj$fn, obj$gr, control = list(maxit = 5000, trace = 3))$par
obj$par <- theta

rep  <- sdreport(obj)
fixed <- summary(rep, "fixed")
random <- summary(rep, "random")

