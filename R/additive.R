# Run 'Prepare data.R'

library(gulf.data)
library(gulf.graphics)
library(TMB)
source("U:/TMB/TMB utilities.R")

clc(); compile("additive.cpp")
dyn.load(dynlib("additive"))

# Categorical variables:
divers <- sort(unique(mm$diver))
transects <- unique(mm$transect.name)
ix <- order(unlist(lapply(strsplit(transects, "[(]"), function(x) x[2])),
            unlist(lapply(strsplit(transects, "[(]"), function(x) x[1])))
transects <- transects[ix]
regions <- sort(unique(mm$region))
lens <- sort(unique(mm$carapace.length))

# Prepare data:
data <- list(z = mm$n,
             len = match(mm$carapace.length, sort(unique(mm$carapace.length)))-1,
             year = year(mm$date) - min(year(mm$date)),
             region = match(mm$region, regions)-1,
             transect = match(mm$transect.name, transects)-1, 
             diver = match(mm$diver, divers)-1, 
             swept_area = mm$swept.area)

# Initialize parameters:
parameters <- list(alpha = log(mean(data$z)),                      # Intercept parameter.
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
                   log_r = 0)                                      # Negative binomial dispersion parameter.

# Build TMB model:
random <- c("length_effect", "year_effect", "region_effect", "transect_effect", "diver_effect")
obj <- MakeADFun(data = data, parameters = parameters, random = random, DLL = "additive")

theta <- optim(obj$par, obj$fn, obj$gr, control = list(maxit = 5000, trace = 3))
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
clc(); compile("interaction.cpp")
dyn.load(dynlib("interaction"))

# Initialize parameters:
parameters <- update.parameters(parameters, fixed, random)
parameters$logit_phi_year <- 0
parameters$logit_phi_year_transect <- 0
parameters$logit_phi_year_region <- 0

parameters$length_year_effect = rep(0, (max(data$len)+1) * (max(data$year)+1))     # Length x year interaction effect.
parameters$length_region_effect = rep(0, (max(data$len)+1) * (max(data$region)+1)) # Length x region interaction effect.
parameters$year_region_effect = rep(0, (max(data$year)+1) * (max(data$region)+1))  # Year x region interaction effect.
parameters$length_diver_effect = rep(0, (max(data$len)+1) * (max(data$diver)+1))   # Length x diver interaction effect.
parameters$diver_year_effect = rep(0, (max(data$year)+1) * (max(data$diver)+1))    # Year x diver interaction effect.
parameters$year_transect_effect = rep(0, (max(data$year)+1) * (max(data$transect)+1))  # Year x region interaction effect.
parameters$length_year_region_effect = rep(0, (max(data$len)+1) * (max(data$year)+1) * (max(data$region)+1))   # Length x year x region interaction effect.

parameters$log_sigma_length_year = -1         # Error for length x year interaction effect.
parameters$log_sigma_length_region = -1       # Error for length x region interaction effect.
parameters$log_sigma_year_region = -1         # Error for year x region interaction effect.
parameters$log_sigma_length_diver = -1        # Error for length x diver interaction effect.
parameters$log_sigma_diver_year = -1          # Error for diver x year interaction effect.
parameters$log_sigma_year_transect <- -1      # Error for year x transect interaction effect. 
parameters$log_sigma_length_year_region = -1  # Error for length x year x region interaction effect.

# Build TMB model:
random <- c("length_effect", "year_effect", "region_effect", "transect_effect", "diver_effect", 
            "length_year_effect", "length_region_effect", "year_region_effect", "length_diver_effect",
            "diver_year_effect", "year_transect_effect", "length_year_region_effect")
obj <- MakeADFun(data = data, parameters = parameters, random = random, DLL = "interaction")

theta <- optim(obj$par, obj$fn, obj$gr, control = list(maxit = 5000, trace = 3))
obj$par <- theta

rep  <- sdreport(obj)
fixed <- summary(rep, "fixed")
random <- summary(rep, "random")

# Plot results:
gbarplot(exp(random[rownames(random) == "length_effect", 1]))
gbarplot(exp(random[rownames(random) == "year_effect", 1]))
gbarplot(exp(random[rownames(random) == "diver_effect", 1]))
hline(1, col = "red", lwd = 2)
text(1:length(divers), 0.5 * exp(random[rownames(random) == "diver_effect", 1]), divers, srt = 90)

diver_effect <- exp(random[rownames(random) == "diver_effect", 1])
names(diver_effect) <- divers
diver_effect <- diver_effect[setdiff(names(diver_effect), c("Asselin, Natalie", "Leblanc, Stepan"))]
gbarplot(diver_effect, grid = TRUE)
hline(1, col = "red", lwd = 2)
text(1:length(diver_effect), 0.5 * diver_effect, names(diver_effect), srt = 90, cex = 1.25)
mtext("Diver", 1, 2.75, cex = 1.5)
mtext("Diver effect", 2, 2.75, cex = 1.5)



gbarplot(exp(random[rownames(random) == "transect_effect", 1]))
gbarplot(exp(random[rownames(random) == "region_effect", 1]))
hline(1, col = "red", lwd = 2)
text(1:length(regions), 0.5 * exp(random[rownames(random) == "region_effect", 1]),regions, srt = 90)


plot(random[rownames(random) == "year_region_effect", 1])


gbarplot(exp(random[rownames(random) == "length_diver_effect", 1]))

# Length x diver effect effect:
length_diver <- random[rownames(random) == "length_diver_effect", 1]
dim(length_diver) <- c(max(data$diver)+1, max(data$len)+1)
dimnames(length_diver) <- list(divers, lens)
length_diver <- length_diver[setdiff(rownames(length_diver), c("Asselin, Natalie", "Leblanc, Stepan")), ]
image(1:nrow(length_diver), 1:ncol(length_diver), length_diver, col = colorRampPalette(c("blue", "white", "red"))(10),
      xaxt = "n", yaxt = "n", zlim = c(-0.5, 0.5))
axis(1, at = 1:nrow(length_diver), labels = rownames(length_diver), las = 2)
axis(2, at = 1:ncol(length_diver), labels = colnames(length_diver), las = 2)

# Diver x year effect:
diver_year_effect <- random[rownames(random) == "diver_year_effect", 1]
dim(diver_year_effect) <- c(max(data$year)+1, max(data$diver)+1)
dimnames(diver_year_effect) <- list(min(year(mm$date)):max(year(mm$date)), divers)
m <- cbind(matrix(0, nrow = 7, ncol = 1), matrix(1, nrow = 7, ncol = 7))
m <- cbind(0, rbind(0, m, 0), 0)
layout(m)
par(mar = c(0,0,0,0))
image(1:nrow(diver_year_effect), 1:ncol(diver_year_effect), diver_year_effect, 
      col = colorRampPalette(c("blue", "white", "red"))(10),
      xaxt = "n", yaxt = "n", zlim = c(-0.1, 0.1))
axis(1, at = 1:nrow(diver_year_effect), labels = rownames(diver_year_effect), las = 2)
axis(2, at = 1:ncol(diver_year_effect), labels = colnames(diver_year_effect), las = 2)

# Length x year effect:
length_year <- random[rownames(random) == "length_year_effect", 1]
dim(length_year) <- c(max(data$len)+1, max(data$year)+1)
dimnames(length_year) <- list(lens, min(year(mm$date)):max(year(mm$date)))
image(1:nrow(length_year), 1:ncol(length_year), length_year, col = colorRampPalette(c("blue", "white", "red"))(10),
      xaxt = "n", yaxt = "n", zlim = c(-0.65, 0.65))
axis(1, at = 1:nrow(length_year), labels = rownames(length_year), las = 2)
axis(2, at = 1:ncol(length_year), labels = colnames(length_year), las = 2)
box(col = "grey50")

# Length x region effect:
length_region <- random[rownames(random) == "length_region_effect", 1]
dim(length_region) <- c(max(data$len)+1, max(data$region)+1)
dimnames(length_region) <- list(lens, regions)
m <- cbind(matrix(0, nrow = 7, ncol = 1), matrix(1, nrow = 7, ncol = 7))
m <- cbind(0, rbind(0, m, 0), 0)
layout(m)
par(mar = c(0,0,0,0))
image(1:nrow(length_region), 1:ncol(length_region), length_region, col = colorRampPalette(c("blue", "white", "red"))(10),
      xaxt = "n", yaxt = "n", zlim = c(-1.5, 1.5))
axis(1, at = 1:nrow(length_region), labels = rownames(length_region), las = 2)
axis(2, at = 1:ncol(length_region), labels = colnames(length_region), las = 2)
box(col = "grey50")

# Year x region effect:
year_region <- random[rownames(random) == "year_region_effect", 1]
dim(year_region) <- c(max(data$year)+1, max(data$region)+1)
dimnames(year_region) <- list(min(year(mm$date)):max(year(mm$date)), regions)
m <- cbind(matrix(0, nrow = 7, ncol = 1), matrix(1, nrow = 7, ncol = 7))
m <- cbind(0, rbind(0, m, 0), 0)
layout(m)
par(mar = c(0,0,0,0))
image(1:nrow(year_region), 1:ncol(year_region), year_region, col = colorRampPalette(c("blue", "white", "red"))(20),
      xaxt = "n", yaxt = "n", zlim = c(-0.85, 0.85))
axis(1, at = 1:nrow(year_region), labels = rownames(year_region), las = 2)
axis(2, at = 1:ncol(year_region), labels = colnames(year_region), las = 2)
box(col = "grey50")

# Predictions:
years <- sort(unique(year(mm$date)))
mu <- obj$report()$mu
dimnames(mu) <- list(length = lens, 
                     year = min(unique(year(mm$date))):max(unique(year(mm$date))),
                     region = regions)

years <- 2010:2021
m <- kronecker(matrix(1:12, ncol = 2), matrix(1, nrow = 5, ncol = 5))
m <- cbind(0, rbind(0, m, 0, 0), 0)
layout(m)
par(mar = c(0,0,0,0))
for (i in 1:length(years)){
   gbarplot(100*mu[,as.character(years[i]), "Caraquet"], 
            xlim = c(0, 100), ylim = c(0, 35), xaxs = "i", yaxs = "i", 
            xaxt = "n", yaxt = "n", grid = TRUE)
   text(par("usr")[1] + 0.8 * diff(par("usr")[1:2]),
        par("usr")[3] + 0.8 * diff(par("usr")[3:4]),
        years[i], cex = 1.25)
   box(col = "grey60")
   if (i %in% 1:6) axis(2)
   if (i %in% c(6, 12)) axis(1)
}

# Variance parameters:
vars <- theta$par[grep("log_sigma_", names(theta$par))]
gbarplot(exp(vars), grid = TRUE)
text(1:length(vars), 0.5 * exp(vars), gsub("log_sigma_", "", names(vars)), srt = 90)



