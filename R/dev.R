library(gulf.data)
library(TMB)

setwd("U:/Lobster/Transect Analysis/TMB")
source("U:/TMB/TMB utilities.R")

compile("dev.cpp")
dyn.load(dynlib("dev"))

# Target subset of data:
# Read SCUBA lobster data:
x <- read.table("C:/Users/SuretteTJ/Desktop/github/lobster-scuba-survey/data/raw/Scuba_Master.csv", header = TRUE, sep = ",", colClasses = "character", stringsAsFactors = FALSE)
x$length <- as.numeric(x$LC_.mm.)
x <- x[!is.na(x$length), ]
x$year <- x$"Année"
names(x) <- tolower(names(x))

# Input data:
data <- aggregate(list(y = x$length), by = x[c("year", "site", "transect")], length)
data$year <- as.numeric(data$year)
data$year <- data$year - min(data$year)
   
# Initialize full parameter vector:
parameters <- list(alpha = 0, 
                   year_effect = rep(0, length(unique(data$year))), 
                   log_sigma_year = 0,
                   log_r = 0)

# Create TMB object:
obj <- MakeADFun(data = data[c("y", "year")], 
                 parameters = parameters,
                 random = "year_effect",
                 DLL = "dev")

theta <- optim(obj$par, obj$fn, obj$gr, control = list(maxit = 5000, trace = 3))
obj$par <- theta$par

rep  <- sdreport(obj)
fixed <- summary(rep, "fixed")
random <- summary(rep, "random")
random <- extract.tmb(random)

dbarplot(exp(fixed["alpha", 1] + random$year_effect), sort(unique(res$year)))
dbarplot(exp(fixed["alpha", 1] + random$region_effect))
dbarplot(exp(fixed["alpha", 1] + random$cohort_effect))

# Create TMB object:
for (i in 1:length(random)){
   if (names(random)[i] %in% names(parameters)){
      parameters[[names(random)[i]]] <- random[[i]]
   }
}
for (i in 1:nrow(fixed)){
   if (rownames(fixed)[i] %in% names(parameters)){
      parameters[[rownames(fixed)[i]]] <- as.numeric(fixed[i,1])
   }
}
obj <- MakeADFun(data = data[data.cpp("NB_lobster_transect_interaction.cpp")], 
                 parameters = parameters[parameters.cpp("NB_lobster_transect_interaction.cpp")],
                 random = c("year_effect", "region_effect", "cohort_effect", "transect_effect", "year_region_effect", "year_cohort_effect", "region_cohort_effect"),
                 DLL = "NB_lobster_transect_interaction")

theta <- optim(obj$par, obj$fn, obj$gr, control = list(maxit = 5000, trace = 3))
theta <- optim(theta$par, obj$fn, obj$gr, control = list(maxit = 5000, trace = 3))

obj$par <- theta$par
rep  <- sdreport(obj)
fixed <- summary(rep, "fixed")
random <- summary(rep, "random")
random <- extract.tmb(random)

dbarplot(random$year_effect)
dbarplot(random$region_effect)
dbarplot(random$cohort_effect)

for (i in 1:length(random)){
   if (names(random)[i] %in% names(parameters)){
      parameters[[names(random)[i]]] <- as.numeric(random[[i]])
   }
}
for (i in 1:nrow(fixed)){
   if (rownames(fixed)[i] %in% names(parameters)){
      parameters[[rownames(fixed)[i]]] <- as.numeric(fixed[i,1])
   }
}

obj <- MakeADFun(data = data[data.cpp("NB_lobster_transect_full.cpp")], 
                 parameters = parameters[parameters.cpp("NB_lobster_transect_full.cpp")],
                 random = c("year_effect", "region_effect", "cohort_effect", "year_region_effect", "year_cohort_effect", 
                            "region_cohort_effect", "year_region_cohort_effect", "transect_effect"),
                 DLL = "NB_lobster_transect_full")

theta <- optim(obj$par, obj$fn, obj$gr, control = list(maxit = 5000, trace = 3))
theta <- optim(theta$par, obj$fn, obj$gr, control = list(maxit = 5000, trace = 3))

obj$par <- theta$par
rep  <- sdreport(obj)
fixed <- summary(rep, "fixed")
random <- summary(rep, "random")
random <- extract.tmb(random)

n <- c(length(unique(data$year)), length(unique(data$region)), length(unique(data$cohort)))
muhat <- rep$value[names(rep$value) == "mu_year_region_cohort"]
dim(muhat) <- n
muhat.sd <- rep$sd[names(rep$value) == "mu_year_region_cohort"]
dim(muhat.sd) <- n

obj$par <- theta$par
rep  <- sdreport(obj, getReportCovariance = FALSE)
fixed <- summary(rep, "fixed")
random <- summary(rep, "random")
ranef <- extract.tmb(random[rownames(random) != "year_region_cohort_effect", ])
ranef.sd <- extract.tmb(random[rownames(random) != "year_region_cohort_effect", ], "sd")

# Year effect:
years <- sort(unique(res$year))
mu_year <- 4 * rep$value[names(rep$value) == "mu_year"]
sigma_year <- 4 * rep$sd[names(rep$value) == "mu_year"]
dbarplot(mu_year, years, xlab = "", ylab = "", ylim = c(0, 70), xaxt = "n")
grid()
dbarplot(mu_year, years, add = TRUE)
for (i in 1:length(years)){
   lines(c(years[i], years[i]), c(mu_year[i]-1.96*sigma_year[i], mu_year[i]+1.96*sigma_year[i]))
   lines(c(years[i]-0.15, years[i]+0.15), c(mu_year[i]-1.96*sigma_year[i], mu_year[i]-1.96*sigma_year[i]))
   lines(c(years[i]-0.15, years[i]+0.15), c(mu_year[i]+1.96*sigma_year[i], mu_year[i]+1.96*sigma_year[i]))
} 
axis(1, at = seq(2003, max(years), by = 2), cex.axis = 1, las = 2)
axis(1, at = seq(2004, max(years), by = 2), cex.axis = 1, las = 2)
box()
mtext("Year", 1, 3.5, cex = 1.4)
mtext(expression(paste("Abundance (number/400 ", m^2, ")")), 2, 2.3, cex = 1.4)

cbind(years, mu_year, mu_year-1.96*sigma_year, mu_year+1.96*sigma_year)

# Region effect:
#regions <- c("PV", "GA1", "AB1", "CR", "PH", "TB", "NG", "RI", "CG", "SH", "RO", "MC", "FH", "TR")
mu_region <- 4*rep$value[names(rep$value) == "mu_region"]
sigma_region <- 4*rep$sd[names(rep$value) == "mu_region"]
par(mar = c(10, 4, 4, 2) + 0.1) # c(bottom, left, top, right)
dbarplot(mu_region, regions, xlab = "", ylab = expression(paste("Abundance (number/400 ", m^2, ")")), ylim = c(0, 45), xaxt = "n")
old <- data.frame(region = c("PV", "CR", "NG", "RI", "CG", "SH", "MC", "FH", "TR"),
                  label = c("Pointe-Verte (23BC)", "Caraquet (23BC)", "Neguac (23G)", "Richibuctou (25N)", "Cocagne (25S)", "Shediac (25S)", 
                            "Murray Corner (25S)", "Fox Harbour (26AD)", "Toney River (26AD)"),
                  mu = c(2.1, 2.3, 4.2, 3.2, 1.3, 0.7, 0.7, 0.1, 0.9),
                  lower = c(1.5, 2.1, 2.7, 2.2, 1.0, 0.5, 0.3,	0.1, 0.5),		
                  upper = c(3.0, 2.7, 6.6, 4.9, 1.7, 0.8, 1.3, 0.2, 1.4))                 
grid()
dbarplot(mu_region, regions, add = TRUE)
# dbarplot(4*old$mu, regions, add = TRUE, col = "red")
lines(match(old$region, regions), 4*old$mu, col = "red", lwd = 2)
for (i in 1:length(regions)){
   lines(c(i, i), c(mu_region[i]-1.96*sigma_region[i], mu_region[i]+1.96*sigma_region[i]))
   lines(c(i-0.15, i+0.15), c(mu_region[i]-1.96*sigma_region[i], mu_region[i]-1.96*sigma_region[i]))
   lines(c(i-0.15, i+0.15), c(mu_region[i]+1.96*sigma_region[i], mu_region[i]+1.96*sigma_region[i]))
   
   lines(c(match(old$region, regions)[i]-0.15, match(old$region, regions)[i]+0.15), c(4*old$lower[i], 4*old$lower[i]), col = "red", lwd = 1)
   lines(c(match(old$region, regions)[i]-0.15, match(old$region, regions)[i]+0.15), c(4*old$upper[i], 4*old$upper[i]), col = "red", lwd = 1)
   lines(c(match(old$region, regions)[i], match(old$region, regions)[i]), c(4*old$lower[i], 4*old$upper[i]), col = "red", lwd = 1)
} 
box() 
axis(1, at = 1:length(regions), label = old$label, cex.axis = 0.90, las = 2)
mtext("Region", 1, 7)
legend("topright", legend = c("2003-2012", "2003-2018"), pch = c(NA, 22), lwd = c(2, NA), col = c("red", "black"), pt.bg = "grey", pt.cex = c(1, 3), bg = "white")

data.frame(regions, mu_region, mu_region-1.96*sigma_region, mu_region+1.96*sigma_region) 

# Cohort effect:
mu_cohort <- 1 * rep$value[names(rep$value) == "mu_cohort"]
sigma_cohort <- 1 * rep$sd[names(rep$value) == "mu_cohort"]
cohorts <- 0:6
dbarplot(mu_cohort, cohorts, xlab = "Cohort", ylab = expression(paste("Abundance (number/400 ", m^2, ")")), ylim = c(0, 12))
grid()
dbarplot(mu_cohort, cohorts, add = TRUE)
for (i in 1:length(cohorts)){
   lines(c(cohorts[i], cohorts[i]), c(mu_cohort[i]-1.96*sigma_cohort[i], mu_cohort[i]+1.96*sigma_cohort[i]))
   lines(c(cohorts[i]-0.15, cohorts[i]+0.15), c(mu_cohort[i]-1.96*sigma_cohort[i], mu_cohort[i]-1.96*sigma_cohort[i]))
   lines(c(cohorts[i]-0.15, cohorts[i]+0.15), c(mu_cohort[i]+1.96*sigma_cohort[i], mu_cohort[i]+1.96*sigma_cohort[i]))
} 
box()

data.frame(cohorts, mu_cohort, mu_cohort-1.96*sigma_cohort, mu_cohort+1.96*sigma_cohort) 


# Year-Region effect:
mu_year_region <- 4*rep$value[names(rep$value) == "mu_year_region"]
sigma_year_region <- 4*rep$sd[names(rep$value) == "mu_year_region"]
temp_mu <- matrix(NA, nrow = length(years), ncol = length(regions))
dimnames(temp_mu) <- list(year = years, region = regions)
temp_sd <- temp_mu
for (i in 1:length(years)){
   for (j in 1:length(regions)){
      temp_mu[i,j] <- mu_year_region[(j-1) * length(years) + i]
      temp_sd[i,j] <- sigma_year_region[(j-1) * length(years) + i]
   }
}
mu_year_region <- temp_mu
sigma_year_region <- temp_sd
plot(c(2003, max(years)), c(0, 140), type = "n", xlab = "Year", ylab = expression(paste("Abundance (n / 400 ", m^2, ")")), yaxs = "i")
grid()
vars <- c("CR", "SH", "FH")
pch = c(21, 22, 24)
col = c("white", "grey", "black")
for (j in 1:length(vars)){
   lines(years, mu_year_region[, vars[j]])
   for (i in 1:length(years)){
      lines(c(years[i], years[i]), c(mu_year_region[i, vars[j]]-1.96*sigma_year_region[i, vars[j]], mu_year_region[i, vars[j]]+1.96*sigma_year_region[i, vars[j]]))
      lines(c(years[i]-0.15, years[i]+0.15), c(mu_year_region[i, vars[j]]-1.96*sigma_year_region[i, vars[j]], mu_year_region[i, vars[j]]-1.96*sigma_year_region[i, vars[j]]))
      lines(c(years[i]-0.15, years[i]+0.15), c(mu_year_region[i, vars[j]]+1.96*sigma_year_region[i, vars[j]], mu_year_region[i, vars[j]]+1.96*sigma_year_region[i, vars[j]]))
   } 
   points(years, mu_year_region[, vars[j]], pch = pch[j], bg = col[j])
}
box() 

# Parse year x region x cohort variable:
cohort <- "1"
mu_year_region_cohort <- 4 * rep$value[names(rep$value) == "mu_year_region_cohort"]
sigma_year_region_cohort <- 4 * rep$sd[names(rep$value) == "mu_year_region_cohort"]
cohorts <- 0:6
temp_mu <- array(NA, dim = c(length(years), length(regions), length(cohorts)))
dimnames(temp_mu) <- list(year = years, region = regions, cohort = cohorts)
temp_sd <- temp_mu
for (i in 1:length(years)){
   for (j in 1:length(regions)){
      for (k in 1:length(cohorts)){
         # year x regions: (j-1) * length(years) + i
         # year x region x cohort : (k-1)*((j-1) * length(years) + i) + k 
         index <- (k-1) * length(years) * length(regions) + (j-1) * length(years) + i
         #index <- (k-1)*((j-1) * length(years) + i-1) + k
         temp_mu[i,j,k] <- mu_year_region_cohort[index]
         temp_sd[i,j,k] <- sigma_year_region_cohort[index]
      }
   }
}
mu_year_region_cohort <- temp_mu
sigma_year_region_cohort <- temp_sd
plot(c(2003, max(years)), c(0,  1.2 * max(mu_year_region_cohort[,vars, cohort])), type = "n", xlab = "Year", ylab = expression(paste("Abundance (n / 400 ", m^2, ")")), yaxs = "i")
grid()
vars <- c("CR", "CG", "SH", "MC", "FH")
pch = c(21, 22, 23, 24, 25)
col = c("white", "grey20", "grey50", "grey70", "black")
col = rainbow(length(vars))
for (j in 1:length(vars)){
   lines(years, mu_year_region_cohort[, vars[j], cohort], col = "grey50", lwd = 4)
   lines(years, mu_year_region_cohort[, vars[j], cohort], col = col[j], lwd = 2)
   for (i in 1:length(years)){
      m <- mu_year_region_cohort[as.character(years[i]), vars[j], cohort]
      s <- sigma_year_region_cohort[as.character(years[i]), vars[j], cohort]
      lines(c(years[i], years[i]), c(m-1.96*s, m+1.96*s), col = col[j])
      lines(c(years[i]-0.15, years[i]+0.15), c(m-1.96*s, m-1.96*s), col = col[j])
      lines(c(years[i]-0.15, years[i]+0.15), c(m+1.96*s, m+1.96*s), col = col[j])
   } 
   points(years,  mu_year_region_cohort[, vars[j], cohort], pch = pch[j], bg = col[j], cex = 1.25)
}
for (j in 1:length(vars)){
   points(years,  mu_year_region_cohort[, vars[j], cohort], pch = pch[j], bg = col[j], cex = 1.25)
}
box() 
title(main = cohort)
legend("topleft", legend = c("Caraquet", "Cocagne", "Shediac", "Murray Corner", "Fox Harbour"), pch = 21:25, pt.cex = 1.5, pt.bg = col, bg = "white")

res <- data.frame()
for (j in 1:length(vars)){
   m <- mu_year_region_cohort[as.character(years), vars[j], "1"]
   s <- sigma_year_region_cohort[as.character(years), vars[j], "1"]
   res <- rbind(res, data.frame(vars[j], years, m, m-1.96*s, m+1.96*s))
}   


# Bubble plot function for showing multi-way effects:
bubble.plot <- function(v, scale = 2){
   plot(range(1:nrow(v)), range(1:ncol(v)), xlab = names(dimnames(v))[1], ylab = names(dimnames(v))[2], type = "n", xaxt = "n", yaxt = "n")
   for (i in 1:nrow(v)){
      for (j in 1:ncol(v)){
         if (v[i,j] > 0) col <- "white" else col <- "grey30"
         points(i, j, cex = scale * sqrt(abs(v[i,j])) / max(abs(v)), pch = 21, bg = col)
         #if (v[i,j] == 0) points(i, j, pch = "x", lwd = 2, cex = 0.5)
      }
   }
   axis(1, at = 1:nrow(v), labels = dimnames(v)[[1]])
   axis(2, at = 1:ncol(v), labels = dimnames(v)[[2]]) 
}


# Parse 3-way effect:
year_region_cohort_effect <- random[rownames(random) == "year_region_cohort_effect", 1]

v <- array(NA, length(random$year_region_cohort_effect))
dim(v) <- c(length(years), length(regions), length(cohorts))

for (i in 1:length(years)){
   for (j in 1:length(regions)){
      for (k in 1:length(cohorts)){
         ii <- ((i-1) * length(regions) + j-1) * length(cohorts) + k
         print(ii)
         v[i,j,k] <- random$year_region_cohort_effect[ii]
      }
   }
}
dimnames(v) <- list(year = years, region = regions, cohort = cohorts)
random$year_region_cohort_effect <- v

clg()

# Bubble plots:
dimnames(random$year_region_effect) <- list(year = years, region = regions)
windows()
bubble.plot(log(mu_year_region_cohort[,,2]), scale = 10)


clg()
for (j in 1:6){
   windows()
   plot(range(years), range(log(mu_year_region_cohort[,,j])), type = "n", ylab = "ln(density)", xlab = "Year")
   grid()
   cols <- rainbow(length(regions))
   for (i in 1:length(regions)){
      lines(years, log(mu_year_region_cohort[,i,j]), col = cols[i], lwd = 2)
   }
   legend("topleft", legend = regions, lwd = 2, col = cols)
   title(main = paste0("Cohort ", j))
}

# Bubble plots:
dimnames(random$year_cohort_effect) <- list(year = years, cohort = cohorts)
windows()
bubble.plot(random$year_cohort_effect, scale = 5)

# Bubble plots:
dimnames(random$region_cohort_effect) <- list(region = regions, cohort = cohorts)
windows()
bubble.plot(random$region_cohort_effect, scale = 5)

r <- summary(rep, "random")
random$year_region_cohort_effect <- r[rownames(r) == "year_region_cohort_effect", 1]
dim(random$year_region_cohort_effect) <- c(length(years), length(regions), length(cohorts))
dimnames(random$year_region_cohort_effect) <- list(year = years, region = regions, cohort = cohorts)

# Bubble plots:
for (i in 1:7){
   windows()
   bubble.plot(random$year_region_cohort_effect[,,i], scale = 5)
   title(main = paste0("Cohort ", cohorts[i]))
}
