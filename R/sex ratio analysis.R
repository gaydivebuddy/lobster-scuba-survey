# Read Scuba data:
x <- read.scuba(2021)

# Keep only completed sections and those not part of diet study sampling:
x$sections <- x$sections[which(x$sections$completed & x$sections$transect.name != ""), ]
x$sections$swept.area <- 5 * x$sections$width.m

# Remove data that have no corresponding transect in the sections table:
x$outings      <- x$outings[x$outings$transect.name %in% unique(x$sections$transect.name), ]
x$transects    <- x$transects[x$transects$name %in% unique(x$sections$transect.name), ]
x$dives        <- x$dives[x$dives$transect.name %in% unique(x$sections$transect.name), ]
x$observations <- x$observations[x$observations$transect.name %in% unique(x$sections$transect.name), ]
x$observations <- x$observations[x$observations$species == "American lobster", ]
x$observations <- x$observations[x$observations$section.id %in% unique(x$sections$section.id), ]
x$observations <- x$observations[!is.na(x$observations$carapace.length.mm), ]

# Calculate size-frequencies by section and sex:
f <- freq(x$observations$carapace.length.mm, by = x$observations[c("section.id", "sex")])
fvars <- names(f)[gsub("[0-9]", "", names(f)) == ""]

# Attach size-frequencies to sections table:
su <- sf <- sm <- x$sections
import(sm, fill = 0) <- f[f$sex == "m", ]
import(sf, fill = 0) <- f[f$sex == "f", ]
import(su, fill = 0) <- f[f$sex == "u", ]

# Standardize size-frequencies by swept area:
sm[fvars] <- sm[fvars] / repvec(sm$swept.area, ncol = length(fvars))
sf[fvars] <- sf[fvars] / repvec(sm$swept.area, ncol = length(fvars))
su[fvars] <- su[fvars] / repvec(sm$swept.area, ncol = length(fvars))

# Average size-frequencies by transect:
rm <- aggregate(sm[fvars], by = sm["transect.name"], mean)
rf <- aggregate(sf[fvars], by = sf["transect.name"], mean)
ru <- aggregate(su[fvars], by = su["transect.name"], mean)

# Global size-frequencies by sex:
gm <- apply(rm[fvars], 2, mean)
gf <- apply(rf[fvars], 2, mean)
gu <- apply(ru[fvars], 2, mean)

# Modified logistic model:
mu <- function(x, theta){
   # Extract model parameters:
   xp <- theta[["xp"]]
   slope <- theta[["slope"]]
   a <- 1 / (1 + exp(-theta[["log_a"]]))
   b <- 1 / (1 + exp(-theta[["log_a"]] - exp(theta[["log_b"]])))
   
   # Bernoulli probabilities:
   p <- (b-a) / (1 + exp(-slope*(x-xp))) + a
}

# Binomial log-likelihood:
loglike <- function(theta, x, y, n){
   # Evaluate modified logistic model:
   p <- mu(x, theta)
   
   # Bernoulli likelihood:
   ll <- n * (y * log(p) + (1-y) * log(1-p))
   
   return(-sum(ll))
}

# Initialize parameter values:
theta <- c(xp = 50, slope = 0.5, log_a = 0, log_b = 0)

# Prepare input data:
n <- apply(f[f$sex == "m", -1:-2], 2, sum) + apply(f[f$sex == "f", -1:-2], 2, sum)
y <- gf / (gm + gf)
n <- n[!is.na(y)]
y <- y[!is.na(y)]

# Fit model:
loglike(theta, x = as.numeric(names(n)), y = y, n = n)
theta <- optim(theta, loglike, x = as.numeric(names(n)), y = y, n = n, control = list(trace = 3))$par

# Plot results:
gbarplot(gf / (gm + gf), ylim = c(0, 1.0), xlim = c(20, 100))
hline(0.5, col = "red", lty = "dashed", lwd = 2)
lines(t, mu(t, theta), lwd = 2, col = "red")

# Create sex assignment vector:
t <- 0:500
p <- mu(t, theta)
names(p) <- t


