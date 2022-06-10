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
theta <- c(xp = 60.93761, slope = 10.14174, log_a = 0.075303, log_b = -0.671195)
p <- mu(as.numeric(fvars), theta) # Proportion of females by size, as per 'sex ratio analysis.R'.

# Parse unassigned sex category to males and females:
fu <- fum <- fuf <- f[f$sex == "u", ]
fum[fvars] <- fu[fvars] * repvec(1-p, nrow = nrow(fu))
fuf[fvars] <- fu[fvars] * repvec(p, nrow = nrow(fu))
fm <- f[f$sex == "m", ]
ff <- f[f$sex == "f", ]
ix <- match(fum$section.id, fm$section.id)
fm[ix[!is.na(ix)], fvars] <- fm[ix[!is.na(ix)], fvars] + fum[!is.na(ix),fvars]
fm <- rbind(fm, fum[is.na(ix),]) # Add sections that had no existing entries.
ix <- match(fuf$section.id, ff$section.id)
ff[ix[!is.na(ix)], fvars] <- ff[ix[!is.na(ix)], fvars] + fuf[!is.na(ix),fvars]
ff <- rbind(ff, fuf[is.na(ix),]) # Add sections that had no existing entries.

# Attach size-frequencies to sections table:
sf <- sm <- x$sections
import(sm, fill = 0) <- fm
import(sf, fill = 0) <- ff

vars <- c("region", "date", "transect.name", "diver", "section.id", "side.display", "interval", "depth.ft", "visibility", "swept.area")
for (i in 1:length(vars)){
   m.tmp <- data.frame(rep(sm[,vars[i]], each = length(fvars)), stringsAsFactors = FALSE)
   f.tmp <- data.frame(rep(sf[,vars[i]], each = length(fvars)), stringsAsFactors = FALSE)
   names(m.tmp) <- vars[i]
   names(f.tmp) <- vars[i]
   if (i == 1){
      rm <- m.tmp
      rf <- f.tmp
   }else{
      rm <- cbind(rm, m.tmp)
      rf <- cbind(rf, f.tmp)
   } 
}
rm$carapace.length <- as.vector(t(repvec(as.numeric(fvars), nrow = nrow(sm))))
rm$n <- as.vector(t(as.matrix(sm[fvars])))
rf$carapace.length <- as.vector(t(repvec(as.numeric(fvars), nrow = nrow(sf))))
rf$n <- as.vector(t(as.matrix(sf[fvars])))

rm <- rm[rm$n > 0, ]
rf <- rf[rf$n > 0, ]

# Standardize size-frequencies by swept area:
#sm[fvars] <- sm[fvars] / repvec(sm$swept.area, ncol = length(fvars))
#sf[fvars] <- sf[fvars] / repvec(sm$swept.area, ncol = length(fvars))

# Linearize data table:

# Parse unassigned sex category to males or females:
x$observations$n <- 1
tmp <- x$observations[x$observations$sex == "u", ]
x$observations <- x$observations[x$observations$sex != "u", ]

# Parse
m <- f <- tmp
m$sex <- "m"
f$sex <- "f"
m$n <- 1-p[as.character(m$carapace.length.mm)]
f$n <- p[as.character(f$carapace.length.mm)]

x$observations <- rbind(x$observations, m, f)

vars <- c("region", "date", "transect.name", "section.id", "diver", "side.display", "interval", "certainty.rating", "sex", "carapace.length.mm")

vars %in% names(x$observations)

tmp <- aggregate(x$observations["n"], by = x$observations[vars], sum)

z <- tmp[tmp$certainty.rating == 0, ]

f <- freq(z$carapace.length.mm, by = z[c("section.id", "diver", "sex")])

fvars <- names(f)[gsub("[0-9]", "", names(f)) == ""]

smooth <- function(v){
   ix <- which(v > 0)
   p <- rep(0, length(fvars))
   for (i in 1:length(ix)){
      p <- p + v[[ix[i]]] * (pnorm(as.numeric(fvars)+0.5, as.numeric(fvars[ix[i]]), 2)-
                                pnorm(as.numeric(fvars)-0.5, as.numeric(fvars[ix[i]]), 2))
      
   }
   p <- sum(v) * p / sum(p)
}

for (i in 1:nrow(f)){
   print(i)
   v <- smooth(f[i,fvars])
   f[i,fvars] <- v
}
gbarplot(apply(f[fvars], 2, sum))
vline(seq(0, 100, by = 5), lty = "dashed", col = "red")


