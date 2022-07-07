library(gulf.data)
library(gulf.graphics)

precision <- 1

# Read Scuba data:
# x <- read.scuba(2000:2021)
load("data/data 2000-2021.rdata")

# Load Comeau exclusions table:
#exclusions <- read.csv(paste0(gsub("/R", "/data/", getwd()), "comeau exclusions.csv"))
exclusions <- read.csv(paste0(getwd(), "/data/", "comeau exclusions.csv"))

names(exclusions) <- c("year", "transect.name", paste0("cohort", 0:6))
for (i in 1:nrow(exclusions)){
   ix <- grep(exclusions$transect.name[i], x$transects$transect.name.old)
   if (length(ix) > 0){
      exclusions$transect.id[i]    <- list(x$transects$transect.id[ix])
      exclusions$transect.names[i] <- list(x$transects$transect.name.old[ix])
   }else{
      exclusions$transect.id[i] <- list(NA)
      exclusions$transect.names[i] <- list(NA)
   }
}

# Define which transects were removed completely and which were the product of UPM seeding:
seeded <- exclusions[apply(exclusions[, grep("cohort", names(exclusions))], 1, sum) < 6, ]
exclusions <- exclusions[apply(exclusions[, grep("cohort", names(exclusions))], 1, sum) >= 6, ]

# Import exclusion flags to dive table:
ix <- match(x$dives$transect.name, x$transects$transect.name)
x$dives$transect.id <- x$transects$transect.id[ix]
x$dives$exclude <- FALSE
for (i in 1:nrow(x$dives)){
   if (x$dives$date[i] != ""){
      tmp <- exclusions[year(x$dives$date[i]) == exclusions$year, ]
      if (nrow(tmp) > 0){
         exclude <- unlist(tmp$transect.id)
         exclude <- sort(exclude[!is.na(exclude)])
         if (x$dives$transect.id[i] %in% exclude) x$dives$exclude[i] <- TRUE 
      } 
   }
}
ix <- match(x$outings$outing.id, x$dives$outing.id)
x$outings$exclude <- x$dives$exclude[ix]
x$outings$exclude[is.na(x$outings$exclude)] <- FALSE
ix <- match(x$sections$dive.id, x$dives$dive.id)
x$sections$exclude <- x$dives$exclude[ix]
x$sections$exclude[is.na(x$sections$exclude)] <- FALSE
ix <- match(x$observations$dive.id, x$dives$dive.id)
x$observations$exclude <- x$dives$exclude[ix]
x$observations$exclude[is.na(x$observations$exclude)] <- FALSE

# Import additional 'was.seeded' flags to dive table:
#x$dives$seeded <- FALSE
for (i in 1:nrow(x$dives)){
   if (x$dives$date[i] != ""){
      tmp <- seeded[year(x$dives$date[i]) == seeded$year, ]
      if (nrow(tmp) > 0){
         was.seeded <- unlist(tmp$transect.id)
         was.seeded <- sort(was.seeded[!is.na(was.seeded)])
         if (x$dives$transect.id[i] %in% was.seeded) x$dives$was.seeded[i] <- TRUE 
      } 
   }
}
ix <- match(x$outings$outing.id, x$dives$outing.id)
x$outings$was.seeded <- x$dives$was.seeded[ix]
x$outings$was.seeded[is.na(x$outings$was.seeded)] <- FALSE
ix <- match(x$sections$dive.id, x$dives$dive.id)
x$sections$was.seeded <- x$dives$was.seeded[ix]
x$sections$was.seeded[is.na(x$sections$was.seeded)] <- FALSE
ix <- match(x$observations$dive.id, x$dives$dive.id)
x$observations$was.seeded <- x$dives$was.seeded[ix]
x$observations$was.seeded[is.na(x$observations$was.seeded)] <- FALSE

# Define regions to remove from analysis:
remove.regions <- c("Grande-Anse", "Pigeon Hill", "Anse-Bleue", "Robichaud", "", "Alberton", "Nine Mile Creek")

# Keep only completed sections and those not part of diet study sampling:
x$sections <- x$sections[which(x$sections$completed & x$sections$transect.name != ""), ]
x$sections$swept.area <- 5 * x$sections$width.m
x$sections$side.display[is.na(x$sections$side.display)] <- ""

# Remove irrelevant regions:
x$outings      <- x$outings[which(!(x$outings$region %in% remove.regions)), ]
x$dives        <- x$dives[which(!(x$dives$region %in% remove.regions)), ]
x$transects    <- x$transects[which(!(x$transects$region %in% remove.regions)), ]
x$sections     <- x$sections[which(!(x$sections$region %in% remove.regions)), ]
x$observations <- x$observations[which(!(x$observations$region %in% remove.regions)), ]

# Remove training dive data:
x$dives        <- x$dives[which(!x$dives$is.training), ]
x$sections     <- x$sections[which(!x$sections$is.training), ]
x$observations <- x$observations[which(!x$observations$is.training), ]

# Remove UPM data:
x$dives        <- x$dives[which(!x$dives$is.upm), ]
x$sections     <- x$sections[which(!x$sections$is.upm), ]
x$observations <- x$observations[which(!x$observations$is.upm), ]

# Remove specific divers:
remove <- c("Godin, M", "Asselin, Natalie", "Landry, Germain", "Leblanc, Stepan", "Landry, Eric", "Cousineau, Maryse")
x$dives        <- x$dives[which(!(x$dives$diver %in% remove)), ]
x$sections     <- x$sections[which(!(x$sections$diver %in% remove)), ]
x$observations <- x$observations[which(!(x$observations$diver %in% remove)), ]

# Clean up dates:
ix <- which(x$dives$date == "")
x$dives$date[ix] <- x$outings$date[match(x$dives$outing.id[ix], x$outings$outing.id)]

# Determine observation precisions:
x$observations$precision <- 1
x$observations$precision[(year(x$observations$date) %in% c(2000:2018)) & (x$observations$certainty.rating == 0) & (x$observations$carapace.length.mm < 30)] <- 2
x$observations$precision[(year(x$observations$date) %in% c(2000:2018)) & (x$observations$certainty.rating == 0) & (x$observations$carapace.length.mm >= 30)] <- 5
x$observations$precision[(year(x$observations$date) %in% c(2019,2021)) & (x$observations$certainty.rating == 0)] <- 5

# Remove excluded transects and data:
x$outings      <- x$outings[!x$outings$exclude, ]
x$dives        <- x$dives[!x$dives$exclude, ]
x$sections     <- x$sections[!x$sections$exclude, ]
x$observations <- x$observations[!x$observations$exclude, ]

# Remove tows with short existence in Scuba survey:
remove <- c(paste0("T-", c(11, 12, 21, 22, 31, 32, 41, 42, 61, 62), " (Neguac)"),
            paste0("T-", c(1:4, 6, 9, 12:42), " (Fox Harbor)"),
            paste0("T-", c(1, 10:13, 15, 17, 2, 22:26, 28:29, 3, 30:34, 4, 8:9), " (Caraquet)"),
            paste0("T-", c(1, 2, 22, 3, 32, 4), " (Pointe-Verte)"), 
            paste0("T-", c(7:9), " (Richibucto)"),
            paste0("T-", c(7, 10:12), " (Shediac)"))
x$outings      <- x$outings[!(x$outings$transect.name %in% remove), ]
x$dives        <- x$dives[!(x$dives$transect.name %in% remove), ]
x$sections     <- x$sections[!(x$sections$transect.name %in% remove), ]
x$observations <- x$observations[!(x$observations$transect.name %in% remove), ]

# Remove 'was.seeded' tows:
x$dives <- x$dives[which(!x$dives$was.seeded), ]
x$sections <- x$sections[which(!is.na(match(x$sections$dive.id, x$dives$dive.id))), ]

# Remove data that have no corresponding transect in the sections table:
x$outings      <- x$outings[x$outings$transect.name %in% unique(x$sections$transect.name), ]
x$transects    <- x$transects[x$transects$transect.name %in% unique(x$sections$transect.name), ]
x$dives        <- x$dives[x$dives$transect.name %in% unique(x$sections$transect.name), ]
x$observations <- x$observations[x$observations$transect.name %in% unique(x$sections$transect.name), ]
x$observations <- x$observations[x$observations$species == "American lobster", ]
x$observations <- x$observations[x$observations$section.id %in% unique(x$sections$section.id), ]
x$observations <- x$observations[!is.na(x$observations$carapace.length.mm), ]

# Attach proportion of sampled sections per transect:
res <- aggregate(list(p = x$sections$interval), by = x$sections[c("region", "date", "transect.name")], function(x) length(unique(x)) / (diff(range(x))+1))
res$p[which(year(res$date) >= 2019)] <- 1
ix <- match(x$sections[c("region", "date", "transect.name")], res[c("region", "date", "transect.name")])
x$sections$proportion.sampled <- res$p[ix]

# Remove years prior to 2003:
x$outings      <- x$outings[year(x$outings$date) >= 2003, ]
x$dives        <- x$dives[year(x$dives$date) >= 2003, ]
x$sections     <- x$sections[year(x$sections$date) >= 2003, ]
x$observations <- x$observations[year(x$observations$date) >= 2003, ]

# Export data to Excel:
# excel(x$transects)
# excel(x$outings)
# excel(x$dives)
# excel(x$sections)
# excel(x$observations)

# Calculate size-frequencies by section, sex and measurement precision:
f <- freq(x$observations$carapace.length.mm, by = x$observations[c("section.id", "sex", "precision")])
fvars <- names(f)[gsub("[0-9]", "", names(f)) == ""]

# Spread-out uncertain size observations:
smooth <- function(v, precision, scale = 0.4){
   ix <- which(v > 0)
   p <- rep(0, length(fvars))
   for (i in 1:length(ix)){
      p <- p + v[[ix[i]]] * (pnorm(as.numeric(fvars)+0.5, as.numeric(fvars[ix[i]]), scale * precision)-
                             pnorm(as.numeric(fvars)-0.5, as.numeric(fvars[ix[i]]), scale * precision))
      
   }
   p <- sum(v) * p / sum(p)
}
for (i in 1:nrow(f)){
   if (f$precision[i] > 1){
      v <- smooth(f[i,fvars], f$precision[i])
      f[i,fvars] <- v
   }
}

# Sum out precision variable:
f <- aggregate(f[fvars], by = f[c("section.id", "sex")], sum)

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
names(p) <- fvars

# Parse unassigned sex category to males and females:
fum <- fuf <- f[f$sex == "u", ]
fum[fvars] <- fum[fvars] * repvec(1-p, nrow = nrow(fum))
fuf[fvars] <- fuf[fvars] * repvec(p, nrow = nrow(fuf))
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

# Linearize size-frequency data:
vars <- c("region", "date", "transect.name", "diver", "section.id", "side.display", "interval", "depth.ft", "proportion.sampled", "visibility", "swept.area")
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

# Reduce data size:
rm <- rm[rm$carapace.length <= 120, ]
rf <- rf[rf$carapace.length <= 120, ]
rm$carapace.length <- precision * round(rm$carapace.length / precision)
rf$carapace.length <- precision * round(rf$carapace.length / precision)
vars <- c("region", "date", "transect.name", "diver", "section.id", "carapace.length", "proportion.sampled")
rm <- aggregate(rm["n"], by = rm[vars], sum)
rf <- aggregate(rf["n"], by = rf[vars], sum)

# Reduce data set to numbers by diver at each transect:
vars <- c("region", "date", "transect.name", "diver", "proportion.sampled")
area <- aggregate(x$sections["swept.area"],  by = x$sections[vars], sum)
mm <- aggregate(rm[c("n")], by = rm[c(vars, "carapace.length")], sum)
ix <- match(mm[vars], area[vars])
mm$swept.area <- area$swept.area[ix]
ff <- aggregate(rf[c("n")], by = rf[c(vars, "carapace.length")], sum)
ix <- match(ff[vars], area[vars])
ff$swept.area <- area$swept.area[ix]

# Attach 'was.seeded' flag:
vars <- c("region", "date", "transect.name", "diver", "proportion.sampled")
u <- unique(x$sections[c(vars, "was.seeded")])
ix <- match(mm[vars], u[vars])
mm$seeded <- u$was.seeded[ix]
ix <- match(ff[vars], u[vars])
ff$seeded <- u$was.seeded[ix]

# Standardize size-frequencies by swept area:
#sm[fvars] <- sm[fvars] / repvec(sm$swept.area, ncol = length(fvars))
#sf[fvars] <- sf[fvars] / repvec(sm$swept.area, ncol = length(fvars))

