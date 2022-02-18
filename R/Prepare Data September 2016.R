library(gulf)

setwd("U:/Lobster/Transect Analysis")
rm(list = ls())
source("U:/Lobster/Transect Analysis/Lobster functions.R")


# Read SCUBA lobster data:
x <- read.table("Lobster Transect Data 2016.csv", header = TRUE, sep = ",", colClasses = "character", stringsAsFactors = FALSE)

# Rename variable fields:
names(x) <- c("region", "marker", "day", "month", "year", "diver", "depth", "site", "transect", "section", 
              "sex", "length", "sandy", "muddy", "hard", "Al", "gravel", "cobble", "Ca", "mallet", "upm")

# Reformat fields:
x$sex <- tolower(x$sex)
x$year <- as.numeric(x$year)
x$length <- as.numeric(x$length)
x$mallet <- as.numeric(x$mallet)                      
x$transect <- paste(x$site, x$transect, sep = "-")

# Compile errors:
index <- which(x$section %in% c("", "116", "81", "21", "i"))
index <- c(index, which((x$day == "") | (x$diver == "")))
index <- sort(index)

# Remove errors:
x <- x[setdiff(1:dim(x)[1], index), ]
#x <- x[x$year < 2013, ]
x <- x[!((x$region == "SH") & (x$year == 2001)), ]
x <- x[!is.na(x$year), ]
x <- x[!(x$region %in% c("AB2", "GA2")), ]

# Classify cohorts:
x$year.class <- year.class(x)
table(x[c("year.class", "mallet")])

# Check 'mallet' indicator:
vars <- c("year", "site", "transect", "year.class")
temp1 <- aggregate(list(n = x$mallet), x[vars], sum, na.rm = TRUE)
temp2 <- aggregate(list(n = x$length), x[vars], function(x) sum(!is.na(x)))
i <- which((temp1$n > 0) & (temp1$n != temp2$n ))  

# Point corrections:
index <- which((x$year == temp1$year[i[1]]) & (x$transect == temp1$transect[i[1]]) & (x$year.class == temp1$year.class[i[1]]) )
x$mallet[index] <- 0
index <- which((x$year == temp1$year[i[2]]) & (x$transect == temp1$transect[i[2]]) & (x$year.class == temp1$year.class[i[2]]) )
x$mallet[index] <- 0

# Determine number of year classes within transects with Mallet indicators:
res <- unique(x[c("year", "site", "transect")])
vars <- c("year", "site", "transect")
for (i in 0:6){
   index <- which(x$year.class == i)
   temp <- aggregate(list(n = x$mallet[index]), x[index, vars], function(x) any(x == 1))
   index <- match(res[vars], temp[vars])
   res[paste0("n", i)] <- temp$n[index]
}
fvars <- paste0("n", 0:6)
temp <- res[fvars]
temp[is.na(temp)] <- FALSE
temp <- temp + 1 - 1
res[fvars] <- temp
res[apply(res[fvars], 1, sum) > 0, ]

# Calculate swept area:
key <- c("year", "region", "transect")
swept.area <- aggregate(list(swept.area = x$section), by = x[key], function(x) length(unique(x))*20)

# Prepare count table:
res <- unique(x[key])
res$year.class <- 0
for (i in 1:6){
   temp <- unique(x[key])
   temp$year.class <- i
   res <- rbind(res, temp)
}

# Calculate counts by year class:
temp <- aggregate(list(n = x$length), by = x[c(key, "year.class")], function(x) sum(!is.na(x)))
res <- merge(res, temp, all.x = TRUE)
res$n[is.na(res$n)] <- 0

# Blank Mallet factor:
temp <- aggregate(list(mallet = x$mallet), by = x[c(key, "year.class")], function(x) sum(x))
res <- merge(res, temp, all.x = TRUE)
res$mallet[is.na(res$mallet)] <- 0
res$n[res$mallet > 0] <- NA

res <- merge(res, swept.area, all.x = TRUE)

winbugs(res$year.label)
winbugs(res$region.label)
winbugs(res$class.label)
winbugs(res$swept.area)
winbugs(res$n)



 
