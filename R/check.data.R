# Check data:


# Check SCUBA section table:
load("section.rdata")
s <- s[order(s$date, s$region, s$transect, s$interval_display), ]

s$updated_at <- substr(s$updated_at, 1, 19)
s$created_at <- substr(s$created_at, 1, 19)

vars <- c("date", "region", "transect", "interval_display")

# Check for sections which have been sampled by more than two divers:
tmp <- aggregate(list(n = s$diver), by = s[vars], length)
tmp <- tmp[order(tmp$date, tmp$region, tmp$transect, tmp$interval_display), ]
tmp[tmp$n > 2, ]

# Check for NA values in indexing variables:
s[which(apply(s[vars], 1, function(x) any(is.na(x)))), c(vars, "diver", "dive", "depth_ft", "updated_at")]

# Define habitat variables:
hvars <- c("percent_algae", "percent_sand", "percent_mud", "percent_hard", "percent_gravel", "percent_cobble", "percent_pebble")

# Show table
aggregate(s[hvars], by = list(year = year(s$date)), sum)

# Show number of sections with no habitat designations by year:
table(year(s$date), apply(s[hvars[]], 1, sum) > 0)

s[which((year(s$date) == 2019) & (apply(s[hvars[]], 1, sum) == 0)), c(vars, "diver", "dive")]

# Percentage of sections with no habitat designation:
aggregate(apply(s[hvars], 1, sum), list(year(s$date)), function(x) round(100*sum(x==0)/length(x)))

# - almost no habitat data for 2007-2018, but a lot of missing sections.
table(year(s$date), sum(apply(s[hvars], 1, sum) == 0) / sum(apply(s[hvars], 1, sum) > 0))


