# Update data from DMApps database:


# Dive log (does not work):
file.copy(from = "http://dmapps/en/scuba/reports/dive.log.xlsx", to = "dive log.xlsx")

# Read in transect section table:
sections <- read.csv("http://dmapps/en/scuba/reports/section/?year=2021")

# Read dives table:
dives <- read.csv("http://dmapps/en/scuba/reports/dive/?year=2021")

# Read outings table:
outings <- read.csv("http://dmapps/en/scuba/reports/outing/?year=2021")

# Read biological data:
observations <- read.csv("http://dmapps/en/scuba/reports/observations/?year=2018")

# Scuba transect table:
read.csv("http://dmapps/en/scuba/reports/scuba_transect/?year=2021")

read.scuba(year, ...)