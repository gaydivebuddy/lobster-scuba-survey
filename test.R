library(TMB)
library(gulf.data)

Sys.setenv(PATH = paste0(Sys.getenv("PATH"), ";", "C:\\rtools40\\mingw64\\bin"))

update.scuba <- function(){
   # Read in transect section table:
   s <- read.csv("http://dmapps/en/scuba/reports/section/?year=2021")

   t <- read.csv("http://dmapps/en/scuba/reports/section/?year")

   # Read biological data:
   b <- read.csv("http://dmapps/en/scuba/reports/observations/?year")
}

# Read transect section table:
load("C:/Users/SuretteTj/Desktop/github/lobster-scuba-project/transect.rdata")

# Load biological data:
load("C:/Users/SuretteTj/Desktop/github/lobster-scuba-project/biological.rdata")

# Calculate swept area:
swept.area <- aggregate(list(sections = s$interval), by = s[c("date", "region", "transect", "diver", "width_m", "depth_ft")], function(x) length(unique(x)))
swept.area$swept.area <-  5 * swept.area$sections * swept.area$width_m

# Compile transect table by year:
transect <- unique(s[c("date", "region", "transect", "width_m", "depth_ft")])
transect <- transect[order(transect$date), ]



compile("test.cpp")

dyn.load(dynlib("NB_lobster_transect_additive"))
dyn.load(dynlib("NB_lobster_transect_interaction"))
dyn.load(dynlib("NB_lobster_transect_full"))



