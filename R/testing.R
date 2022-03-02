# Read transect section table:
load("transect.rdata")

# Load biological data:
load("biological.rdata")

# Calculate swept area:
swept.area <- aggregate(list(sections = s$interval), by = s[c("date", "region", "transect", "diver", "width_m", "depth_ft")], function(x) length(unique(x)))
swept.area$swept.area <-  5 * swept.area$sections * swept.area$width_m

tmp <- aggregate(swept.area$width_m, by = swept.area[c("date", "region", "transect", "diver")], length)


# Compile transect table by year:
vars <- c("date", "region", "transect", "diver")
transect <- unique(s[vars])
transect <- transect[order(transect$date, transect$region, transect$transect, transect$diver), ]

res <- matrix(NA, nrow = 20, ncol = nrow(transect))
for (i in 1:nrow(transect)){
   ix <- rep(TRUE, nrow(s))
   for (j in 1:length(vars)) ix <- ix & (transect[i, vars[j]] == s[, vars[j]])
   ix <- which(ix)
   if (length(ix) > 0){
      res[s$interval[ix], i] <- 1
   }
}

image(t(res))

clg()
years <- sort(unique(year(s$date)))
for (i in 1:length(years)){
   ix <- which(year(transect$date) == years[i]) 
   image(1:length(ix), 1:20, t(res[, ix]), xlab = "", ylab = "", xaxt = "n")
   mtext("Section", 2, 2.5, cex = 1.5)
   mtext("Transect", 1, 1.5, cex = 1.5)
   
   tmp <- transect[ix, ]
   
   vars <- c("date", "region", "transect")
   k <- 1
   for (j in 2:nrow(tmp)){
      if (!any(is.na(tmp[j-1,vars])) & !any(is.na(tmp[j,vars]))){
         if (!all(tmp[j,vars] == tmp[j-1,vars])){
            vline(j - 0.5)
            k <- k + 1
            if ((k %% 10) == 0) mtext(k, 1, at = j + 2)
         } 
      } 
   }
 
   mtext(years[i], 3, 1.5, cex = 1.5)
}
  
 



