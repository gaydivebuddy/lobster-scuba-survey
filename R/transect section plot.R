# Read transect section table:
load("section.rdata")
s$year <- year(s$date)

# Load biological data:
load("biological.rdata")

# Calculate swept area:
swept.area <- aggregate(list(sections = s$interval), by = s[c("date", "region", "transect", "diver", "width_m", "depth_ft")], function(x) length(unique(x)))
swept.area$swept.area <-  5 * swept.area$sections * swept.area$width_m

tmp <- aggregate(swept.area$width_m, by = swept.area[c("date", "region", "transect", "diver")], length)


# Compile transect table by year:
vars <- c("year", "date", "region", "transect", "diver")
transect <- unique(s[vars])
transect <- transect[order(transect$year, transect$region, transect$date, transect$transect, transect$diver), ]

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
   pdf(file = paste0("Transect sampling ", years[i], ".pdf"), height = 8.5, width = 11)
   
   ix <- which(year(transect$date) == years[i]) 
   
   y <- t(res[, ix])
   u <- unique(transect$region[ix])
   u <- u[!is.na(u)]
   y <- y * repvec(match(transect$region[ix], u), ncol = 20)
   
   image(1:length(ix), 1:20, y, xlab = "", ylab = "", xaxt = "n", col = hcl.colors(10, palette = "dark3"))
   mtext("Section", 2, 2.5, cex = 1.5)
   mtext("Transect", 1, 2.0, cex = 1.5)
   
   tmp <- transect[ix, ]
   
   vars <- c("date", "region", "transect")
   k <- 1
   for (j in 2:nrow(tmp)){
      if (!any(is.na(tmp[j-1,vars])) & !any(is.na(tmp[j,vars]))){
         if (!all(tmp[j,vars] == tmp[j-1,vars])){
            vline(j - 0.5)
            k <- k + 1
           # if ((k %% 10) == 0) mtext(k, 1, at = j + 2)
         } 
      } 
   }
 
   mtext(years[i], 3, 1.75, cex = 1.5)

   for (j in 1:19){
      hline(j+0.5, lwd = 0.5, lty = "dashed", col = "grey70")
   }
   
   for (j in 1:length(u)){
      tmp <- which(transect$region[ix] == u[j])
      mtext(u[j], 3, 0.2 + 0.7*(j %% 2), at = mean(tmp), cex = 0.75)
      vline(c(min(tmp) - 0.5, max(tmp) + 0.5), lwd = 3)
   }
   
   # Add transect labels:
   tmp <- aggregate(list(at = 1:length(ix)), by = list(region = transect$region[ix], transect = transect$transect[ix]), mean)
   tmp <- tmp[order(tmp$region, tmp$transect), ]
   for (j in 1:nrow(tmp)){
      mtext(tmp$transect[j], 1, 0.1 + 0.7 * (j %% 2), at = tmp$at[j], cex = 0.75)
   }
   
   box()
   
   dev.off()
}
  
 



