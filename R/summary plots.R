#load("data 2010-2021.rdata")

# Diver time-line:
tab <- table(year(x$sections$date), x$sections$diver)
tab <- tab[, rev(1:ncol(tab))]
years <- 2001:2021
mm <- matrix(0, nrow = length(setdiff(as.character(years), rownames(tab))), ncol = ncol(tab))
rownames(mm) <- setdiff(as.character(years), rownames(tab))
colnames(mm) <- colnames(tab)
tab <- rbind(tab, mm)
tab <- tab[order(as.numeric(rownames(tab))), ]
m <- cbind(matrix(0, nrow = 7, ncol = 1), matrix(1, nrow = 7, ncol = 7))
m <- cbind(0, rbind(0, m, 0), 0)
layout(m)
par(mar = c(0,0,0,0))
image(1:nrow(tab), 1:ncol(tab), tab > 0, col = c("white", "grey70"), xaxt = "n", yaxt = "n")
axis(1, at = 1:nrow(tab), labels = rownames(tab), cex.axis = 1.25)
axis(2, at = 1:ncol(tab), labels = colnames(tab), las = 2, cex.axis = 1.25)
vline(0.5:21, col = "grey50")
hline(0.5:21, col = "grey50")
box(col = "grey50")
mtext("Diver", 2, 10, cex = 1.5)
mtext("Year", 1, 3.5, cex = 1.5)

# Diver pairings within transects:
divers <- sort(unique(x$sections$diver))
tab <- unique(x$sections[c("diver", "transect.name", "date")])
mm <- matrix(0, nrow = length(setdiff(as.character(years), rownames(tab))), ncol = ncol(tab))
rownames(mm) <- setdiff(as.character(years), rownames(tab))
colnames(mm) <- colnames(tab)
tab <- rbind(tab, mm)
tab <- tab[order(as.numeric(rownames(tab))), ]
r <- matrix(0, nrow = length(divers), ncol = length(divers))
dimnames(r) <- list(diver = divers, diver = divers)
u <- unique(tab[c("transect.name", "date")])
for (i in 1:nrow(u)){
   ix <- which(tab$date == u$date[i] & tab$transect.name == u$transect.name[i])
   if (length(ix) >= 2){
      ix <- expand.grid(ix, ix)
      ix <- ix[ix[,1] != ix[,2], ]
      for (j in 1:nrow(ix)){
         r[tab$diver[ix[j,1]], tab$diver[ix[j,2]]] <- r[tab$diver[ix[j,1]], tab$diver[ix[j,2]]] + 1
      }
   }
}
image(1:length(divers), 1:length(divers), r > 0, col = c("white", "grey70"), xaxt = "n", yaxt = "n")
axis(1, at = 1:length(divers), labels = divers, las = 3, cex.axis = 1.25)
axis(2, at = 1:length(divers), labels = divers, las = 2, cex.axis = 1.25)
vline(0.5:20, col = "grey50")
hline(0.5:20, col = "grey50")
box(col = "grey50")
mtext("Diver", 2, 10, cex = 1.5)

# Region time-line:
tab <- table(year(x$sections$date), x$sections$region)
years <- 2001:2021
tab <- tab[, rev(1:ncol(tab))]
mm <- matrix(0, nrow = length(setdiff(as.character(years), rownames(tab))), ncol = ncol(tab))
rownames(mm) <- setdiff(as.character(years), rownames(tab))
colnames(mm) <- colnames(tab)
tab <- rbind(tab, mm)
tab <- tab[order(as.numeric(rownames(tab))), ]
m <- cbind(matrix(0, nrow = 7, ncol = 1), matrix(1, nrow = 7, ncol = 7))
m <- cbind(0, rbind(0, m, 0), 0)
layout(m)
par(mar = c(0,0,0,0))
image(1:nrow(tab), 1:ncol(tab), tab > 0, col = c("white", "grey70"), xaxt = "n", yaxt = "n")
axis(1, at = 1:nrow(tab), labels = rownames(tab), cex.axis = 1.25)
axis(2, at = 1:ncol(tab), labels = colnames(tab), las = 2, cex.axis = 1.25)

vline(0.5:20, col = "grey50")
hline(0.5:20, col = "grey50")
box(col = "grey50")
mtext("Region", 2, 10, cex = 1.5)
mtext("Year", 1, 3.5, cex = 1.5)

# Transect time-lines:
years <- 2001:2021
for (i in 1:length(regions)){
   xx <- x$sections[x$sections$region == regions[i], ]
   
   tab <- table(year(xx$date), xx$transect.name)
   tab <- tab[, rev(1:ncol(tab))]
   tab <- as.matrix(tab)
   mm <- matrix(0, nrow = length(setdiff(as.character(years), rownames(tab))), ncol = ncol(tab))
   rownames(mm) <- setdiff(as.character(years), rownames(tab))
   colnames(mm) <- colnames(tab)
   tab <- rbind(tab, mm)
   tab <- tab[order(as.numeric(rownames(tab))), ]
   
   
   m <- cbind(matrix(0, nrow = 7, ncol = 1), matrix(1, nrow = 7, ncol = 7))
   m <- cbind(0, rbind(0, m, 0), 0)
   layout(m)
   par(mar = c(0,0,0,0))
   image(1:nrow(tab), 1:ncol(tab), tab > 0, col = c("white", "grey70"), xaxt = "n", yaxt = "n")
   axis(1, at = 1:nrow(tab), labels = rownames(tab), cex.axis = 1.25)
   axis(2, at = 1:ncol(tab), labels = colnames(tab), las = 2, cex.axis = 1.25)
   
   vline(0.5:nrow(tab), col = "grey50")
   hline(0.5:ncol(tab), col = "grey50")
   box(col = "grey50")
   mtext("Region", 2, 10, cex = 1.5)
   mtext("Year", 1, 3.5, cex = 1.5)
   mtext(regions[i], 3, 1.5, cex = 1.5)
}



for (i in 1:length(regions)){
   xx <- x$sections[x$sections$region == regions[i], ]
   tab <- table(xx$date, xx$trasect)
   
   
}


