#load("data 2010-2021.rdata")

# Diver time-line:
language <- language("en")
tiff(file = paste0("figures/diver timeline - ", language, ".tiff"), compression = "lzw", units = "in", res = 300, height = 5.5, width = 7.5)
m <- rbind(0, cbind(0, 0, 0, 0, matrix(1, nrow = 10, ncol = 10), 0), 0, 0)
layout(m)
par(mar = c(0,0,0,0))
tab <- table(year(x$sections$date), x$sections$diver)
tab <- tab[, rev(1:ncol(tab))]
years <- 2001:2021
mm <- matrix(0, nrow = length(setdiff(as.character(years), rownames(tab))), ncol = ncol(tab))
rownames(mm) <- setdiff(as.character(years), rownames(tab))
colnames(mm) <- colnames(tab)
tab <- rbind(tab, mm)
tab <- tab[order(as.numeric(rownames(tab))), ]
image(1:nrow(tab), 1:ncol(tab), tab > 0, col = c("white", "grey70"), xaxt = "n", yaxt = "n")
axis(1, at = 1:nrow(tab), labels = rownames(tab), cex.axis = 1.25)
axis(2, at = 1:ncol(tab), labels = colnames(tab), las = 2, cex.axis = 1.25)
vline(0.5:21, col = "grey50")
hline(0.5:21, col = "grey50")
box(col = "grey50")
mtext("Diver", 2, 10, cex = 1.5)
mtext("Year", 1, 3.5, cex = 1.5)
dev.off()

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
language <- language("en")
tiff(file = paste0("figures/region timeline - ", language, ".tiff"), compression = "lzw", units = "in", res = 300, height = 4.5, width = 7.5)
m <- rbind(0, cbind(0, 0, 0,  matrix(1, nrow = 10, ncol = 10), 0), 0, 0)
layout(m)
par(mar = c(0,0,0,0))
tab <- table(year(x$sections$date), x$sections$region)
years <- 2001:2021
tab <- tab[, rev(1:ncol(tab))]
mm <- matrix(0, nrow = length(setdiff(as.character(years), rownames(tab))), ncol = ncol(tab))
rownames(mm) <- setdiff(as.character(years), rownames(tab))
colnames(mm) <- colnames(tab)
tab <- rbind(tab, mm)
tab <- tab[order(as.numeric(rownames(tab))), ]
image(1:nrow(tab), 1:ncol(tab), tab > 0, col = c("white", "grey70"), xaxt = "n", yaxt = "n")
axis(1, at = 1:nrow(tab), labels = rownames(tab), cex.axis = 1.25)
axis(2, at = 1:ncol(tab), labels = colnames(tab), las = 2, cex.axis = 1.25)
vline(0.5:20, col = "grey50")
hline(0.5:20, col = "grey50")
box(col = "grey50")
mtext("Region", 2, 9, cex = 1.5)
mtext("Year", 1, 3.5, cex = 1.5)
dev.off()

# Transect time-lines:
years <- 2003:2021
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
   
   
   language <- language("en")
   tiff(file = paste0("figures/transect timeline - ", regions[i], " - ", language, ".tiff"), compression = "lzw", units = "in", res = 300, height = 4.5, width = 7.5)
   m <- rbind(0, cbind(0, 0, 0,  matrix(1, nrow = 10, ncol = 10), 0), 0, 0)
   layout(m)
   par(mar = c(0,0,0,0))
   
   image(1:nrow(tab), 1:ncol(tab), tab > 0, col = c("white", "grey70"), xaxt = "n", yaxt = "n")
   axis(1, at = 1:nrow(tab), labels = rownames(tab), cex.axis = 1.25)
   axis(2, at = 1:ncol(tab), labels = colnames(tab), las = 2, cex.axis = 1.25)
   
   vline(0.5:nrow(tab), col = "grey50")
   hline(0.5:ncol(tab), col = "grey50")
   box(col = "grey50")
   # mtext("Region", 2, 10, cex = 1.5)
   mtext("Year", 1, 3.5, cex = 1.5)
   #mtext(regions[i], 3, 1.5, cex = 1.5)
   dev.off()
}

#res <- aggregate(list(n = apply(sm[fvars], 1, sum)), by = sm[c("region", "date", "transect.name", "interval")], sum)
#res$year <- year(res$date)
#res <- res[res$year >= 2019, ]

#aa <- aggregate(list(mu = res$n), by = res[c("region", "year", "transect.name")], mean)
#aa$sigma <- aggregate(list(sigma = res$n), by = res[c("region", "year", "transect.name")], sd)$sigma
#aa$n <- aggregate(list(n = res$n), by = res[c("region", "year", "transect.name")], length)$n

#plot(aa$sigma / aa$mu)

# Plot of maximum section sampled by year:
res <- aggregate(list(max = sm$interval), by = sm[c("region", "date", "transect.name")], max)
res <- aggregate(list(max = sm$interval), by = sm[c("region", "date", "transect.name")], function(x) diff(range(x))+1)
res$year <- year(res$date)
a <- aggregate(list(n = res$max), by = res["year"], mean)
a$n <- round(a$n, 1)
a$lower <- aggregate(list(n = res$max), by = res["year"], quantile, p = 0.25)$n
a$upper <- aggregate(list(n = res$max), by = res["year"], quantile, p = 0.75)$n
gbarplot(a$n, a$year, ylim = c(0, 20), yaxs = "i")
error.bar(a$year, lower = a$lower, upper = a$upper)
mtext("Year", 1, 2.5)
mtext("Maximum section sampled", 2, 2.5)

# Plot of average number of sections by year:
res <- aggregate(list(max = sm$interval), by = sm[c("region", "date", "transect.name")], function(x) length(unique(x)))
res$year <- year(res$date)
a <- aggregate(list(n = res$max), by = res["year"], mean)
a$n <- round(a$n, 1)
a$lower <- aggregate(list(n = res$max), by = res["year"], quantile, p = 0.25)$n
a$upper <- aggregate(list(n = res$max), by = res["year"], quantile, p = 0.75)$n
gbarplot(a$n, a$year)
error.bar(a$year, lower = a$lower, upper = a$upper)
mtext("Year", 1, 2.5)
mtext("Number of sections sampled", 2, 2.5)

# Plot of average proportion of retained sections by year:
res <- aggregate(list(p = sm$interval), by = sm[c("region", "date", "transect.name")], function(x) length(unique(x)) / max(x))
res$year <- year(res$date)
a <- aggregate(list(n = res$p), by = res["year"], mean)
a$n <- round(a$n, 3)
a$lower <- aggregate(list(n = res$p), by = res["year"], quantile, p = 0.25)$n
a$upper <- aggregate(list(n = res$p), by = res["year"], quantile, p = 0.75)$n
gbarplot(100 * a$n, a$year)
error.bar(a$year, lower = 100 * a$lower, upper = 100 * a$upper)
mtext("Year", 1, 2.5)
mtext("Percentage of sections sampled before transect halted", 2, 2.5)
mtext("Proportion of sections retained among sampled sections", 3, 0.5, cex = 1.25)

# Start date variations by region:
years <- min(year(unique(x$section$date))):max(year(unique(x$section$date)))
regions <- sort(unique(x$sections$region))
plot(range(years), c(170, 240), type = "n", xlab = "", ylab = "")
grid()
tab <- unique(x$sections[c("region", "date")])
tab$julian <- julian(as.POSIXct(tab$date))
regions <- sort(unique(x$sections$region))
for (i in 1:length(regions)){
   tt <- tab[tab$region == regions[i], ]
   tt$year <- year(tt$date)
   tt <- aggregate(tt["julian"], by = tt["year"], mean)
   lines(tt$year, tt$julian, col = rainbow(length(regions))[i], lwd = 2)
}
legend("topleft", legend = regions, lwd = 2, col = rainbow(length(regions)))
mtext("Day of year", 2, 2.5, cex = 1.5)
mtext("Year", 1, 2.5, cex = 1.5)

#
res <- aggregate(x$sections["interval"], by = x$sections[c("region", "date", "transect.name")], max)

# Last digit plot (global):
language <- language("en")
tiff(file = paste0("figures/CL last digit - ", language, ".tiff"), compression = "lzw", units = "in", res = 300, height = 5.0, width = 7.5)
m <- rbind(0, cbind(0, 0, matrix(1, nrow = 10, ncol = 10), 0), 0, 0)
layout(m)
par(mar = c(0,0,0,0))
y <- as.numeric(substr(x$observations$carapace.length.mm, nchar(x$observations$carapace.length.mm), nchar(x$observations$carapace.length.mm)))
gbarplot(table(y), grid = TRUE, xaxt = "n")
axis(1, at = seq(0, 10, by = 2), cex.axis = 1.25)
axis(1, at = seq(1, 10, by = 2), cex.axis = 1.25)
mtext("Last digit of carapace length", 1, 3.0, cex = 1.25)
mtext("Frequency", 2, 3, cex = 1.25)
box(col = "grey50")
dev.off()

# Last digit plot (by year):
language <- language("en")
tiff(file = paste0("figures/CL last digit by year - ", language, ".tiff"), compression = "lzw", units = "in", res = 300, height = 5.0, width = 7.5)
m <- rbind(0, cbind(0, 0, matrix(1, nrow = 10, ncol = 10), 0), 0, 0)
layout(m)
par(mar = c(0,0,0,0))
ix <- which(x$observations$certainty.rating == 1)
n <- nchar(x$observations$carapace.length.mm)
y <- as.numeric(substr(x$observations$carapace.length.mm[ix], n[ix], n[ix]))
year <- year(x$observations$date[ix])
tab <- table(year, y)
tab <- tab / repvec(apply(tab, 1, sum), ncol = 10)
image(1:nrow(tab), (1:ncol(tab))-1, tab, col = hcl.colors(20, "YlOrRd", rev = TRUE), xaxt = "n", yaxt = "n", xlab = "", ylab = "")
axis(1, at = 1:nrow(tab), labels = rownames(tab), cex.axis = 1.25)
axis(2, at = seq(0, 10, by = 2), cex.axis = 1.25)
axis(2, at = seq(1, 10, by = 2), cex.axis = 1.25)
hline((0:9) + 0.5, col = "grey70")
vline((1:nrow(tab)) + 0.5, col = "grey70")
mtext("Last digit of carapace length", 2, 3.0, cex = 1.25)
mtext("Year", 1, 3.5, cex = 1.25)
box(col = "grey50")
dev.off()

# Last digit plot (by diver):
language <- language("en")
tiff(file = paste0("figures/CL last digit by diver - ", language, ".tiff"), compression = "lzw", units = "in", res = 300, height = 5.0, width = 7.5)
m <- rbind(0, cbind(0, 0, matrix(1, nrow = 10, ncol = 10), 0), 0, 0, 0, 0, 0)
layout(m)
par(mar = c(0,0,0,0))
ix <- which(x$observations$certainty.rating == 1)
n <- nchar(x$observations$carapace.length.mm)
y <- as.numeric(substr(x$observations$carapace.length.mm[ix], n[ix], n[ix]))
diver <- x$observations$diver[ix]
tab <- table(diver, y)
tab <- tab / repvec(apply(tab, 1, sum), ncol = 10)
image(1:nrow(tab), (1:ncol(tab))-1, tab, col = hcl.colors(20, "YlOrRd", rev = TRUE), xaxt = "n", yaxt = "n", xlab = "", ylab = "")
axis(1, at = 1:nrow(tab), labels = rownames(tab), cex.axis = 1.25, las = 2)
axis(2, at = seq(0, 10, by = 2), cex.axis = 1.25)
axis(2, at = seq(1, 10, by = 2), cex.axis = 1.25)
hline((0:9) + 0.5, col = "grey70")
vline((1:nrow(tab)) + 0.5, col = "grey70")
mtext("Last digit of carapace length", 2, 3.0, cex = 1.25)
mtext("Diver", 1, 9.5, cex = 1.25)
box(col = "grey50")
dev.off()
