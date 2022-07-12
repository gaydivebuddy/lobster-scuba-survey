# Year effect:
gbarplot(exp(ranef[rownames(ranef) == "year_effect", 1]))

# Diver effect:
diver_effect <- exp(ranef[rownames(ranef) == "diver_effect", 1])
names(diver_effect) <- divers
gbarplot(diver_effect, grid = TRUE)
hline(1, col = "red", lwd = 2)
text(1:length(diver_effect), 0.5 * diver_effect, names(diver_effect), srt = 90, cex = 1.25)
mtext("Diver", 1, 2.75, cex = 1.5)
mtext("Diver effect", 2, 2.75, cex = 1.5)

# Transect effect:
gbarplot(exp(ranef[rownames(ranef) == "transect_effect", 1]), grid = TRUE)

# Region effect:
region_effect <- exp(ranef[grep("region_effect", rownames(ranef)), 1])[match(regions, regions.ordered)]
names(region_effect) <- regions.ordered
gbarplot(region_effect)
hline(1, col = "red", lwd = 2)
text(1:length(regions), 0.5 * region_effect, names(region_effect), srt = 90)

# Length x diver effect effect:
length_diver <- ranef[rownames(ranef) == "length_diver_effect", 1]
dim(length_diver) <- c(max(data$diver)+1, max(data$len)+1)
dimnames(length_diver) <- list(divers, lens)
length_diver <- length_diver[setdiff(rownames(length_diver), c("Asselin, Natalie", "Leblanc, Stepan")), ]
image(1:nrow(length_diver), 1:ncol(length_diver), length_diver, col = colorRampPalette(c("blue", "white", "red"))(10),
      xaxt = "n", yaxt = "n", zlim = c(-0.5, 0.5))
axis(1, at = 1:nrow(length_diver), labels = rownames(length_diver), las = 2)
axis(2, at = 1:ncol(length_diver), labels = colnames(length_diver), las = 2)

# Diver x year effect:
diver_year_effect <- ranef[rownames(ranef) == "diver_year_effect", 1]
dim(diver_year_effect) <- c(max(data$year)+1, max(data$diver)+1)
dimnames(diver_year_effect) <- list(min(year(mm$date)):max(year(mm$date)), divers)
m <- cbind(matrix(0, nrow = 7, ncol = 1), matrix(1, nrow = 7, ncol = 7))
m <- cbind(0, rbind(0, m, 0), 0)
layout(m)
par(mar = c(0,0,0,0))
image(1:nrow(diver_year_effect), 1:ncol(diver_year_effect), diver_year_effect, 
      col = colorRampPalette(c("blue", "white", "red"))(10),
      xaxt = "n", yaxt = "n", zlim = c(-0.1, 0.1))
axis(1, at = 1:nrow(diver_year_effect), labels = rownames(diver_year_effect), las = 2)
axis(2, at = 1:ncol(diver_year_effect), labels = colnames(diver_year_effect), las = 2)

# Length x year effect:
language <- language("en")
tiff(file = paste0("figures/length x year effect - ", language, ".tiff"), compression = "lzw", units = "in", res = 300, height = 10, width = 7.5)
m <- rbind(0, cbind(0, 0, matrix(1, nrow = 10, ncol = 10), 0), 0)
layout(m)
par(mar = c(0,0,0,0))
length_year <- as.numeric(ranef[rownames(ranef) == "length_year_effect", 1])
dim(length_year) <- c(max(data$year)+1, max(data$len)+1)
dimnames(length_year) <- list(year = min(year(mm$date)):max(year(mm$date)), length = lens)
image(1:nrow(length_year), 1:ncol(length_year), length_year, 
      col = colorRampPalette(c("blue", "white", "red"))(20),
      xaxt = "n", yaxt = "n", zlim = c(-1.50, 1.50))
axis(1, at = 1:nrow(length_year), labels = rownames(length_year), las = 2, cex.axis = 1.25)
axis(2, at = seq(0, 120, by = 10), cex.axis = 1.25)
mtext("Year", 1, 4.5, cex = 1.25)
mtext("Carapace length (mm)", 2, 3, cex = 1.25)
box(col = "grey50")
dev.off()

# Length x region effect:
language <- language("en")
tiff(file = paste0("figures/length x region effect - ", language, ".tiff"), compression = "lzw", units = "in", res = 300, height = 10, width = 7.5)
m <- rbind(0, cbind(0, 0, matrix(1, nrow = 10, ncol = 10), 0), 0, 0)
layout(m)
par(mar = c(0,0,0,0))
length_region <- as.numeric(ranef[rownames(ranef) == "length_region_effect", 1])
dim(length_region) <- c(max(data$region)+1, max(data$len)+1)
dimnames(length_region) <- list(region = regions, length = lens)
image(1:nrow(length_region), 1:ncol(length_region), length_region, 
      col = colorRampPalette(c("blue", "white", "red"))(20),
      xaxt = "n", yaxt = "n", zlim = c(-1.6, 1.6))
axis(1, at = 1:nrow(length_region), labels = rownames(length_region), las = 2, cex.axis = 1.25)
axis(2, at = seq(0, 120, by = 10), cex.axis = 1.25)
mtext("Region", 1, 8, cex = 1.25)
mtext("Carapace length (mm)", 2, 3, cex = 1.25)
box(col = "grey50")
dev.off()

# Year x region effect:
# region_year_effect[region[i] * n_year + year[i]] +
tiff(file = paste0("figures/region x year effect - ", language, ".tiff"), compression = "lzw", units = "in", res = 300, height = 5, width = 7.5)
m <- cbind(matrix(0, nrow = 7, ncol = 1), matrix(1, nrow = 7, ncol = 7))
m <- cbind(0, rbind(0, m, 0, 0), 0)
layout(m)
par(mar = c(0,0,0,0))
region_year <- as.numeric(ranef[rownames(ranef) == "region_year_effect", 1])
dim(region_year) <- c(max(data$year)+1, max(data$region)+1)
dimnames(region_year) <- list(year = min(year(mm$date)):max(year(mm$date)), region = regions)

image(1:nrow(region_year), 1:ncol(region_year), region_year, 
      col = colorRampPalette(c("blue", "white", "red"))(20),
      xaxt = "n", yaxt = "n", zlim = c(-3.50, 3.50))
axis(1, at = 1:nrow(region_year), labels = rownames(region_year), las = 2, cex.axis = 1.25)
axis(2, at = 1:ncol(region_year), labels = colnames(region_year), las = 2, cex.axis = 1.25)
mtext("Year", 1, 4.5, cex = 1.25)
mtext("Region", 2, 8.5, cex = 1.25)
box(col = "grey50")
dev.off()

# Predictions:
mu <- obj$report()$mu
dimnames(mu) <- list(length = lens, 
                     year = years,
                     region = regions)

bounds <- list("Caraquet" = c(0, 6, 1),
               "Cocagne" = c(0, 4, 0.5),
               "Fox Harbor" = c(0, 0.02, 0.005), 
               "Murray Corner" = c(0, 0.5, 0.1),
               "Neguac" = c(0, 1.5, 0.3),
               "Pointe-Verte" = c(0, 4, 0.5),
               "Richibucto" = c(0, 3.0, 0.5),
               "Shediac" = c(0, 1.5, 0.3),
               "Toney River" = c(0, 0.25, 0.05))
for (j in 1:length(regions)){
   tiff(file = paste0("figures/length frequencies ", regions[j], " - ", language, ".tiff"), compression = "lzw", units = "in", res = 300, height = 10, width = 7.5)
   m <- kronecker(matrix(1:12, ncol = 2), matrix(1, nrow = 5, ncol = 5))
   m <- cbind(0, rbind(0, m, 0, 0), 0)
   layout(m)
   par(mar = c(0,0,0,0))
   for (i in 8:length(years)){
      ylim = bounds[regions[j]][[1]]
      gbarplot(100*mu[,as.character(years[i]), regions[j]], 
               xlim = c(0, 100), ylim = ylim[1:2], xaxs = "i", yaxs = "i", 
               xaxt = "n", yaxt = "n", grid = TRUE)
      text(par("usr")[1] + 0.8 * diff(par("usr")[1:2]),
           par("usr")[3] + 0.8 * diff(par("usr")[3:4]),
           years[i], cex = 1.25)
      
      if ((i-7) %in% 3) mtext("Density (# /m2)", 2, 3.0, at = 0)
      if ((i-7) %in% 1) axis(2, at = seq(ylim[1], ylim[2], by = ylim[3]))
      if ((i-7) %in% 2:6) axis(2, at = seq(ylim[1], ylim[2]-ylim[3], by = ylim[3]))
      if ((i-7) == 6)  axis(1, at = 10 * seq(0, 8, by = 2))
      if ((i-7) == 12){
         mtext("Carapace length (mm)", 1, 3.0, at = 0)
         axis(1, at = 10 * seq(0, 10, by = 2))
      } 
      #if ((i-7) == 7)  mtext(regions[j], 3, 0.5, cex = 1.5)
      box(col = "grey60")
   }
   dev.off()
}


# Variance parameters:
tiff(file = paste0("figures/variance parameters - ", language, ".tiff"), compression = "lzw", units = "in", res = 300, height = 7, width = 7.5)
vars <- theta[grep("log_sigma_", names(theta))]
gbarplot(exp(vars), grid = TRUE, xaxt = "n", ylim = c(0, 2.25))
text((1:length(vars))-0.25, 1.000 * exp(vars) + 0.02, gsub("_", " x ", gsub("log_sigma_", "", names(vars))), 
     srt = 90, pos = 4, font = 2)
#axis(1, at = 1:length(vars), labels = gsub("log_sigma_", "", names(vars)), las = 2)
mtext("Random effect", 1, 1.0, cex = 1.25)
mtext("Effect scale", 2, 2.5, cex = 1.25)
box(col = "grey50")
dev.off()

# Plot indices:
language <- language("en")
tiff(file = paste0("figures/abundance 0-40mm CL - ", language, ".tiff"), compression = "lzw", units = "in", res = 300, height = 3.75, width = 7.5)
m <- rbind(0, cbind(0, 0, matrix(1, nrow = 10, ncol = 10), 0), 0, 0)
layout(m)
par(mar = c(0,0,0,0))
L40 <- obj$report()$abundance_L40
dimnames(L40) <- list(year = years, region = regions)
image(1:nrow(L40), 1:ncol(L40), L40, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
axis(1, at = 1:nrow(L40), labels = rownames(L40))
axis(2, at = 1:ncol(L40), labels = colnames(L40), las = 2)
hline((1:ncol(COM))-0.5, col = "grey80") 
vline((1:nrow(COM))-0.5, col = "grey80") 
if (language == "english"){
   mtext("Year", 1, 3.25 , cex = 1.5)
}else{
   mtext("Année", 1, 3.25 , cex = 1.5)
}
box(col = "grey70")
dev.off()

# Plot indices:
language <- language("en")
tiff(file = paste0("figures/abundance 20-40mm CL - ", language, ".tiff"), compression = "lzw", units = "in", res = 300, height = 3.75, width = 7.5)
m <- rbind(0, cbind(0, 0, matrix(1, nrow = 10, ncol = 10), 0), 0, 0)
layout(m)
par(mar = c(0,0,0,0))
G20L40 <- obj$report()$abundance_G20L40
dimnames(G20L40) <- list(year = years, region = regions)
image(1:nrow(G20L40), 1:ncol(G20L40), G20L40, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
axis(1, at = 1:nrow(G20L40), labels = rownames(G20L40))
axis(2, at = 1:ncol(G20L40), labels = colnames(G20L40), las = 2)
hline((1:ncol(COM))-0.5, col = "grey80") 
vline((1:nrow(COM))-0.5, col = "grey80") 
if (language == "english"){
   mtext("Year", 1, 3.25 , cex = 1.5)
}else{
   mtext("Année", 1, 3.25 , cex = 1.5)
}
box(col = "grey70")
dev.off()

# Plot indices:
language <- language("en")
tiff(file = paste0("figures/abundance commercial - ", language, ".tiff"), compression = "lzw", units = "in", res = 300, height = 3.75, width = 7.5)
m <- rbind(0, cbind(0, 0, matrix(1, nrow = 10, ncol = 10), 0), 0, 0)
layout(m)
par(mar = c(0,0,0,0))
COM <- obj$report()$com
dimnames(COM) <- list(year = years, region = regions)
image(1:nrow(COM), 1:ncol(COM), COM, xlab = "", ylab = "", xaxt = "n", yaxt = "n")
axis(1, at = 1:nrow(COM), labels = rownames(COM))
axis(2, at = 1:ncol(COM), labels = colnames(COM), las = 2)
hline((1:ncol(COM))-0.5, col = "grey80") 
vline((1:nrow(COM))-0.5, col = "grey80") 
if (language == "english"){
   mtext("Year", 1, 3.25 , cex = 1.5)
}else{
   mtext("Année", 1, 3.25 , cex = 1.5)
}
box(col = "grey70")
dev.off()
