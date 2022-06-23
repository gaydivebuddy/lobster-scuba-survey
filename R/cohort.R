cohort <- function(x, year, region){
   # Collate to data frame:
   x <- data.frame(length = x,
                   year = year, 
                   region = region)
   
   x$region <- tolower(x$region)
   res <- rep(NA, nrow(x))
   
   regions <- sort(unique(x$region))
   years <- sort(unique(x$year))
   t <- read.table("data/cohort.bounds.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)
   for (i in 1:length(regions)){
      for (j in 1:length(years)){
         index <- (x$region == regions[i]) & (x$year == years[j])
         if (sum(index) > 0){
            tt <- t[(t$region == regions[i]) & (t$year == years[j]), ]
            if (nrow(tt) == 0){
               if (regions[i] %in% c("Cocagne", "Fox Harbor", "Murray Corner", "Robichaud", "Toney River", "Shediac")) tt <- t[(t$region == "Shediac") & (t$year == years[j]), ]
               if (regions[i] %in% c("Caraquet", "Anse Bleue", "Grande-Anse", "Neguac", "Pigeon Hill", "Pointe-Verte", "Tabusintac")) tt <- t[(t$region == "Caraquet") & (t$year == years[j]), ]
            }
            if (nrow(tt) == 0) tt <- t[(t$region == "Caraquet") & (t$year == years[j]), ]  
            print(years[j])
            print(regions[i])
            print(tt)
            xx <- repvec(x$length[index], ncol = nrow(tt))
            ll <- repvec(tt$lower, nrow = nrow(xx))
            uu <- repvec(tt$upper, nrow = nrow(xx))
            fun <- function(x){
               index <- which(x)
               if (length(index) == 0) index <- NA
               return(index)
            } 
            res[index] <- apply((xx >= ll) & (xx <= uu), 1, fun) - 1
         }
      }
   }
   
   return(res)
}

