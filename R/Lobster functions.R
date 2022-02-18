region.str <- function(x){
   str <- c("AB1", "Anse Bleue 1",
            "AB2", "Anse Bleue 2",
            "CG",  "Cocagne",	
            "CR",  "Caraquet",
            "FH",  "Fox Harbour",	
            "GA1", "Grande-Anse 1",
            "GA2", "Grande-Anse 2",
            "MC",  "Murray Corner",
            "NG",  "Neguac",
            "PH",  "Pigeon Hill",	
            "PV",  "Pointe-Verte",	
            "RI",  "Richibucto",		
            "RO",  "Robichaud",
            "SH",  "Shediac",	
            "TB",  "Tabusintac",	
            "TR",  "Toney River")

   n <- round(length(str) / 2)
   tab <- data.frame(code = str[seq(1, length(str), by = 2)],
                     name = str[seq(2, length(str), by = 2)],
                     stringsAsFactors = FALSE)
                     
   v <- tab$name[match(toupper(x), tab[, 1]) ]
   
   return(v)
}

year.class <- function(x){
   # YEAR.CLASS - Return age class depending on size and year.

   x$region <- tolower(x$region)
   res <- rep(NA, nrow(x))

   regions <- sort(unique(x$region))
   years <- sort(unique(x$year))
   t <- read.table("Data/Year-Class.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)
   for (i in 1:length(regions)){
      for (j in 1:length(years)){
         index <- (x$region == regions[i]) & (x$year == years[j])
         if (sum(index) > 0){
            tt <- t[(t$region == regions[i]) & (t$year == years[j]), ]
            if (nrow(tt) == 0){
               if (regions[i] %in% c("cg", "fh", "mc", "ro", "tr", "sh")) tt <- t[(t$region == "sh") & (t$year == years[j]), ]
               if (regions[i] %in% c("ab1", "ab2", "ga1", "ga2", "ng", "ph", "pv", "tb")) tt <- t[(t$region == "cr") & (t$year == years[j]), ]
            }
            if (nrow(tt) == 0) tt <- t[(t$region == "cr") & (t$year == years[j]), ]  
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

swept.area <- function(x, by = c("year", "month", "region", "site", "transect")){
   # Calculate swept area by specified grouping.

   by <- setdiff(by, "year.class")

   key <- unique(c(by, "year", "month", "region", "site", "transect", "diver", "section"))
   temp <- unique(x[key])
   temp <- aggregate(data.frame(swept.area = temp$year), by = temp[unique(c(by, "year", "region", "site", "transect", "diver", "section")), drop = FALSE], length)
   temp$swept.area <- temp$swept.area * 10
   temp <- aggregate(temp["swept.area"], by = temp[by, drop = FALSE], sum)
   temp <- sort(temp, by = by)

   return(temp)
}

exclude <- function(x){
   # Add Michel's exclusions:
   
   file <- paste("Data/", "CR exclusion.csv", sep = "")

   t <- read.table(file , sep = ",", header = TRUE, stringsAsFactors = FALSE)
   names(t) <- tolower(names(t))
   temp <- names(t)
   index <- grep("x", names(t))
   temp[index] <- gsub("x", "n", temp[index])
   names(t) <- temp
   t$site <- gsub("-", "", t$site)
   
   ivars <- c("year", "region", "site", "transect")
   vars <- setdiff(names(t), ivars)
   
   if (!all(names(t) %in% names(x))) stop("Some variables missing to apply exclusions.")
   
   temp <- x
   temp[vars] <- NA
   temp <- gulf::merge(temp, t, by = c("year", "region", "site", "transect"), all.x = TRUE, overwrite = TRUE)
   temp[is.na(temp)] <- 0
   index <- (temp[vars] == 1)

   temp <- x[vars]
   temp[index] <- NA
   x[vars] <- temp
   
   return(x)
}

frequency <- function(x, index = NULL, by = c("year", "region", "site", "transect"),
                      sex = NULL, year.class = NULL, berried = NULL,
                      region = "cr"){
   # FREQUENCY - Calculate observed frequencies by specified category.

   # Parse category arguements:
   if (!is.null(sex)) sex <- tolower(sex)
   if (!is.null(berried)) berried <- tolower(berried)

   # Initialize index variable:
   if (is.null(index)) index <- !is.na(x$length)

   # Build index variable:
   if (!is.null(sex)) index <- index & (data$sex %in% sex)
   if (!is.null(berried)) index <- index & (data$berried %in% berried)
   if (!is.null(year.class)) index <- index & (year.class(data) %in% year.class)

   # Define 'year.class' as a variable:
   if ("year.class" %in% by){
      x$year.class <- year.class(x, region = region)
   }

   # Calculate number of specimens using 'by':
   temp <- aggregate(data.frame(n = x$year[index]), by = x[index, by, drop = FALSE], length)

   # Merge results into swept area table:
   temp <- merge(swept.area(x, by = by), temp, by = by, all.x = TRUE)
   temp$n[is.na(temp$n)] <- 0
   temp$n[is.na(temp$swept.area)] <- NA
   temp <- sort(temp,  by = by)

   return(temp)
}
