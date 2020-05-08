# Read data file:
x <- read.csv("data/raw/Scuba_Master.csv", fileEncoding = "Windows-1252")

# Chnage a few column names:
str <- names(x)
str[str == "Depth_.ft."] <-"depth"
str[str == "LC_.mm."] <- "length"
str[str == "Année"] <- "year"
str[str == "Jour"] <- "day"
str[str == "Mois"] <- "month"
str[str == "Plongeur"] <- "diver"
names(x) <- str

# Correct diver names:
x$diver <- gsub(" ", "", x$diver) # Remove all spaces.
x$diver[x$diver == "MComeau"] <- "Mcomeau"
x$diver[x$diver == "gpaulin"] <- "Gpaulin"

# Create diver table:
divers <- c("Mcomeau; Michel Comeau", 
            "Bcomeau; Bruno Comeau",
            "Rdoucette; Renelle Doucette",
            "Gpaulin; Gilles Paulin")
divers <- strsplit(divers, "; ")
divers <- data.frame(short = unlist(lapply(divers, function(x) x[1])), 
                     long = unlist(lapply(divers, function(x) x[2])), 
                     stringsAsFactors = FALSE)
   
# Replace divers with full names:
index <- match(x$diver, divers$short)
ii <- which(!is.na(index))
x$diver[ii] <- divers$long[index[ii]]

# Transect table:
vars <- c("year", "Region", "Transect")
transects <- unique(x[vars])



