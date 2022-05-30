# Read transect section table:
load("section.rdata")
s$year <- year(s$date)

# Load biological data:
load("biological.rdata")
b <- b[b$species == "American lobster", ]

delete <- c('id', "species", 'created_at', "section_id", 'updated_at', 'created_by', 'updated_by', 'created_by_id', 'updated_by_id', "species_id", "interval_display", "interval", "sample_id")
names(b) <- gsub("side_display", "side", names(b))
names(b) <- gsub("carapace_length_mm", "length", names(b))
names(b) <- gsub("certainty_rating", "certainty", names(b))


vars <- c("sample", "date", "region", "transect", "side", "section", "length", "certainty")

b <- b[, c(vars, setdiff(names(b), c(vars, delete)))]

b <- b[year(b$date) == 2021, ]
