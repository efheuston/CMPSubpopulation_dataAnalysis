# data-raw/data.R
# Import data lists

# retrieve mouseIDTable path
mouseIDTable.file <- system.file(
	"extdata",
	"mouseIDTable.txt",
	package = "mouseSingleCell"
)


mouseIDTable <- read.csv(mouseIDTable.file, header = TRUE, sep = "\t")

usethis::use_data(mouseIDTable, overwrite = TRUE)
