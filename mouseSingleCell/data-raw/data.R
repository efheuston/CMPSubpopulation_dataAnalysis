# data-raw/data.R
# Import data lists

#' @title mouseIDTable
#' @description Local table containing mouse gene identifier conversion information. Table generated 2021.04.20 with bioconductor biomaRt
#' @format A data frame with 24117 rows and 5 columns
#' \describe{
#' \item{\code{Mouse.Gene_ID}}{double mouse ensembl gene ID}
#' \item{\code{MGI}}{double MGI gene ID}
#' \item{\code{Entrez.Gene_ID}}{double entrez gene ID}
#' \item{\code{Human.Gene_ID}}{double human homologue ensembl gene ID}
#' \item{\code{HGNC}}{double HGNC gene ID}
#' }
#' @source \url{http://www.ensembl.org}

# retrieve mouseIDTable path
mouseIDTable.file <- system.file(
	"extdata",
	"mouseIDTable.txt",
	package = "mouseSingleCell"
)


mouseIDTable <- read.csv(mouseIDTable.file, header = TRUE, sep = "\t")

usethis::use_data(mouseIDTable, overwrite = TRUE)
