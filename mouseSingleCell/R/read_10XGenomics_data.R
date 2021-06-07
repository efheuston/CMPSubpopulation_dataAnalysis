#' read_10XGenomics_data
#'
#' @param sample.list list of sample names that contain gene, barcode, and matrix files the directory
#'
#' @return named character of sample names and directory paths
#' @export
#'
read_10XGenomics_data<-function(sample.list){
	data <- sapply(sample.list, FUN = function(x) {x = dirname(list.files(path = x, pattern = 'barcodes.tsv', full.names = TRUE, recursive = TRUE, no.. = TRUE))})
	names(data) <- sample.list
	return(data)

}

