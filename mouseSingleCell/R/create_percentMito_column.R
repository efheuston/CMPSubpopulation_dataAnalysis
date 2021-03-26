#' create_percentMito_column
#'
#' @param my_object object of Seurat class <=v2.4.3
#'
#' @return
#'
#' @examples
create_percentMito_column<-function(my_object){
	mito.genes<-grep(pattern = "^MT-", x = rownames(x=my_object@data), value = TRUE, ignore.case = TRUE)
	percent.mito <-Matrix::colSums(my_object@raw.data[mito.genes,])/Matrix::colSums(my_object@raw.data)
	my_object<-AddMetaData(my_object, metadata = percent.mito, col.name = "percent.mito")

}
