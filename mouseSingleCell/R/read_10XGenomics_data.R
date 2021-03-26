#' read_10XGenomics_data
#'
#' @param x directory containing gene, barcode, and matrix files along the path /outs/filtered_gene_bc_matrices/<genome>
#'
#' @return
#'
#' @examples
read_10XGenomics_data<-function(x){
	print(paste("reading", paste(x, "/outs/filtered_gene_bc_matrices/",option_arguments$genome,"/", sep=""), sep = " "))
	x = paste(x, "/outs/filtered_gene_bc_matrices/",option_arguments$genome,"/", sep="")
}

