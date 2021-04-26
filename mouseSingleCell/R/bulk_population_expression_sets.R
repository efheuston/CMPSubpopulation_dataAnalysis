#' bulk_population_expression_sets
#'
#' @description generate expression lists for gene set enrichment analysis
#'
#' @param bulk.contrasts set of contrasts to make
#' @import edgeR
#'
#' @return
#' @export
#'
bulk_population_expression_sets <- function(bulk.contrasts){
	for(i in 1:(ncol(bulk.contrasts))){
		outname <- names(bulk.contrasts[1,])[i]
		print(outname)
		lrt <- edgeR::glmLRT(fit, contrast = bulk.contrasts[,i])
		# Create rankings
		all.topTags <- topTags(lrt, n = nrow(lrt$table))$table
		all.topTags <- all.topTags[all.topTags$FDR <= 0.05,]
		all.topTags$ranking <- all.topTags$logFC * -log10(all.topTags$FDR)
		all.topTags<-all.topTags$Symbol
		data[,ncol(data) + 1] <- new
	}
}
