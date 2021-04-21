#' gsea_ranking
#'
#' @description generate ranked subpopulation gene lists
#' @import edgeR
#'
#' @param lrt edgeR lrt object
#' @param outname prefix for generated files
#' @param FDRmax max FDR to allow
#'
#' @return
#' @export
#'
#' @examples
gsea_ranking <- function(lrt, outname, FDRmax = 0.05){
	all.topTags <- topTags(lrt, n = nrow(lrt$table))$table
	all.topTags <- all.topTags[all.topTags$FDR <= FDRmax,]
	all.topTags<-all.topTags[with(all.topTags, order(-logFC)),]
	write.table(all.topTags[, c('HsHmlg', 'logFC')],
							file = paste(outname, '-rankedGSEA.rnk', sep = ""),
							sep = '\t',
							row.names = FALSE,
							col.names = TRUE,
							quote = FALSE,
							na = "not_mapped")
}
