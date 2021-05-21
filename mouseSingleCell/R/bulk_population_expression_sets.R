#' bulk_population_expression_sets
#'
#' @description generate expression lists for gene set enrichment analysis
#'
#' @param file.name name of gmt file
#' @param bulk.contrasts set of contrasts to make
#'
#' @import edgeR
#'
#' @return
#' @export
#'
bulk_population_expression_sets <- function(bulk.contrasts, file.name){
	if(!grepl(pattern = 'gmt$', file.name)){
		file.name <- paste0(file.name, '.gmt')
	}
	bulk.list <- list()
	for(i in 1:(ncol(bulk.contrasts))){
		outname <- names(bulk.contrasts[1,])[i]
		print(outname)
		lrt <- edgeR::glmLRT(fit, contrast = bulk.contrasts[,i])
		# Create rankings
		all.topTags <- topTags(lrt, n = nrow(lrt$table))$table
		all.topTags <- all.topTags[all.topTags$FDR <= 0.05,]
		all.topTags$ranking <- all.topTags$logFC * -log10(all.topTags$FDR)
		bulk.list[[outname]]<-all.topTags$Symbol
		# data[,ncol(data) + 1] <- new
	}
	cell.types <- lapply(sapply(names(mylist), strsplit, split = 'v'), `[[`,1) # extract first named cell type from each comparison
	cell.types <- unique(unname(unlist(cell.types))) # get unique list of cell types

	gene.list <- list()
	for(cell.pop in cell.types){
		cell.set <- list()
		for(cell.comparison in names(mylist)){
			# extract list names to be compared
			if(grepl(paste0('^', cell.pop, 'v'), cell.comparison) == TRUE){
				cell.set <- c(cell.set, toString(cell.comparison))
			}
		}
		cell.set <- unlist(cell.set)
		gene.set <- mylist[[cell.set[1]]]
		if(length(cell.set) >1){
			for(i in 2:length(cell.set)){
				gene.set <- intersect(gene.set, mylist[[cell.set[i]]])
			}
		}
		gene.list <- c(gene.list, cell.pop = list(gene.set))
	}

	names(gene.list) <- cell.types
	gene.list <- data.frame(lapply(gene.list, "length<-", max(lengths(gene.list))))
	gene.list <- t(gene.list)
	write.table(gene.list, file = file.name, row.names = TRUE, col.names = FALSE, quote = FALSE, sep = '\t', na = '')
}


