#' drawmyplot
#'
#' @param geneList List of genes on which to base violin and feature plots
#' @param tsne.obj Seurat2 object
#' @param name FilePrefex
#' @importFrom Seurat VlnPlot FeaturePlot
#'
#' @return
#'
#' @examples
drawmyplot<-function(geneList, tsne.obj, name){
	for(gene in geneList){
		png(filename=paste("ViolinPlots/", name, "-Vln_",gene, ".png", sep=""), width=800, heigh=800)
		try(print(VlnPlot(tsne.obj, gene, do.return = T, point.size.use = 0, cols.use = my_palette)))
		dev.off()
		png(filename=paste("FeaturePlots/", name, "-Ftr_",gene, ".png", sep=""), width=800, height=800)
		try(FeaturePlot(tsne.obj, gene,cols.use=c("grey", "blue"), pt.size = 1, do.return = T))
		dev.off()
	}
}
