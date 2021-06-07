#' create_multi_object_list
#'
#' @param x directory containing barcode, gene, and matrix files along the path /outs/filtered_feature_bc_matrix/
#' @import Seurat
#' @return a multiobject list
#' @export
#'
#'
#'
#'



create_multi_object_list<-function(x){
	if(!any(grepl(pattern="gz$", x = list.files(paste(x, "/outs/filtered_feature_bc_matrix/", sep = ""))))){
		print("Unzipping cellranger zips")
		cellranger_zips<-paste(x, "/outs/filtered_feature_bc_matrix/",
													 grep(pattern="gz$", x = list.files(paste(x, "/outs/filtered_feature_bc_matrix/", sep = "")), value = TRUE),
													 sep="")
		lapply(cellranger_zips, gunzip, remove=FALSE)
		file.rename(from = paste(x, "/outs/filtered_feature_bc_matrix/features.tsv", sep=""), to = paste(x, "/outs/filtered_feature_bc_matrix/genes.tsv", sep=""))
	}

	print(paste("reading", paste(x, "/outs/filtered_feature_bc_matrix/", sep=""), sep = " "))
	cellranger_files<-paste(x, "/outs/filtered_feature_bc_matrix/", sep="")
	my_object.data<-Read10X(cellranger_files)
	my_object<- CreateSeuratObject(raw.data = my_object.data, min.cells = 3, min.genes = 200, project = projectName)
	my_object<-RenameCells(my_object, add.cell.id = x)
	print(head(my_object@cell.names))
	my_object<- create_percentMito_column(my_object)
	my_object<-NormalizeData(my_object, normalization.method = "LogNormalize", scale.factor = 10000)
	my_object<-FindVariableGenes(my_object,mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
	my_object<-ScaleData(my_object, vars.to.regress = vars_to_regress, do.par = TRUE, num.cores = future::availableCores())
	my_object<-CellCycleScoring(my_object, s.genes = s.genes, g2m.genes = g2m.genes, set.ident = FALSE)
	my_object@meta.data$orig.ident<-x

	multi_object_list$x<-my_object

}

