#' runMonocleAnalysis
#'
#' @param args list of arguments passed to runMonocleAnalysis function (seurat_object_filename, num_dim, cca_variables)
#'
#' @return
#' @export
#' @import Seurat
#' @import monocle
#' @importFrom future availableCores
#' @import R.utils
#'

runMonocleAnalysis <- function(args){

	# Process commandLine input -----------------------------------------------

	option_arguments<-read_monocle_args(args)
	option_arguments$num_dim<-as.numeric(option_arguments$num_dim)
	option_arguments$max_components<-as.numeric(option_arguments$max_components)


	print(paste("Setting UMI bounded filtering to", option_arguments$UMI_bounded_filtering, sep = " "))
	print(paste("Setting correction variables as", paste(option_arguments$cca_variables, collapse = " "), sep = " "))
	print(paste("Will process", option_arguments$seurat_object_filename, "using", option_arguments$num_dim, "dimensions and", option_arguments$max_components, "components", sep = " "))

	output_prefix<-gsub(x=option_arguments$seurat_object_filename, pattern = ".rds", replacement = "")
	output_prefix<-gsub(x = output_prefix, pattern =  "_tsne", replacement =  "")



	# Load Seurat object as CDS -----------------------------------------------

	if(grepl(x = option_arguments$seurat_object_filename, pattern = "rds$", ignore.case = T)){
		seurat_object<-readRDS(option_arguments$seurat_object_filename)
	} else if (grepl(x = option_arguments$seurat_object_filename, pattern = "Robj$", ignore.case = T)){
		seurat_object<-get(load(option_arguments$seurat_object_filename))
	} else {
		print("I don't recognize this type of file")
	}

	seurat_varGenes<-seurat_object@var.genes

	if(option_arguments$color_by_seurat_res == TRUE){
		seurat_res<-colnames(seurat_object@meta.data)[grep(x=colnames(seurat_object@meta.data), pattern = "^res")]
	}



	# Initialize monocle object -----------------------------------------------

	monocle_object<-importCDS(seurat_object, import_all = FALSE)


	if(option_arguments$color_by_seurat_res == TRUE){
		adjust_palette_size(length(unique(seurat_object@ident)), basic_color_palette = basic_color_palette)
	} else{
		basic_color_palette = basic_color_palette
	}

	remove(seurat_object)
	gc()

	monocle_object<-estimateSizeFactors(monocle_object)
	monocle_object<-estimateDispersions(monocle_object)


	# Filter low-quality cells ------------------------------------------------

	if(option_arguments$perform_expression_filtering == TRUE){
		monocle_object<-detectGenes(monocle_object, min_expr = 0.1)
		expressed_genes<-row.names(subset(fData(monocle_object), num_cells_expressed >=10))
	}


	pData(monocle_object)$Total_nUMI<-Matrix::colSums(exprs(monocle_object))
	monocle_object<-monocle_object[,pData(monocle_object)$Total_nUMI < 1e6]

	lower_bound<-10^(mean(log10(pData(monocle_object)$Total_nUMI)) - 2*sd(log10(pData(monocle_object)$Total_nUMI)))
	upper_bound<-10^(mean(log10(pData(monocle_object)$Total_nUMI)) + 2*sd(log10(pData(monocle_object)$Total_nUMI)))

	# Plot UMI distributions and filter --------------------------------------------------

	png_plotFunction(qplot(Total_nUMI, data = pData(monocle_object), color = orig.ident, geom = "density") +
									 	geom_vline(xintercept = upper_bound) +
									 	geom_vline(xintercept = lower_bound),
									 filename = paste(output_prefix, "_nUMI_2SD-byID.png", sep = ""))

	try(
		if(option_arguments$color_by_seurat_res == TRUE){
			lapply(seurat_res[1:length(seurat_res)], plotting_function = "qplot", the_object = monocle_object, cycle_plot_param)
		}
	)

	if(grep("upper", option_arguments$UMI_bounded_filtering, ignore.case = TRUE)){
		monocle_object<-monocle_object[,pData(monocle_object)$Total_nUMI < upper_bound]
	} else if(grep("both", option_arguments$UMI_bounded_filtering, ignore.case = TRUE)){
		monocle_object<-monocle_object[,pData(monocle_object)$Total_nUMI < upper_bound & pData(monocle_object)$Total_nUMI > lower_bound]
	} else if(is.null(option_arguments$UMI_bounded_filtering) | grep("none", option_arguments$UMI_bounded_filtering, ignore.case = TRUE)){
		print("Not filtering based on UMI counts")
	} else {
		print("No boundries requested")
	}

	monocle_object<-detectGenes(monocle_object, min_expr = 0.1)



	# Cluster cells without marker genes --------------------------------------

	disp_table<-dispersionTable(monocle_object)

	monocle_clustering_genes<-subset(disp_table, mean_expression >= 0.1)
	monocle_object<-setOrderingFilter(monocle_object, ordering_genes = monocle_clustering_genes$gene_id)

	png_plotFunction(plot_ordering_genes(monocle_object), filename = paste(output_prefix, "_PlotOrderingGenes.png", sep = ""))
	png_plotFunction(plot_pc_variance_explained(monocle_object, return_all = FALSE), filename = paste(output_prefix, "_PCVarXplained.png", sep = ""))

	# Perform clustering ------------------------------------------------------

	monocle_object<-reduceDimension(monocle_object,
																	max_components = option_arguments$max_components,
																	num_dim = as.numeric(option_arguments$num_dim),
																	reduction_method = 'tSNE',
																	residualModelFormulaStr = paste(option_arguments$cca_variables, collapse = " "),
																	verbose = TRUE)
	monocle_object<-clusterCells(monocle_object)

	png_plotFunction(plot_cell_clusters(monocle_object, 1, 2, color = "orig.ident"),
									 filename = paste(output_prefix, "_clstr-orig.ident.png", sep = ""),
									 height = 1600,
									 width = 1600)


	try(
		if(option_arguments$color_by_seurat_res == TRUE){
			lapply(seurat_res[cycle_parameter = 1:length(seurat_res)], plotting_function = "cluster", the_object = monocle_object, cycle_plot_param)
		}
	)

	save.image(file=paste(output_prefix, ".RData", sep = ""))


	# Perform trajectory analysis with monocle-defined genes ------------------

	clustering_DEG_monocle<- differentialGeneTest(monocle_object[expressed_genes,], fullModelFormulaStr='~Cluster')
	#current default is to take top 1000 clustering_DEG_monocle genes

	monocleBased_orderedgenes<-row.names(clustering_DEG_monocle)[order(clustering_DEG_monocle$qval)][1:1000]
	monocle_object<-setOrderingFilter(monocle_object, ordering_genes=monocleBased_orderedgenes)
	monocle_object<-reduceDimension(monocle_object, method="DDRTree")
	monocle_object<-orderCells(monocle_object)



	# Plot cell trajectories --------------------------------------------------



	try(
		png_plotFunction(my_plot_cell_trajectory(monocle_object,
																						 cell_size = 1,
																						 point_colors = basic_color_palette),
										 filename = paste(output_prefix, "_trajectory.png", sep = ""))
	)
	try(
		png_plotFunction(my_plot_cell_trajectory(monocle_object,
																						 color_by = "orig.ident",
																						 cell_size = 1,
																						 point_colors = basic_color_palette),
										 filename = paste(output_prefix, "_trajectory-orig.ident.png", sep = ""))
	)
	try(
		png_plotFunction(my_plot_cell_trajectory(monocle_object,
																						 color_by = "orig.ident",
																						 cell_size = 1,
																						 point_colors = basic_color_palette) +
										 	facet_wrap(~orig.ident),
										 filename = paste(output_prefix, "_trajectory-orig.identFACET.png", sep = ""))
	)

	try(
		if(option_arguments$cca_variables == TRUE){
			lapply(seurat_res[cycle_parameter = 1:length(seurat_res)], plotting_function = "trajectory", the_object = monocle_object, cycle_plot_param)
		}
	)


	saveRDS(monocle_object, file = paste(output_prefix, "_UnsupClustMonocle.rds", sep = ""))
	save.image(file=paste(output_prefix, ".RData", sep = ""))


	# Basic differential expression analysis ----------------------------------

	# Hard-coding to use Seurat vargenes as markers... ?
	marker_genes <-row.names(subset(fData(monocle_object), gene_short_name %in% seurat_varGenes))
	diff_test_res <-differentialGeneTest(monocle_object[marker_genes,],
																			 cores = future::availableCores(),
																			 fullModelFormulaStr = "~Cluster") # ASSUMING CLUSTER IS TEH CORRECT FULL MODEL FORMULA STRING
	sig_genes<-subset(diff_test_res, qval < 0.1) # Set for FDR < 10%
	sig_genes[,c("gene_short_name", "pval", "qval")]
	write.table(sig_genes, file = paste(output_prefix, "_differentiallyExpressed_seuratVar.txt", sep = ""), quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)



}
