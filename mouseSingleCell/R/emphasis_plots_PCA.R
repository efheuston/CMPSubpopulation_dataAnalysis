#' emphasis_plots_PCA
#'
#' @param object seurat object
#' @param file_suffix name to add to end of file
#' @param iterationSlot what to iterate over (i.e., dimensions)
#' @param background_color color to give non-emphasized pops
#' @param emphasis_color color to give emphasized pops
#' @param ... parameters passed on to png_plotFunction
#'
#' @return
#'
emphasis_plots_PCA<-function(object, file_suffix, iterationSlot, background_color = "black", emphasis_color = "red", ...){
	iteration_vector<-(unique(FetchData(object = object, vars.all = iterationSlot)[,1]))
	for(x in 1:length(iteration_vector)){
		color_vector <-c(rep(background_color, length(iteration_vector)))
		color_vector[x]<-emphasis_color
		unique_id<-iteration_vector[x]
		png_plotFunction(PCAPlot(object = object, cols.use = color_vector, group.by = iterationSlot, pt.size = 1), filename = paste(option_arguments$projectName, "_", unique_id,"_", file_suffix,".png", sep = ""), ...)
	}
}
