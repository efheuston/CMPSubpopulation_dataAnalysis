#' read_monocle_args
#'
#' @param args
#'
#' @return
#'
#' @examples
read_monocle_args <-function(args){
	args<-paste(unlist(args), collapse = " ")
	args<-unlist(strsplit(args, "--"))[-1]
	option_arguments<-sapply(args, function(x){
		unlist(strsplit(x, " "))[-1]
	})
	option_names<-sapply(args, function(x){
		unlist(strsplit(x, " "))[1]
	})
	names(option_arguments) <- unlist(option_names)
	required_keys<-c("seurat_object_filename", "num_dim")

	if(!(all(required_keys %in% names(option_arguments)))){
		print(paste("Missing:", paste(required_list[!(required_list%in%names(option_arguments))], collapse = " "), sep = " "))
		stop("Required arguments:
           --seurat_object_filename: name of seurat object to import (note: file must have .Robj or .rds extension)
           --num_dim: dimensionality of the reduced space
           Additional arguments:
           --max_components: Dimensions to plot (default = 2) note: currently I don't think this does anything
           --perform_expression_filtering: should analyzed RNAs be limited to expression of at least 0.1 and in at least 10 cells (default = TRUE)
           --color_by_seurat_res: (default = TRUE)
           --order_by_seurat_varGenes: (default = FALSE)
           --UMI_bounded_filtering: (default = \"upper\", can be \"upper\", \"lower\", \"both\", or \"none\")
           --cca_variables: (default = \"~nUMI + nGene\")")
	} else{

		if(!("max_components" %in% names(option_arguments))){
			print("Not filtering expression values")
			option_arguments$max_components<-2
		}
		if(!("perform_expression_filtering" %in% names(option_arguments))){
			print("Not filtering expression values")
			option_arguments$perform_expression_filtering<-TRUE
		}
		if(!("UMI_bounded_filtering" %in% names(option_arguments))){
			print("Using upper UMI_ounded filtering")
			option_arguments$UMI_bounded_filtering<-"upper"
		}
		if(!("order_by_seurat_varGenes" %in% names(option_arguments))){
			print("Not ordering by seurat_varGenes")
			option_arguments$order_by_seurat_varGenes<-FALSE
		}
		if(!("color_by_seurat_res" %in% names(option_arguments))){
			print("Coloring by seurat resolutions")
			option_arguments$color_by_seurat_res<-TRUE
		}
		if(!("cca_variables" %in% names(option_arguments))){
			print("Not ordering by seurat_varGenes")
			option_arguments$cca_variables<- "~nUMI + nGene"
		}

		return(option_arguments)

	}
}
