#' read_seurat_args
#'
#' @return
#'
#'
read_seurat_args<-function(data, projectName, genome, max_pcs, resolutionList, vars_to_regress){
	args<-paste(unlist(args), collapse = " ")
	args<-unlist(strsplit(args, "--"))[-1]
	option_arguments<-sapply(args, function(x){
		unlist(strsplit(x, " "))[-1]
	})
	option_names<-sapply(args, function(x){
		unlist(strsplit(x, " "))[1]
	})
	names(option_arguments) <- unlist(option_names)
	required_keys<-c("data", "projectName", "genome", "max_pcs", "resolutionList")

	if(!(all(required_keys %in% names(option_arguments)))){
		print(paste("Missing:", paste(required_keys[!(required_keys%in%names(option_arguments))], collapse = " "), sep = " "))
		stop("Required arguments:
           --data: folders containing outs files (space seprated)
           --projectName: name to append to printed files
           --genome: genome (mm10 or GRCh38)
           --max_pcs: dimensionality of reduced space
           --resolutionList:
          Additional arguments:
           --mitoFilter: whether or not to correct for cell cycle stage (default = FALSE)
           --vars_to_regress: variables to regress on (default = c(\"nUMI\", \"nGene\"))
           --perform_cca: whether or not to perform multi canonical clustering alignment
         Sample Input:
           R --vanilla < SeuratPipeline.R --args --data DBA526B DBA526C DBA677 --projectName 20181207_pbDBA --genome GRCh38 --max_pcs 20 --resolutionList 1.0 --perform_cca TRUE")
	} else{
		if(!("mitoFilter" %in% names(option_arguments))){
			print("Not filtering based on mt-gene expression")
			option_arguments$mitoFilter<-FALSE
		}else{
			print("Filtering based on mt-gene expression")
		}
		if(!("vars_to_regress" %in% names(option_arguments))){
			option_arguments$vars_to_regress<-c("nGene", "nUMI")
		}
		if(!("perform_cca" %in% names(option_arguments))){
			option_arguments$perform_cca<-FALSE
		}
		return(option_arguments)
	}
}
