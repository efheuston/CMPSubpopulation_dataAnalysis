#' cycle_plot_param
#'
#' @param plotting_function oneof type "cluster", "trajectory", "qplot"
#' @param cycle_parameter the number of colors needed
#' @param the_object monocle object
#'
#' @return
#'
cycle_plot_param<-function(plotting_function, cycle_parameter, the_object){
	# Create color pallete
	basic_color_palette<-basic_color_palette[1:25]
	new_length <-0
	if(length(cycle_parameter) > length(basic_color_palette)){
		new_length<-length(unique(my_object@ident)) - length(basic_color_palette)
	}
	new_length
	basic_color_palette<-c(basic_color_palette, primary.colors(new_length))

	if(plotting_function == "cluster"){
		png_plotFunction(my_plot_cell_clusters(the_object, 1, 2, color = cycle_parameter, point_colors = basic_color_palette) +
										 	ggtitle(paste(output_prefix, "-", cycle_parameter, sep = "")),
										 filename = paste(output_prefix, "_clstr-", cycle_parameter,".png", sep = ""),
										 height = 1600,
										 width = 1600)
		png_plotFunction(my_plot_cell_clusters(the_object, 1, 2, color = cycle_parameter, point_colors = basic_color_palette) +
										 	facet_wrap(as.vector(cycle_parameter)) +
										 	ggtitle(paste(paste(output_prefix, "-", cycle_parameter, sep = ""))),
										 filename = paste(paste(output_prefix, "_clstr-", cycle_parameter,"FACET.png", sep = "")),
										 height = 1600,
										 width = 1600)
		png_plotFunction(plot_cell_clusters(the_object, 1, 2, color = "orig.ident") +
										 	facet_wrap(as.vector(cycle_parameter)) +
										 	ggtitle(paste(output_prefix, "-orig.identFACET", sep = "")),
										 filename = paste(output_prefix, "_clstr-orig.identFACET.png", sep = ""),
										 height = 1600,
										 width = 1600)

	}
	if(plotting_function == "trajectory"){
		png_plotFunction(my_plot_cell_trajectory(the_object,
																						 color_by = cycle_parameter,
																						 cell_size = 1,
																						 point_colors = basic_color_palette) + ggtitle(paste(output_prefix, "-", cycle_parameter, sep = "")),
										 filename = paste(output_prefix, "_trajectory-", cycle_parameter,".png", sep = ""),
										 height = 1600,
										 width = 1600)
		png_plotFunction(my_plot_cell_trajectory(the_object,
																						 color_by = cycle_parameter,
																						 cell_size = 1,
																						 point_colors = basic_color_palette) +
										 	facet_wrap(as.vector(cycle_parameter)) + ggtitle(paste(output_prefix, "-", cycle_parameter, sep = "")),
										 filename = paste(output_prefix, "_trajectory-", cycle_parameter,"FACET.png", sep = ""),
										 height = 1600,
										 width = 1600)
		png_plotFunction(my_plot_cell_trajectory(the_object,
																						 color_by = cycle_parameter,
																						 cell_size = 1,
																						 point_colors = basic_color_palette) +
										 	facet_wrap(~orig.ident) + ggtitle(paste(output_prefix, "-", cycle_parameter,"FACETorig.ident", sep = "")),
										 filename = paste(output_prefix, "_trajectory-", cycle_parameter,"FACETorig.ident.png", sep = ""),
										 height = 1600,
										 width = 1600)
		png_plotFunction(my_plot_cell_trajectory(the_object,
																						 color_by = "orig.ident",
																						 cell_size = 1,
																						 point_colors = basic_color_palette) +
										 	facet_wrap(as.vector(cycle_parameter)) + ggtitle(paste(output_prefix, "-orig.identFACET", cycle_parameter, sep = "")),
										 filename = paste(output_prefix, "_trajectory-orig.identFACET", cycle_parameter,".png", sep = ""),
										 height = 1600,
										 width = 1600)
	}
	if(plotting_function == "qplot"){
		png_plotFunction(qplot(Total_nUMI, data = pData(the_object), color = as.factor(cycle_parameter), geom = "density") +
										 	geom_vline(xintercept = upper_bound) +
										 	geom_vline(xintercept = lower_bound) +
										 	ggtitle(paste(output_prefix, "-", cycle_parameter, sep = "")),
										 filename = paste(output_prefix, "_nUMI_2SD-by", cycle_parameter,".png", sep = ""))

	}
}
