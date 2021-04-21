#' adjust_palette_size
#'
#' @param object_length Num of parameters to color
#'
#' @return
#' @export
#'
adjust_palette_size <- function(object_length, basic_color_palette = c("#cb4bbe",
																																			 "lightskyblue",
																																			 "grey37",
																																			 "#53ba48",
																																			 "moccasin",
																																			 "#dcca44",
																																			 "#502d70",
																																			 "#afd671",
																																			 "#cb4471",
																																			 "#69da99",
																																			 "#d44732",
																																			 "#6fd3cf",
																																			 "#5d2539",
																																			 "#cfd0a0",
																																			 "blue",
																																			 "#d2883b",
																																			 "#6793c0",
																																			 "#898938",
																																			 "#c98ac2",
																																			 "yellow",
																																			 "#c4c0cc",
																																			 "#7d3d23",
																																			 "#00a5ff",
																																			 "#d68d7a",
																																			 "#a2c1a3")){
	if(length(unique(object_length)) > length(basic_color_palette)){
		new_length <- length(unique(object_length)) - length(basic_color_palette)
		my_palette <- c(basic_color_palette, primary.colors(new_length))
	} else {
		my_palette <- basic_color_palette
	}
}

