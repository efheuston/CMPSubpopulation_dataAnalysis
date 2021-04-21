#' png_plotFunction
#'
#' @param plot_to_make plot function
#' @param filename filename to for saved plot
#' @param height png height (px)
#' @param width png width (px)
#'
#' @return
#' @export
#'
png_plotFunction<-function(plot_to_make = x, filename = y, height = 800, width = 800){ # default window is 800x800
	png(filename = filename, height = height, width = width)
	plot_to_make
	dev.off()
}
