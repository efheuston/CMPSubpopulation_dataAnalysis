#' png_plotFunction
#'
#' @param plot_to_make
#' @param filename
#' @param height
#' @param width
#'
#' @return
#' @export
#'
#' @examples
png_plotFunction<-function(plot_to_make = x, filename = y, height = 800, width = 800){ # default window is 800x800
	png(filename = filename, height = height, width = width)
	plot_to_make
	dev.off()
}
