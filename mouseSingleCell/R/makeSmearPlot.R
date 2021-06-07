#' makeSmearPlot
#'
#' @param lrt dge object
#' @param nameA comparison string to describe A
#' @param nameB comparison string to describe  B
#' @param outname name to use for saved image
#' @param fc fold-change threshold
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggrepel geom_label_repel
#' @return smear plot
#' @export
#'
makeSmearPlot <- function(lrt, nameA = 'nameA', nameB = 'CMP', fc = FC, outname = outname){

	outname <- paste(nameA, 'vs', nameB, sep = '')

	# Code to overlay IDs on SmearPlot for IDs within top n largest logFC of the topTags list
	ids <- row.names(topTags(lrt, n = 50, sort.by = 'logFC')$table)
	gene.labels <- lrt$table[row.names(lrt$table) %in% ids,]

	gene.labels <- lrt$table[row.names(lrt$table) %in% ids,]
	gene.labels$Symbol <- lrt$genes$Symbol[row.names(lrt$table) %in% row.names(gene.labels)]

	# de.summary <- summary(decideTests(lrt,adjust.method = 'BH', p.value = 0.05))
	de.summary <- summary(decideTests(lrt,adjust.method = 'BH', p.value = 0.05, lfc = log2(fc)))

	# Create plot
	my.plot <- ggplot2::ggplot(data = lrt$table, mapping = aes(x = logCPM,  y = logFC)) +
		geom_point(alpha = 0.3) +
		geom_point(data = gene.labels, aes(x = logCPM, y = logFC), color = 'red') +
		geom_text_repel(data = gene.labels, mapping = aes(x = logCPM, y = logFC, label = Symbol), size = 7) +
		scale_x_continuous(limits = c(-0.3, max(lrt$table$logCPM) + 1), ) +
		scale_y_continuous(limits = c(min(lrt$table$logFC) - 1, max(lrt$table$logFC) + 1), labels = abs) +
		annotate('rect', xmin = -0.05, xmax = 0.05, ymin = 0, ymax = max(lrt$table$logFC) + 1, color = 'red', fill = 'red') +
		annotate('rect', xmin = -0.05, xmax = 0.05, ymin = 0, ymax = min(lrt$table$logFC) - 1, color = 'blue', fill = 'blue') +
		annotate(geom = 'text', x = -0.3, y = max(lrt$table$logFC)/2, label = paste('Higher in', nameA), angle = 90, cex = 7) +
		annotate(geom = 'text', x = -0.3, y = min(lrt$table$logFC)/2, label = paste('Higher in', nameB), angle = 90, cex = 7) +
		# text(x=gene.labels$logCPM, y=gene.labels$logFC, labels=lrt$genes$Symbol[rownames(lrt$table) %in% row.names(gene.labels)], cex=.7, pos = 4, col = 'blue')
		labs(title = outname, subtitle = paste('Down: ', de.summary[1,1], '; NotSig: ', de.summary[2,1], '; Up: ', de.summary[3,1], sep = '')) +
		theme(axis.title = element_text(size = 24, face = 'bold'), axis.text = element_text(size = 18, face = 'bold')) +
		theme(plot.title = element_text(face = 'bold', size = 24), plot.subtitle = element_text(face = 'bold', size = 20))



	# Save plot
	png(filename = outname, height = 1600, width = 1600)
	plot(my.plot)
	dev.off()

}

