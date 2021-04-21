#' runEdgeRAnalysis
#'
#' @return
#' @importFrom tximport tximport
#' @importFrom DESeq2 DESeqDataSetFromMatrix rlog plotPCA
#' @import TxDb.Mmusculus.UCSC.mm10.knownGene
#' @import edgeR
#' @import colorRamps
#' @import org.Mm.eg.db
#' @import ggplot2
#' @import biomaRt
#' @import ggrepel
#' @param file_list Directory containing RSEM files in the format "<file>.genes.results"
#' @param dgelist_groups Character vector of sample names (for factor grouping)
#' @param ref_sample Element from dgelist_groups to define as reference sample
#'
#' @export

runEdgeRAnalysis <- function(file_list, dgelist_groups, ref_sample){

	# Tximport -------------------------------------------------------------------

	# import rsem.genes.results files

	mstxi <- tximport::tximport(files = file_list, type = 'rsem', txIn = FALSE, txOut = FALSE, geneIdCol = 3)


	# SumarizedExperiment object and CPM calculations -------------------------------------------------------------------

	cts <- mstxi$counts
	colnames(cts)

	# Set reference sample =-----------------------------------------------------------


	dgelist_groups <- relevel(dgelist_groups, ref = ref_sample)
	levels(dgelist_groups)

	dgelist <- DGEList(counts = cts, group = dgelist_groups)

	apply(dgelist$counts, 2, sum)

	dim(dgelist)
	keep <- rowSums(cpm(dgelist) > 5) >=2 # only keep genes that have a count of at least 5 in at least 2 samples
	dgelist <- dgelist[keep,]



	# Add gene synonyms -------------------------------------

	mmusculus_genes <- biomaRt::getBM(attributes = c("ensembl_gene_id", "mgi_symbol", "entrezgene_id"),
													 mart = useMart("ensembl", dataset = "mmusculus_gene_ensembl"),
													 useCache = FALSE)


	# create the conversion table
	Mouse2HumanTable <- Mouse2Human(MouseGenes = mmusculus_genes$mgi_symbol)


	# Convert symbol column ot dgelist

	symbol = rownames(dgelist$counts)
	dgelist$genes <- data.frame(Symbol = symbol)

	# Add Hs homologue to dgelist
	HsHmlg <- with(Mouse2HumanTable, Mouse2HumanTable$HGNC[match(row.names(dgelist), Mouse2HumanTable$MGI)])
	dgelist$genes$HsHmlg <- HsHmlg
	dgelist$genes$HsHmlg[dgelist$genes$HsHmlg == ""] <- NA
	dgelist$genes[465:467,]


	# Convert row names to entrez gene ID
	entrezid <- with(Mouse2HumanTable, Mouse2HumanTable$Entrez.Gene_ID[match(row.names(dgelist), Mouse2HumanTable$MGI)])
	entrezid[is.na(entrezid)] <- 'not_mapped'
	rownames(dgelist$counts) <- entrezid
	dgelist <- dgelist[!(row.names(dgelist) %in% 'not_mapped'), ] # get rid of unmapped entrez ids
	head(dgelist)


	# Calculate library size --------------------------------------------------


	dgelist$samples$lib.size <- colSums(dgelist$counts)
	dgelist <- calcNormFactors(dgelist)




	# PCA plot ------------------------------------------------------------

	points <- c(rep(19, 26))
	colors <- c("#080606",
							"#00e6ec",
							"#ef8517",
							"#7352c7",
							"#832c31",
							"#009900",
							"#e2e3a8",
							"#f5001f",
							"#888787",
							"#4b643a",
							"#c747a4",
							"#3908b5",
							"#99e44b")



	pca.design.matrix <- model.matrix(~ 0 + dgelist$samples$group)
	colnames(pca.design.matrix) <- levels(dgelist$samples$group)
	colors <- c('darkviolet', 'cyan2', 'sienna1', 'seagreen3', 'brown4', 'violetred2', 'steelblue2', 'red1', 'darkgreen', 'black', 'royalblue')


	edger.dds <- DESeqDataSetFromMatrix(countData = round(dgelist$counts), colData = dgelist$samples, design = pca.design.matrix)
	transform.edger.dds <- rlog(edger.dds, blind = TRUE)
	pca.dds <- plotPCA(transform.edger.dds, intgroup = c("group"), returnData = TRUE)
	percentVar <- round(100*attr(pca.dds, "percentVar"))

	png(filename = 'RNASeq_pcaPlot.png', height = 1600, width = 1600, bg = 'transparent')
	ggplot(pca.dds, aes(PC1, PC2, color = group)) +
		geom_point(size = 10) +
		xlab(paste0("PC1: ", percentVar[1], "% variance")) +
		ylab(paste0("PC2: ", percentVar[2], "% variance")) +
		coord_fixed() +
		scale_color_manual(values = colors)+
		theme_linedraw()
	dev.off()



	# Plot MDS 3D -------------------------------------------------------------



	scatterplot3d::scatterplot3d(mds.dim,
															 color = colors[dgelist_groups],
															 pch = points[dgelist_groups],
															 lwd = 4,
															 cex.symbols = 8,
															 grid = T, box = F,
															 type = 'h',
															 xlab = 'Coordinate1', ylab = 'Coordinate2', zlab = 'Coordinate3',
															 main = 'RNA_noMatGran_3dPCA',
															 angle = 215)
	legend("topright", legend = levels(dgelist_groups), pch = points, col = colors, ncol = 1, pt.cex = 5, cex = 4)

	coordinatesq<- scatterplot3d::scatterplot3d(mds.dim,
																							color = colors[dgelist_groups],
																							pch = points[dgelist_groups],
																							lwd = 4,
																							cex.symbols = 8,
																							grid = T, box = F,
																							type = 'h',
																							xlab = 'Coordinate1', ylab = 'Coordinate2', zlab = 'Coordinate3',
																							main = 'RNA_noMatGran_3dPCA',
																							angle = 215)


	png(filename = 'RNA_3dPCA-noLegend.png', width = 1600, height = 1600, bg = 'transparent')
	scatterplot3d::scatterplot3d(mds.dim,
															 color = colors[dgelist_groups],
															 pch = points[dgelist_groups],
															 lwd = 4,
															 cex.symbols = 8,
															 grid = T, box = F,
															 type = 'h',
															 xlab = 'Coordinate1', ylab = 'Coordinate2', zlab = 'Coordinate3',
															 main = 'RNA_3dPCA',
															 angle = 215)
	dev.off()



	# Export CPM values for heat Figure 5b -----------------------------------------

	csvcmp.cpm <-cpmByGroup(dgelist)

	# Convert row names to entrez gene ID
	symbolid <- with(Mouse2HumanTable, Mouse2HumanTable$MGI[match(row.names(csvcmp.cpm), Mouse2HumanTable$Entrez.Gene_ID)])
	rownames(csvcmp.cpm) <- symbolid
	write.table(csvcmp.cpm, file = "avgCPM_forCXvCMP.txt", quote = FALSE, sep = '\t')



	# Estimate dispersion (GLM) - WITH CMP Intercept------------------------------------------------------------

	designMat <- model.matrix(~dgelist$samples$group)
	colnames(designMat) <- c('(Intercept)', levels(dgelist$samples$group)[2:length(levels(dgelist$samples$group))])
	designMat

	# Apply GLM dispersion ------------------------------------------------------------

	dgelist <- estimateGLMCommonDisp(dgelist, designMat)
	dgelist <- estimateGLMTrendedDisp(dgelist, designMat, method = 'bin.spline')
	dgelist <- estimateGLMTagwiseDisp(dgelist, designMat)

	# GLM differential expression testing ------------------------------------------------------------
	fit <- glmFit(dgelist, designMat)

	lrt.c10 <-  glmLRT(fit, coef = 2)
	de.c10 <- decideTestsDGE(lrt.c10, adjust.method = 'BH', p.value = 0.05)
	detags.c10 <- rownames(dgelist)[as.logical(de.c10)]
	plotSmear(lrt.c10, de.tags = detags.c10)
	abline(h =  c(-2, 2), col = 'blue')
	goana.table <- goana(lrt.c10[row.names(lrt.c10) %in% detags.c10,], species = "Mm")
	topgo.table <- topGO(goana.table, ontology = c("BP", "MF"), number = Inf)
	write.table(topgo.table, file ='C10vsIntercept_edgeR-GOAll.txt', quote = F, row.names = T, sep = '\t')


	lrt.c11 <-  glmLRT(fit, coef = 3)
	de.c11 <- decideTestsDGE(lrt.c11, adjust.method = 'BH', p.value = 0.05)
	detags.c11 <- rownames(dgelist)[as.logical(de.c11)]
	plotSmear(lrt.c11, de.tags = detags.c11)
	abline(h =  c(-2, 2), col = 'blue')
	goana.table <- goana(lrt.c11[row.names(lrt.c11) %in% detags.c11,], species = "Mm")
	topgo.table <- topGO(goana.table, ontology = c("BP", "MF"), number = Inf)
	write.table(topgo.table, file ='C11vsIntercept_edgeR-GOAll.txt', quote = F, row.names = T, sep = '\t')


	lrt.c17 <-  glmLRT(fit, coef = 4)
	de.c17 <- decideTestsDGE(lrt.c17, adjust.method = 'BH', p.value = 0.05)
	detags.c17 <- rownames(dgelist)[as.logical(de.c17)]
	plotSmear(lrt.c17, de.tags = detags.c17)
	abline(h =  c(-2, 2), col = 'blue')
	goana.table <- goana(lrt.c17[row.names(lrt.c17) %in% detags.c17,], species = "Mm")
	topgo.table <- topGO(goana.table, ontology = c("BP", "MF"), number = Inf)
	write.table(topgo.table, file ='C17vsIntercept_edgeR-GOAll.txt', quote = F, row.names = T, sep = '\t')


	lrt.c3 <-  glmLRT(fit, coef = 5)
	de.c3 <- decideTestsDGE(lrt.c3, adjust.method = 'BH', p.value = 0.05,)
	detags.c3 <- rownames(dgelist)[as.logical(de.c3)]
	plotSmear(lrt.c3, de.tags = detags.c3)
	abline(h =  c(-2, 2), col = 'blue')
	goana.table <- goana(lrt.c3[row.names(lrt.c3) %in% detags.c3,] , species = "Mm")
	topgo.table <- topGO(goana.table, ontology = c("BP", "MF"), number = Inf)
	write.table(topgo.table, file ='c3vsIntercept_edgeR-GOAll.txt', quote = F, row.names = T, sep = '\t')



	# Get genes uniquely upregulated in each cluster ------------------------------------------------------------

	c10unique <- setdiff(detags.c10, detags.c11)
	c10unique <- setdiff(c10unique, detags.c17)
	c10unique <- setdiff(c10unique, detags.c3)
	dgelist.c10unique <- dgelist[(row.names(dgelist) %in% c10unique), ]
	lrt.c10.unique <- lrt.c10[(row.names(lrt.c10) %in% c10unique),]
	goana.table <- goana(lrt.c10.unique, species = "Mm")
	topgo.table <- topGO(goana.table, ontology = c("BP", "MF"), number = Inf)
	write.table(topgo.table, file ='C10vsIntercept_edgeR_unq-GOAll.txt', quote = F, row.names = T, sep = '\t')

	c11unique <- setdiff(detags.c11, detags.c10)
	c11unique <- setdiff(c11unique, detags.c17)
	c11unique <- setdiff(c11unique, detags.c3)
	dgelist.c11unique <- dgelist[(row.names(dgelist) %in% c11unique), ]
	lrt.c11.unique <- lrt.c11[(row.names(lrt.c11) %in% c11unique),]
	goana.table <- goana(lrt.c11.unique, species = "Mm")
	topgo.table <- topGO(goana.table, ontology = c("BP", "MF"), number = Inf)
	write.table(topgo.table, file ='C11vsIntercept_edgeR_unq-GOAll.txt', quote = F, row.names = T, sep = '\t')

	c17unique <- setdiff(detags.c17, detags.c11)
	c17unique <- setdiff(c17unique, detags.c10)
	c17unique <- setdiff(c17unique, detags.c3)
	dgelist.c17unique <- dgelist[(row.names(dgelist) %in% c17unique), ]
	lrt.c17.unique <- lrt.c17[(row.names(lrt.c17) %in% c17unique),]
	goana.table <- goana(lrt.c17.unique, species = "Mm")
	topgo.table <- topGO(goana.table, ontology = c("BP", "MF"), number = Inf)
	write.table(topgo.table, file ='C17vsIntercept_edgeR_unq-GOAll.txt', quote = F, row.names = T, sep = '\t')

	c3unique <- setdiff(detags.c3, detags.c11)
	c3unique <- setdiff(c3unique, detags.c17)
	c3unique <- setdiff(c3unique, detags.c10)
	dgelist.c3unique <- dgelist[(row.names(dgelist) %in% c3unique), ]
	lrt.c3.unique <- lrt.c3[(row.names(lrt.c3) %in% c3unique),]
	goana.table <- goana(lrt.c3.unique, species = "Mm")
	topgo.table <- topGO(goana.table, ontology = c("BP", "MF"), number = Inf)
	write.table(topgo.table, file ='c3vsIntercept_edgeR_unq-GOAll.txt', quote = F, row.names = T, sep = '\t')



	# GLM differential expression testing ------------------------------------------------------------
	fit <- glmFit(dgelist, designMat)

	FC = 1.5

	lrt.c10 <-  glmLRT(fit, coef = 2)
	de.c10 <- decideTestsDGE(lrt.c10, adjust.method = 'BH', p.value = 0.05, lfc = log2(FC))
	detags.c10 <- rownames(dgelist)[as.logical(de.c10)]
	plotSmear(lrt.c10, de.tags = detags.c10)
	abline(h =  c(-2, 2), col = 'blue')
	makeSmearPlot(lrt = lrt.c10, nameA = 'c10', nameB = 'CMP', fc = FC)
	goana.table <- goana(lrt.c10[(row.names(lrt.c10) %in% detags.c10),], species = "Mm")
	topgo.table <- topGO(goana.table, ontology = c("BP", "MF"), number = Inf)
	write.table(topgo.table, file =paste('C10vsIntercept_edgeR-FC', FC, '-GOAll.txt', sep = ''), quote = F, row.names = T, sep = '\t')


	lrt.c11 <-  glmLRT(fit, coef = 3)
	de.c11 <- decideTestsDGE(lrt.c11, adjust.method = 'BH', p.value = 0.05, lfc = log2(FC))
	detags.c11 <- rownames(dgelist)[as.logical(de.c11)]
	plotSmear(lrt.c11, de.tags = detags.c11)
	abline(h =  c(-2, 2), col = 'blue')
	makeSmearPlot(lrt = lrt.c11, nameA = 'c11', nameB = 'CMP', fc = FC)
	goana.table <- goana(lrt.c11[(row.names(lrt.c11) %in% detags.c11),], species = "Mm")
	topgo.table <- topGO(goana.table, ontology = c("BP", "MF"), number = Inf)
	write.table(topgo.table, file =paste('C11vsIntercept_edgeR-FC', FC, '-GOAll.txt', sep = ''), quote = F, row.names = T, sep = '\t')


	lrt.c17 <-  glmLRT(fit, coef = 4)
	de.c17 <- decideTestsDGE(lrt.c17, adjust.method = 'BH', p.value = 0.05, lfc = log2(FC))
	detags.c17 <- rownames(dgelist)[as.logical(de.c17)]
	plotSmear(lrt.c17, de.tags = detags.c17)
	abline(h =  c(-2, 2), col = 'blue')
	makeSmearPlot(lrt = lrt.c17, nameA = 'c17', nameB = 'CMP', fc = FC)
	goana.table <- goana(lrt.c17[(row.names(lrt.c17) %in% detags.c17),], species = "Mm")
	topgo.table <- topGO(goana.table, ontology = c("BP", "MF"), number = Inf)
	write.table(topgo.table, file =paste('C17vsIntercept_edgeR-FC', FC, '-GOAll.txt', sep = ''), quote = F, row.names = T, sep = '\t')


	lrt.c3 <-  glmLRT(fit, coef = 5)
	de.c3 <- decideTestsDGE(lrt.c3, adjust.method = 'BH', p.value = 0.05, lfc = log2(FC))
	detags.c3 <- rownames(dgelist)[as.logical(de.c3)]
	plotSmear(lrt.c3, de.tags = detags.c3)
	abline(h =  c(-2, 2), col = 'blue')
	makeSmearPlot(lrt = lrt.c3, nameA = 'c3', nameB = 'CMP', fc = FC)
	goana.table <- goana(lrt.c3[(row.names(lrt.c3) %in% detags.c3),], species = "Mm")
	topgo.table <- topGO(goana.table, ontology = c("BP", "MF"), number = Inf)
	write.table(topgo.table, file =paste('c3vsIntercept_edgeR-FC', FC, '-GOAll.txt', sep = ''), quote = F, row.names = T, sep = '\t')





	# Write counts files ----------------------------------------------
	counts.frame <- as.matrix(dgelist$counts)


	for(i in 1:ncol(counts.frame)){
		write.table(x = counts.frame[,i, drop = FALSE], file = paste0(colnames(counts.frame)[i], '_edgeRcounts.txt'), row.names = TRUE, quote = FALSE, sep = '\t', col.names = FALSE)
	}

	# Get genes uniquely upregulated in each cluster ------------------------------------------------------------

	c10unique <- setdiff(detags.c10, detags.c11)
	c10unique <- setdiff(c10unique, detags.c17)
	c10unique <- setdiff(c10unique, detags.c3)
	dgelist.c10unique <- dgelist[(row.names(dgelist) %in% c10unique), ]
	lrt.c10.unique <- lrt.c10[(row.names(lrt.c10) %in% c10unique),]
	goana.table <- goana(lrt.c10.unique, species = "Mm")
	topgo.table <- topGO(goana.table, ontology = c("BP", "MF"), number = Inf)
	write.table(topgo.table, file ='C10vsIntercept_edgeR_unq-FC15-GOAll.txt', quote = F, row.names = T, sep = '\t')

	c11unique <- setdiff(detags.c11, detags.c10)
	c11unique <- setdiff(c11unique, detags.c17)
	c11unique <- setdiff(c11unique, detags.c3)
	dgelist.c11unique <- dgelist[(row.names(dgelist) %in% c11unique), ]
	lrt.c11.unique <- lrt.c11[(row.names(lrt.c11) %in% c11unique),]
	goana.table <- goana(lrt.c11.unique, species = "Mm")
	topgo.table <- topGO(goana.table, ontology = c("BP", "MF"), number = Inf)
	write.table(topgo.table, file ='C11vsIntercept_edgeR_unq-FC15-GOAll.txt', quote = F, row.names = T, sep = '\t')

	c17unique <- setdiff(detags.c17, detags.c11)
	c17unique <- setdiff(c17unique, detags.c10)
	c17unique <- setdiff(c17unique, detags.c3)
	dgelist.c17unique <- dgelist[(row.names(dgelist) %in% c17unique), ]
	lrt.c17.unique <- lrt.c17[(row.names(lrt.c17) %in% c17unique),]
	goana.table <- goana(lrt.c17.unique, species = "Mm")
	topgo.table <- topGO(goana.table, ontology = c("BP", "MF"), number = Inf)
	write.table(topgo.table, file ='C17vsIntercept_edgeR_unq-FC15-GOAll.txt', quote = F, row.names = T, sep = '\t')

	c3unique <- setdiff(detags.c3, detags.c11)
	c3unique <- setdiff(c3unique, detags.c17)
	c3unique <- setdiff(c3unique, detags.c10)
	dgelist.c3unique <- dgelist[(row.names(dgelist) %in% c3unique), ]
	lrt.c3.unique <- lrt.c3[(row.names(lrt.c3) %in% c3unique),]
	goana.table <- goana(lrt.c3.unique, species = "Mm")
	topgo.table <- topGO(goana.table, ontology = c("BP", "MF"), number = Inf)
	write.table(topgo.table, file ='c3vsIntercept_edgeR_unq-FC15-GOAll.txt', quote = F, row.names = T, sep = '\t')





	# Automate pairwise comparisons -------------------------------------------
	cluster.contrasts <- makeContrasts(c10vsc11 = C10-C11,
																		 c10vsc17 = C10-C17,
																		 c10vsc3 = C10-c3,

																		 c11vsc10 = C11-C10,
																		 c11vsc17 = C11-C17,
																		 c11vsc3 = C11-c3,

																		 c17vsc10 = C17-C10,
																		 c17vsc11 = C17-C11,
																		 c17vsc3 = C17-c3,

																		 c3vsc10 = c3-C10,
																		 c3vsc11 = c3-C11,
																		 c3vsc17 = c3-C17,
																		 levels = designMat
	)

	for(i in 1:(ncol(cluster.contrasts))){
		outname <- names(cluster.contrasts[1,])[i]
		print(outname)
		lrt <- glmLRT(fit, contrast = cluster.contrasts[,i])
		# Create rankings
		all.topTags <- topTags(lrt, n = nrow(lrt$table))$table
		all.topTags <- all.topTags[all.topTags$FDR <= 0.05,]
		all.topTags$ranking <- all.topTags$logFC * -log10(all.topTags$FDR)
		all.topTags<-all.topTags[with(all.topTags, order(-abs(ranking))),]

		# Define ID set to plot based on ranking
		ids <- c()
		for (x in 1:nrow(all.topTags)){
			if(x <=50){
				ids <- c(ids, rownames(all.topTags)[x])
			}
		}

		gene.labels <- lrt$table[row.names(lrt$table) %in% ids,]
		gene.labels$Symbol <- lrt$genes$Symbol[row.names(lrt$table) %in% row.names(gene.labels)]

		# de.summary <- summary(decideTests(lrt,adjust.method = 'BH', p.value = 0.05))
		de.summary <- summary(decideTests(lrt,adjust.method = 'BH', p.value = 0.05, lfc = log2(1.5)))

		# Create plot
		my.plot <- ggplot(data = lrt$table, mapping = aes(x = logCPM,  y = logFC)) +
			geom_point(alpha = 0.3) +
			geom_point(data = gene.labels, aes(x = logCPM, y = logFC), color = 'red') +
			geom_text_repel(data = gene.labels, mapping = aes(x = logCPM, y = logFC, label = Symbol), size = 7) +
			scale_x_continuous(limits = c(-0.3, max(lrt$table$logCPM) + 1), ) +
			scale_y_continuous(limits = c(min(lrt$table$logFC) - 1, max(lrt$table$logFC) + 1), labels = abs) +
			annotate('rect', xmin = -0.05, xmax = 0.05, ymin = 0, ymax = max(lrt$table$logFC) + 1, color = 'red', fill = 'red') +
			annotate('rect', xmin = -0.05, xmax = 0.05, ymin = 0, ymax = min(lrt$table$logFC) - 1, color = 'blue', fill = 'blue') +
			annotate(geom = 'text', x = -0.3, y = max(lrt$table$logFC)/2, label = paste('Higher in', row.names(cluster.contrasts)[cluster.contrasts[,i] == 1]), angle = 90, cex = 7) +
			annotate(geom = 'text', x = -0.3, y = min(lrt$table$logFC)/2, label = paste('Higher in', row.names(cluster.contrasts)[cluster.contrasts[,i] == -1]), angle = 90, cex = 7) +
			labs(title = outname, subtitle = paste('Down: ', de.summary[1,1], '; NotSig: ', de.summary[2,1], '; Up: ', de.summary[3,1], sep = '')) +
			theme(axis.title = element_text(size = 24, face = 'bold'), axis.text = element_text(size = 18, face = 'bold')) +
			theme(plot.title = element_text(face = 'bold', size = 24), plot.subtitle = element_text(face = 'bold', size = 20))



		# Save plot
		png(filename = paste(outname, '-Intercept.png', sep = ''), height = 1600, width = 1600)
		plot(my.plot)
		dev.off()


		# goana/topGo
		goana.table <- goana(lrt, species = "Mm")
		topgo.table <- topGO(goana.table, ontology = c("BP", "MF"), number = Inf)
		write.table(topgo.table, file =paste(outname, '_edgeR-GOAll.txt', sep = ''), quote = F, row.names = T, sep = '\t')

		# Export Symbol/ranking
		all.topTags<-all.topTags[with(all.topTags, order(-ranking)),]
		write.table(all.topTags[, c('HsHmlg', 'ranking')],
								file = paste(outname, 'Intercept-rankedSymbol.txt', sep = ""),
								sep = '\t',
								row.names = FALSE,
								col.names = FALSE,
								quote = FALSE)

		# Export for GSEA
		all.topTags<-all.topTags[with(all.topTags, order(-ranking)),]
		write.table(all.topTags[, c('HsHmlg', 'ranking')],
								file = paste(outname, 'Intercept-rankedGSEA.rnk', sep = ""),
								sep = '\t',
								row.names = FALSE,
								col.names = TRUE,
								quote = FALSE,
								na = "not_mapped")

		# all.topTags<-all.topTags[with(all.topTags, order(ranking)),]
		# write.table(all.topTags[all.topTags$ranking <=-1.2, c('Symbol', 'ranking')],
		# 						file = paste(outname, 'noIntercept-rankedGSEA_neg.txt', sep = ""),
		# 						sep = '\t',
		# 						row.names = FALSE,
		# 						col.names = FALSE,
		# 						quote = FALSE)

	}
}


