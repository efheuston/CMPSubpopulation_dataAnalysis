#' runSeuratAnalysis
#'
#' @param args list of arguments (data, project name, genome, max_pcs, resolutionList, perform_cca, vars_to_regress)
#'
#' @return
#' @export
#' @import Seurat
#' @import dplyr
#' @import colorRamps
#' @import parallel
#' @import future
#' @import R.utils
#'
#' \dontrun{
#' args <- c('--data CMP-m1 CMP-m2 --projectName CMP12k_vs_CMP3k --genome mm10 --max_pcs 25 --resolutionList 0.6 --perform_cca FALSE --vars_to_regress nUMI nGene orig.ident')
#' runSeuratAnalysis(args)
#'	}

runSeuratAnalysis<- function(args){

	option_arguments<-read_seurat_args()

	print(paste("Performing clustering on", paste(option_arguments$data, collapse = " "), sep = " "))
	print(paste("Creating files with the output name", option_arguments$projectName, sep = " "))
	print(paste("Running seurat with genome", option_arguments$genome, "on", option_arguments$max_pcs, "principle components at resolutions", option_arguments$resolutionList, "with regression vars", paste(option_arguments$vars_to_regress, collapse = " "), sep = " "))



	sink(paste0(option_arguments$projectName, '_sessionInfo.txt'))
	sessionInfo()
	sink()

	# Set cca-dependent variables ---------------------------------------------------
	if(option_arguments$perform_cca == TRUE){
		option_arguments$projectName<-paste(option_arguments$projectName, "_cca", sep = "")
		reduction_type <- "cca.aligned"
		tsne_reduction_input <- "cca.aligned"
	} else{
		reduction_type <-"pca"
	}

	# Set mt_filtering-dependent variables
	if(option_arguments$mitoFilter == TRUE){
		option_arguments$data<-paste(option_arguments$projectName, "_mt", sep = "")
		option_arguments$vars_to_regress <- c(option_arguments$vars_to_regress, "percent.mito")
	}


	# define PCA plotting function that colors one variable red and all others black, and cycles through ##can change colors, too

	# define function to change case of cell cycle genes to match mouse gene nomenclature
	if(option_arguments$genome == "mm10"){
		s.genes<-c("Mcm4", "Exo1", "Slbp", "Gmnn", "Cdc45", "Msh2", "Mcm6", "Rrm2", "Pold3", "Blm", "Ubr7", "Mcm5", "Clspn", "Hells", "Nasp", "Rpa2", "Rad51ap1", "Tyms", "Rrm1", "Rfc2", "Prim1", "Brip1", "Usp1", "Ung", "Pola1", "Mcm2", "Fen1", "Tipin", "Pcna", "Cdca7", "Uhrf1", "Casp8ap2", "Cdc6", "Dscc1", "Wdr76", "E2f8", "Dtl", "Ccne2", "Atad2", "Gins2", "Chaf1b", "Pcna-ps2")
		g2m.genes<-c("Nuf2", "Psrc1", "Ncapd2", "Ccnb2", "Smc4", "Lbr", "Tacc3", "Cenpa", "Kif23", "Cdca2", "Anp32e", "G2e3", "Cdca3", "Anln", "Cenpe", "Gas2l3", "Tubb4b", "Cenpf", "Dlgap5", "Hjurp", "Cks1brt", "Gtse1", "Bub1", "Birc5", "Ube2c", "Rangap1", "Hmmr", "Ect2", "Tpx2", "Ckap5", "Cbx5", "Nek2", "Ttk", "Cdca8", "Nusap1", "Ctcf", "Cdc20", "Cks2", "Mki67", "Tmpo", "Ckap2l", "Aurkb", "Kif2c", "Cdk1", "Kif20b", "Top2a", "Aurka", "Ckap2", "Hmgb2", "Cdc25c", "Ndc80", "Kif11")  # list of mouse cell cycle genes provided by leonfodoulian from https://github.com/satijalab/seurat/issues/462; ##is a converted set of the seurat cc.genes list
	} else if(option_arguments$genome == "GRCh38") {
		s.genes<-Seurat::cc.genes$s.genes
		g2m.genes<-Seurat::cc.genes$g2m.genes
	} else{
		print ("No cell cycle genes match this your option_arguments$genome")
	}
	# Load data ---------------------------------------------------------------

	if(option_arguments$perform_cca == TRUE){
		multi_object_list<-list()
		multi_object_list<-lapply(1:length(option_arguments$data), function(x) create_multi_object_list(option_arguments$data[x]))
		names(multi_object_list)<-option_arguments$data
		multi_object_list

		genes.use<-c()
		for (i in 1:length(multi_object_list)){
			if(length(multi_object_list[[i]]@var.genes)>=1000){
				genes.use<-c(genes.use, head(rownames(multi_object_list[[i]]@hvg.info), 1000))
				print(length(paste("Using", genes.use, "genes to align subspace", sep = " ")))
			} else {
				genes.use<-c(genes.use, head(rownames(multi_object_list[[i]]@hvg.info)))
				print(length(paste("Using", genes.use, "genes to align subspace", sep = " ")))
			}
		}
		genes.use<-names(which(table(genes.use)>1))
		for (i in 1:length(multi_object_list)){
			genes.use<-genes.use[genes.use %in% rownames(multi_object_list[[i]]@scale.data)]
		}

		if(length(multi_object_list) >=3){
			my_object<-RunMultiCCA(multi_object_list, genes.use = genes.use, num.ccs = as.numeric(option_arguments$max_pcs))
		} else if(length(multi_object_list) == 2){
			my_object<- RunCCA(object = multi_object_list[[1]], object2 = multi_object_list[[2]], num.cc = as.numeric(option_arguments$max_pcs))
		}
		p1<-DimPlot(my_object, reduction.use="cca", group.by="orig.ident", do.return=TRUE)
		p2<-VlnPlot(my_object, features.plot="CC1", group.by = "orig.ident", do.return = TRUE)

		png(filename = paste(option_arguments$projectName,"_dim",option_arguments$max_pcs, "_CCFit.png", sep = ""), height = 800, width = 800)
		plot_grid(p1, p2)
		dev.off()

		png(filename = paste(option_arguments$projectName,"_dim",option_arguments$max_pcs, "_BicorPlot.png", sep = ""), height = 800, width = 800)
		MetageneBicorPlot(my_object, grouping.var = "orig.ident", dims.eval = 1:as.numeric(option_arguments$max_pcs))
		dev.off()

		png(filename = paste(option_arguments$projectName,"_dim",option_arguments$max_pcs, "_Heatmap.png", sep = ""), height = 2400, width = 800)
		DimHeatmap(my_object, reduction.type = "cca", cells.use = 500, dim.use = 1:as.numeric(option_arguments$max_pcs), do.balanced = TRUE)
		dev.off()


		my_object<-CalcVarExpRatio(my_object, reduction.type = "pca", grouping.var = "orig.ident", dims.use = 1:as.numeric(option_arguments$max_pcs))
		my_object<-SubsetData(my_object, subset.name = "var.ratio.pca", accept.low = 0.5)
		my_object<-AlignSubspace(my_object, reduction.type = "cca", grouping.var = "orig.ident", dims.align = 1:as.numeric(option_arguments$max_pcs))

	} else {

		data_files.list<-sapply(option_arguments$data[1:length(option_arguments$data)], function(x) read_10XGenomics_data(x))
		my_object.data<-Read10X(data_files.list)
		my_object<- CreateSeuratObject(raw.data = my_object.data, min.cells = 3, min.genes = 200, project = option_arguments$projectName)
		my_object<- create_percentMito_column(my_object)

		par(mfrow=c(1,1))
		png(paste(option_arguments$projectName, "_VlnStats.png", sep = ""), width = 800, height = 800)
		VlnPlot(my_object, features.plot = c('nGene', 'nUMI'), nCol = 2)
		dev.off()

		if(option_arguments$mitoFilter == TRUE){
			png(paste(option_arguments$projectName, "_GenePlots.png", sep = ""), width = 800, height = 800)
			par(mfrow = c(1,2))
			GenePlot(my_object, gene1 = "nUMI", gene2 = "percent.mito")
			GenePlot(my_object, gene1 = "nUMI", gene2 = "nGene")
			dev.off()
			par(mfrow=c(1,1))
		} else {
			png(paste(option_arguments$projectName, "_GenePlots.png", sep = ""), width = 800, height = 800)
			GenePlot(my_object, gene1 = "nUMI", gene2 = "nGene")
			dev.off()
		}

		# Normalize and scale data ----------------------------------------------------------

		my_object<-NormalizeData(my_object, normalization.method = "LogNormalize", scale.factor = 10000)


		#Variable gene detection helps control for relationship between variability and average expression
		#Set to mark visual outliers, potentially ~2000

		my_object<-FindVariableGenes(my_object,mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
		# VariableGenePlot(my_object)

		png(paste(option_arguments$projectName, "_VarGenePlot.png", sep = ""), width = 800, height = 800)
		VariableGenePlot(my_object)
		dev.off()


		my_object<-ScaleData(my_object, vars.to.regress = option_arguments$vars_to_regress)
		my_object<-CellCycleScoring(my_object, s.genes = s.genes, g2m.genes = g2m.genes, set.ident = FALSE)

		# Run PCA -----------------------------------------------------------------
		my_object<-RunPCA(my_object, pc.genes = my_object@var.genes, pcs.compute = as.numeric(option_arguments$max_pcs), do.print = TRUE, pcs.print = 1:5, genes.print = 5)

		plot_title<-paste(option_arguments$projectName, "dim", option_arguments$max_pcs, sep = "")

		PrintPCA(my_object, pcs.print = 1:5, genes.print = 5, use.full = FALSE) #Set use.full to TRUE to see projected PCA
		VizPCA(my_object, pcs.use = 1:5, use.full = FALSE, font.size = 1)
		PCAPlot(my_object, dim.1 = 1, dim.2 = 2, use.full = FALSE)
		PCHeatmap(my_object, pc.use = 1:10, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE, do.return=TRUE)

		my_plot<-PCAPlot(my_object, dim.1 = 1, dim.2 = 2, use.full = FALSE, do.return=TRUE)
		my_plot<-my_plot+ggtitle(plot_title)
		png(paste(option_arguments$projectName,"dim",option_arguments$max_pcs, "_PCAPlot.png", sep = ""), width = 800, height = 800)
		plot(my_plot)
		dev.off()

		VizPCA(my_object, pcs.use = 1:5, use.full = FALSE, font.size = 1.5)
		png(paste(option_arguments$projectName, "dim",option_arguments$max_pcs, "_VizPCA.png", sep = ""), width = 800, height = 1400)
		par(mar=c(5,10,5,5))
		VizPCA(my_object, pcs.use = 1:5, use.full = FALSE, font.size = 1.5)
		dev.off()
		PCHeatmap(my_object, pc.use = 1:10, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE, do.return=FALSE)
		png(paste(option_arguments$projectName,"dim",option_arguments$max_pcs, "_PCAHeatmap.png", sep = ""), width = 800, height = 1400)
		par(mar=c(5,5,5,10))
		PCHeatmap(my_object, pc.use = 1:10, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE, do.return = FALSE)
		dev.off()
		par(mar=c(0,0,0,0))


		# Find statistically signficant principal components -----------------------

		# my_object<-JackStraw(my_object, num.replicate = 100, num.pc = as.numeric(option_arguments$max_pcs), display.progress = TRUE, do.par = TRUE, num.cores = future::availableCores())

		PCElbowPlot(my_object)
		# JackStrawPlot(my_object, PCs = 1:as.numeric(option_arguments$max_pcs))


		png(paste(option_arguments$projectName,"dim",option_arguments$max_pcs, "_ElbowPlot.png", sep = ""), width = 800, height = 800)
		PCElbowPlot(my_object)
		dev.off()
		# png(paste(option_arguments$projectName,"dim",option_arguments$max_pcs, "_JackStraw.png", sep = ""), width = 800, height = 800)
		# JackStrawPlot(my_object, PCs = 1:as.numeric(option_arguments$max_pcs))
		# dev.off()
		#
	}

	# Write object for monocle-import

	# try(saveRDS(my_object, file = paste(option_arguments$projectName, "_monocleObject.rds", sep = "")), silent = FALSE)


	# create folders for feature plots

	suppressWarnings(dir.create("ViolinPlots"))
	suppressWarnings(dir.create("FeaturePlots"))

	# Clustering --------------------------------------------------------------

	for(resolution in option_arguments$resolutionList){
		resolution<-as.numeric(resolution)
		my_object<-FindClusters(my_object, reduction.type = reduction_type, dims.use = 1:as.numeric(option_arguments$max_pcs), resolution = resolution, print.output = 0, save.SNN = TRUE)
		PrintFindClustersParams(my_object)
		my_object<-RunTSNE(my_object, dims.use = 1:as.numeric(option_arguments$max_pcs), do.fast = TRUE, reduction.use = reduction_type)

		my_palette<-adjust_palette_size(my_object@ident, basic_color_palette)

		# TSNEPlot(my_object, colors.use = my_palette)
		# TSNEPlot(my_object, colors.use = my_palette, do.label = TRUE, label.size = 10, group.by = "orig.ident")


		plot_title<-paste(option_arguments$projectName, "dim",option_arguments$max_pcs, "res",resolution, sep = "")
		my_tsne_plot<-TSNEPlot(my_object, colors.use = my_palette, do.return=TRUE)
		my_tsne_plot<-my_tsne_plot + ggtitle(plot_title)
		png(paste(option_arguments$projectName, "dim",option_arguments$max_pcs, "res",resolution,"_tsne.png", sep = ""), width = 800, height = 800)
		try(plot(my_tsne_plot))
		dev.off()
		my_tsne_plot<-TSNEPlot(my_object, colors.use = my_palette, do.label = TRUE, label.size = 10, do.return=TRUE)
		my_tsne_plot<-my_tsne_plot + ggtitle(plot_title)
		png(paste(option_arguments$projectName, "dim",option_arguments$max_pcs, "res",resolution,"_tsne-Labeled.png", sep = ""), width = 800, height = 800)
		try(plot(my_tsne_plot))
		dev.off()

		if(resolution == as.numeric(option_arguments$resolutionList[1])){
			my_tsne_plot<-TSNEPlot(my_object, colors.use = my_palette, group.by = "orig.ident", do.return=TRUE)
			my_tsne_plot<-my_tsne_plot + ggtitle(plot_title)
			png(paste(option_arguments$projectName, "dim",option_arguments$max_pcs,"_tsnebyID.png", sep = ""), width = 800, height = 800)
			try(plot(my_tsne_plot))
			dev.off()
			my_tsne_plot<-TSNEPlot(my_object, colors.use = my_palette, do.label = TRUE, label.size = 10, group.by = "orig.ident", do.return=TRUE)
			my_tsne_plot<-my_tsne_plot + ggtitle(plot_title)
			png(paste(option_arguments$projectName, "dim",option_arguments$max_pcs,"_tsnebyID-Labeled.png", sep = ""), width = 800, height = 800)
			try(plot(my_tsne_plot))
			dev.off()
			my_tsne_plot<-TSNEPlot(my_object, colors.use = my_palette, group.by = "Phase", do.return=TRUE)
			my_tsne_plot<-my_tsne_plot + ggtitle(plot_title)
			png(paste(option_arguments$projectName, "dim",option_arguments$max_pcs,"_tsnebyPhase.png", sep = ""), width = 800, height = 800)
			try(plot(my_tsne_plot))
			dev.off()
			my_tsne_plot<-TSNEPlot(my_object, colors.use = my_palette, do.label = TRUE, label.size = 10, group.by = "Phase", do.return=TRUE)
			my_tsne_plot<-my_tsne_plot + ggtitle(plot_title)
			png(paste(option_arguments$projectName, "dim",option_arguments$max_pcs,"_tsnebyPhase-Labeled.png", sep = ""), width = 800, height = 800)
			try(plot(my_tsne_plot))
			dev.off()
		}



		# saveRDS(my_object, file = paste(option_arguments$projectName, "dim",option_arguments$max_pcs, "res",resolution,"_tsne.rds", sep = ""))
		# print(paste("Saved ",option_arguments$projectName, "dim",option_arguments$max_pcs, "res",resolution,"_tsne.rds", sep = ""))
		# Plot cell cycle states --------------------------------------------------------------

		ridgeplot_genes<- c("Pcna", "Top2a", "Mcm6", "Mki67")

		# Finding markers ---------------------------------------------------------

		my_object_markers<-FindAllMarkers(my_object, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)

		my_object_markers_top100 <- my_object_markers %>% group_by(cluster) %>% top_n(100, avg_logFC)
		write.table(my_object_markers_top100, file = paste(option_arguments$projectName, "dim",option_arguments$max_pcs, "res",resolution,"_markers_top100.txt", sep = ""), sep="\t", quote=F)
		my_object_markers_all <- my_object_markers %>% group_by(cluster)
		write.table(my_object_markers_all, file = paste(option_arguments$projectName, "dim",option_arguments$max_pcs, "res",resolution,"_markers_all.txt", sep = ""), sep="\t", quote=F)


		# Print gene plots --------------------------------------------------------

		geneList = c("Nfe2", "Gata1", "Gata2", "Hbg1", "Hbg2", "Zfpm1", "Hbb", "Hba1","Hba2", "Hbd", "Hbe1", "Klf1", "Fli1", "Meis1", "Kit", "Vwf", "Pf4", "Mpo", "Runx1", "Csf1", "Tfr2", "Cnrip1", "Myc" ,"Tk1", "Rrm2", "Itga2b", "Lmna", "Nfia", "Gp1bb", "Plek", "Pbx1", "Hes6", "E2f4","Dntt", "Vpreb1", "Id3", "Atf3", "Jchain", "Cd79a", "Satb1", "Sp140", "Tgfbi", "Lgmn", "Irf8", "Irf7", "Tcf4",  "Batf3", "Tcf19", "Sell", "Cd52", "Hoxb5", "Gata3", "Eno1", "Mpo", "Atp8b4", "Spi1", "Mafk", "Hdc", "Prg2", "Lmo4","Ctsg","Elane", "Cebpa", "Lgals1", "Fosb", "Prtn3",  "Tfrc", "Mpl", "Flt3", "Ca1", "Cd177", "Cd180","Cd244","Cd24","Cd27","Cd34","Cd37","Cd47","Cd48","Cd52","Cd53","Cd63","Cd68","Cd69","Cd72","Cd74","Cd81","Cd82","Cd84","Cd9","Cd93", "Mt1a", "Mt1g", "Mt1f", "Mt1e", "Mt1x")
		housekeepingGenes<-c("Actb", "Gapdh", "Rn18s", "Ppia", "Rpl13a", "Rplp0","B2m", "Hmbs", "Pgk1", "Alas1", "Gusb", "Sdha", "Tbp", "Tubb", "Ywhaz")
		if(option_arguments$genome == "GRCh38"){
			geneList <- toupper(geneList)
			housekeepingGenes<-toupper(housekeepingGenes)
		}

		suppressWarnings(drawmyplot(geneList, my_object, name = paste(option_arguments$projectName, "dim",option_arguments$max_pcs, "res",resolution, sep = "")))
		suppressWarnings(drawmyplot(geneList, my_object, name = paste(option_arguments$projectName, "dim",option_arguments$max_pcs, "res", resolution, "-HKgene", sep = "")))

		# Tabulate data -----------------------------------------------------------

		#Number of cells in each cluster
		table(my_object@ident)
		write.table(table(my_object@ident), file = paste(option_arguments$projectName, "dim",option_arguments$max_pcs, "res",resolution,"_popCounts.txt", sep = ""), row.names = FALSE, quote = FALSE, sep = "\t")


		#Number of cells in each ID
		table(my_object@meta.data$orig.ident)
		write.table(table(my_object@meta.data$orig.ident), file = paste(option_arguments$projectName, "dim",option_arguments$max_pcs, "res",resolution,"_popCounts.txt", sep = ""), row.names = FALSE, quote = FALSE, sep = "\t", append = TRUE)

		#Proportion of cells in each cluster
		prop.table(x = table(my_object@ident))
		write.table(prop.table(x = table(my_object@ident)), file = paste(option_arguments$projectName, "dim",option_arguments$max_pcs, "res",resolution,"_popCounts.txt", sep = ""), row.names = FALSE, quote = FALSE, sep = "\t", append = TRUE)

		#Number of cells from each ID in each cluster
		table(my_object@ident, my_object@meta.data$orig.ident)
		write.table(table(my_object@ident, my_object@meta.data$orig.ident), file = paste(option_arguments$projectName, "dim",option_arguments$max_pcs, "res",resolution,"_popCounts.txt", sep = ""), row.names = FALSE, quote = FALSE, sep = "\t", append = TRUE)

		#Proportion of cells of each ID in each cluster
		prop.table(x = table(my_object@ident, my_object@meta.data$orig.ident))
		write.table(prop.table(x = table(my_object@ident, my_object@meta.data$orig.ident)), file = paste(option_arguments$projectName, "dim",option_arguments$max_pcs, "res",resolution,"_popCounts.txt", sep = ""), row.names = FALSE, quote = FALSE, sep = "\t", append = TRUE)

		#Average expression of each gene in each cluster
		avgexp<-AverageExpression(my_object)
		write.table(avgexp, file = paste(option_arguments$projectName, "dim",option_arguments$max_pcs, "res",resolution,"_AvgXprsn.txt", sep = ""), row.names = TRUE, quote = FALSE, sep = "\t")

	}

	try(saveRDS(my_object, file = paste(option_arguments$projectName, "dim",option_arguments$max_pcs, ".rds", sep = "")), silent = FALSE)
	try(save.image(file = paste(option_arguments$projectName, "dim",option_arguments$max_pcs, ".RData", sep = "")), silent = FALSE)







}
