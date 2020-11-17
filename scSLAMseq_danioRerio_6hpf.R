##this script has been run on a system with the following specifications: ####
#R version 3.6.0 (2019-04-26)
#Platform: x86_64-redhat-linux-gnu (64-bit)
#Running under: CentOS Linux 7 (Core)

#Matrix products: default
#BLAS/LAPACK: /usr/lib64/R/lib/libRblas.so

#locale:
# [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

#attached base packages:
#[1] grid      stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
# [1] deconvolveR_1.1   mixtools_1.2.0    gridExtra_2.3     cowplot_1.0.0     ggrepel_0.8.1     reshape2_1.4.3   
# [7] viridis_0.5.1     viridisLite_0.3.0 Seurat_3.1.2      seqinr_3.6-1      ggpubr_0.3.0      ggplot2_3.2.1    

#loaded via a namespace (and not attached):
#  [1] TH.data_1.0-10      Rtsne_0.15          colorspace_1.4-1    ggsignif_0.6.0      rio_0.5.16          ellipsis_0.3.0     
#  [7] ggridges_0.5.1      leiden_0.3.1        listenv_0.8.0       npsurv_0.4-0        mvtnorm_1.0-11      codetools_0.2-16   
# [13] splines_3.6.0       R.methodsS3_1.7.1   mnormt_1.5-5        lsei_1.2-0          TFisher_0.2.0       ade4_1.7-15        
# [19] jsonlite_1.6        broom_0.5.6         ica_1.0-2           kernlab_0.9-29      cluster_2.0.8       png_0.1-7          
# [25] R.oo_1.23.0         uwot_0.1.5          sctransform_0.2.1   compiler_3.6.0      httr_1.4.1          backports_1.1.5    
# [31] assertthat_0.2.1    Matrix_1.2-17       lazyeval_0.2.2      htmltools_0.4.0     tools_3.6.0         rsvd_1.0.2         
# [37] igraph_1.2.4.2      gtable_0.3.0        glue_1.3.1          RANN_2.6.1          dplyr_0.8.3         rappdirs_0.3.1     
# [43] Rcpp_1.0.3          carData_3.0-4       Biobase_2.44.0      cellranger_1.1.0    vctrs_0.3.0         multtest_2.40.0    
# [49] gdata_2.18.0        ape_5.3             nlme_3.1-139        gbRd_0.4-11         lmtest_0.9-37       stringr_1.4.0      
# [55] globals_0.12.5      openxlsx_4.1.5      lifecycle_0.2.0     irlba_2.3.3         gtools_3.8.1        rstatix_0.5.0      
# [61] future_1.15.1       MASS_7.3-51.4       zoo_1.8-6           scales_1.1.0        hms_0.5.3           parallel_3.6.0     
# [67] sandwich_2.5-1      RColorBrewer_1.1-2  curl_4.3            reticulate_1.14     pbapply_1.4-2       segmented_1.1-0    
# [73] stringi_1.4.3       mutoss_0.1-12       plotrix_3.7-7       caTools_1.17.1.3    BiocGenerics_0.30.0 zip_2.0.4          
# [79] bibtex_0.4.2.2      Rdpack_0.11-1       SDMTools_1.1-221.2  rlang_0.4.6         pkgconfig_2.0.3     bitops_1.0-6       
# [85] lattice_0.20-38     ROCR_1.0-7          purrr_0.3.3         htmlwidgets_1.5.1   tidyselect_1.1.0    RcppAnnoy_0.0.14   
# [91] plyr_1.8.5          magrittr_1.5        R6_2.4.1            generics_0.0.2      gplots_3.0.1.1      multcomp_1.4-11    
# [97] haven_2.3.0         foreign_0.8-71      pillar_1.4.3        withr_2.1.2         sn_1.5-4            fitdistrplus_1.0-14
#[103] abind_1.4-5         survival_2.44-1.1   tibble_3.0.1        future.apply_1.3.0  tsne_0.1-3          car_3.0-8          
#[109] crayon_1.3.4        KernSmooth_2.23-15  plotly_4.9.1        readxl_1.3.1        data.table_1.12.8   forcats_0.5.0      
#[115] metap_1.2           digest_0.6.23       tidyr_1.0.0         numDeriv_2016.8-1.1 R.utils_2.9.2       RcppParallel_4.4.4 
#[121] stats4_3.6.0        munsell_0.5.0      

##load environment ####
library(Seurat)
library(viridis)
library(ggplot2)
library(reshape2)
library(ggrepel)
library(cowplot)
library(grid)
library(gridExtra)
library(mixtools)

##sanity checks for our protocol of scSLAM-seq in live embryos, code by @anikaneuschulz####
#calculate T-to-C mutation per gene
#the output of the script mutation_eff_per_gene_v4_MTabs_final.py produces the input for the following plot:
per_gene <- read.csv("/your/file/path/possorted_genome_bam_actual_cells_MD_MTgene_specific_mutation_rates.csv",
                     header = TRUE, sep = "\t")
per_gene <- subset(per_gene, readcount > 99)
per_gene$condition <- "labeled"

per_gene_unlabeled <- read.csv("/your/file/path/Aligned.out_MT.bam.featureCountsgene_specific_mutation_rates.csv",
                               header = TRUE, sep ="\t")
per_gene_unlabeled <- subset(per_gene_unlabeled, readcount > 99)
per_gene_unlabeled$condition <- "unlabeled"

per_gene_both <- rbind(per_gene,per_gene_unlabeled)

#the following code produces the plot in Fig. 3c
per_gene_plot <- ggplot(data = per_gene_both, mapping = aes(x=per_gene_both$TC*100, fill = condition, alpha = condition)) +
  geom_histogram(bins = 75, position="identity") +
  xlab("T to C frequency [%]") +
  ylab("number of genes") +
  scale_fill_manual(values = c("#DDAA33", "#004488")) +
  scale_alpha_manual(values = c(1,0.5)) +
  theme_linedraw()

#the output of script 4 produces the input for the following code chuck, that we used to produce the plot in Figure 3b (base mutation frequencies):
mutations <- "/your/file/path/unlabeled_labeled_2mutations_Q20_sorted_mutation_occurences_20_0N.txt"
nucleotides <- "/your/file/path/unlabeled_labeled_2mutations_Q20_sorted_nucleotide_counts_Q20_C0N.txt"

mutations_unlabeled <- "/your/file/path/Aligned.out_MT_mutation_occurences_20_0N.txt"
nucleotides_unlabeled <- "/your/file/path/Aligned.out_MT_nucleotide_counts_Q20_C0N.txt"

mutation_occurences <-read.table(file = mutations, sep = '\t', header = FALSE)
nucleotide_numbers <- read.table(file = nucleotides, sep = "\t", header = FALSE)

mutation_occurences_ul <- read.table(file = mutations_unlabeled, sep = '\t', header = FALSE)
nucleotide_numbers_ul <- read.table(file = nucleotides_unlabeled, sep = "\t", header = FALSE)#yes, this is idiotic. It didn't give errors though

mutationTableA <- subset(mutation_occurences, V1=="AC")
mutationTableA <- rbind(mutationTableA, subset(mutation_occurences, V1=="AG"))
mutationTableA <- rbind(mutationTableA, subset(mutation_occurences, V1=="AT"))
mutationTableT <- subset(mutation_occurences, V1=="TC")
mutationTableT <- rbind(mutationTableT, subset(mutation_occurences, V1=="TG"))
mutationTableT <- rbind(mutationTableT, subset(mutation_occurences, V1=="TA"))
mutationTableC <- subset(mutation_occurences, V1=="CT")
mutationTableC <- rbind(mutationTableC, subset(mutation_occurences, V1=="CG"))
mutationTableC <- rbind(mutationTableC, subset(mutation_occurences, V1=="CA"))
mutationTableG <- subset(mutation_occurences, V1=="GT")
mutationTableG <- rbind(mutationTableG, subset(mutation_occurences, V1=="GC"))
mutationTableG <- rbind(mutationTableG, subset(mutation_occurences, V1=="GA"))

mutationTableA$Freq <- (mutationTableA$V2/subset(nucleotide_numbers, V1 == "A")$V2)*100
mutationTableT$Freq <- (mutationTableT$V2/subset(nucleotide_numbers, V1 == "T")$V2)*100
mutationTableC$Freq <- (mutationTableC$V2/subset(nucleotide_numbers, V1 == "C")$V2)*100
mutationTableG$Freq <- (mutationTableG$V2/subset(nucleotide_numbers, V1 == "G")$V2)*100
mutationFrequencyTable <- rbind(mutationTableA, mutationTableC, mutationTableG, mutationTableT)
mutationFrequencyTable$condition <- "Q20 labeled"

# unlabeled data
mutationTableA_ul <- subset(mutation_occurences_ul, V1=="AC")
mutationTableA_ul <- rbind(mutationTableA_ul, subset(mutation_occurences_ul, V1=="AG"))
mutationTableA_ul <- rbind(mutationTableA_ul, subset(mutation_occurences_ul, V1=="AT"))
mutationTableT_ul <- subset(mutation_occurences_ul, V1=="TC")
mutationTableT_ul <- rbind(mutationTableT_ul, subset(mutation_occurences_ul, V1=="TG"))
mutationTableT_ul <- rbind(mutationTableT_ul, subset(mutation_occurences_ul, V1=="TA"))
mutationTableC_ul <- subset(mutation_occurences_ul, V1=="CT")
mutationTableC_ul <- rbind(mutationTableC_ul, subset(mutation_occurences_ul, V1=="CG"))
mutationTableC_ul <- rbind(mutationTableC_ul, subset(mutation_occurences_ul, V1=="CA"))
mutationTableG_ul <- subset(mutation_occurences_ul, V1=="GT")
mutationTableG_ul <- rbind(mutationTableG_ul, subset(mutation_occurences_ul, V1=="GC"))
mutationTableG_ul <- rbind(mutationTableG_ul, subset(mutation_occurences_ul, V1=="GA"))

mutationTableA_ul$Freq <- (mutationTableA_ul$V2/subset(nucleotide_numbers_ul, V1 == "A")$V2)*100
mutationTableT_ul$Freq <- (mutationTableT_ul$V2/subset(nucleotide_numbers_ul, V1 == "T")$V2)*100
mutationTableC_ul$Freq <- (mutationTableC_ul$V2/subset(nucleotide_numbers_ul, V1 == "C")$V2)*100
mutationTableG_ul$Freq <- (mutationTableG_ul$V2/subset(nucleotide_numbers_ul, V1 == "G")$V2)*100
mutationFrequencyTable_ul <- rbind(mutationTableA_ul, mutationTableC_ul, mutationTableG_ul, mutationTableT_ul)
mutationFrequencyTable_ul$condition <- "Q20 unlabeled"

mutationFrequencyTable_all <- rbind(mutationFrequencyTable_ul, 
                                    mutationFrequencyTable)

#here comes Figure 3c:
ggplot(data = mutationFrequencyTable_all, aes(x=mutationFrequencyTable_all$V1, y=mutationFrequencyTable_all$Freq, color = condition)) +
  geom_point() +
  xlab("")+ 
  ylab("mutation frequency in % of base") +
  scale_color_manual( values = c("#ddaa33","#004488" )) +
  theme_linedraw()

#for plot S3c, run the script UMIs_per_cell+genes_per_cell_final.py and use the output for the following ggplot
umi_gene_stats <- read.csv("/your/file/path/...UMIs+genes_per_cell.csv", header=T, sep = "\t")

ggplot(data = umi_gene_stats, mapping = aes(x = genes)) +
  geom_histogram(bins = 200) +
  ylab("number of cells") +
  theme_linedraw()
#for plot in S3d, use the output of script 3
UMI_stats_rep1 <- read.table("/your/file/path/possorted_genome_bam_actual_cells_MD_MT_umi_stats.csv",
                             header = TRUE,
                             sep = "\t")

ggplot(data = subset(UMI_stats_rep1, UMI_stats_rep1$nReads < 50), mapping = aes(x = nReads)) +
  geom_histogram(binwidth = 1) +
  xlab("Reads per UMI") +
  ylab("number of occurences") +
  theme_linedraw()

##Now let's have a look at the single cell data! ####
#mitochondrial genes with annotation of dR11 ####
mito.genes <- read.table("/your/file/path/mito.genes.vs.txt",sep = ",")
mito.genes <- mito.genes$V3
mito.genes <- as.character(mito.genes)

#prepare seurat-object for the unlabeled data
unlabeled <- Read10X(data.dir = paste0("/your/file/path/Solo.out_/"))
unlabeled <- CreateSeuratObject(counts = unlabeled,
                                min.cells = 3, min.features = 150,
                                project = "unlabeled")
mito.genes.use <- setdiff(mito.genes,setdiff(mito.genes,rownames(unlabeled[["RNA"]])))
unlabeled[["percent.mito"]] <- PercentageFeatureSet(object = unlabeled, features = mito.genes.use)
unlabeled <- subset(x = unlabeled, subset = percent.mito < 25 & nFeature_RNA < 3500)

unlabeled <- NormalizeData(object = unlabeled, normalization.method = "LogNormalize", scale.factor = 1e4)

unlabeled <- FindVariableFeatures(object = unlabeled,selection.method = 'vst', nfeatures = 2000)

unlabeled <- ScaleData(unlabeled, verbose = TRUE, vars.to.regress = c("nFeature_RNA","nCount_RNA","percent.mito")) #, features = all.genes)

unlabeled <- RunPCA(unlabeled, npcs = 100, verbose = TRUE)

unlabeled <- RunUMAP(unlabeled, reduction = "pca", dims = 1:20)

unlabeled <- FindNeighbors(unlabeled, reduction = "pca", dims = 1:20)

unlabeled <- FindClusters(unlabeled, resolution = 0.5)

## prepare Seurat-Object for labeled reads ####
labeled <- Read10X(data.dir = paste0("/your/file/path/Solo.out_/"))
labeled <- CreateSeuratObject(counts = labeled,
                              min.cells = 3, min.features = 150,
                              project = "labeled")
mito.genes.use <- setdiff(mito.genes,setdiff(mito.genes,rownames(labeled[["RNA"]])))
labeled[["percent.mito"]] <- PercentageFeatureSet(object = labeled, features = mito.genes.use)
labeled <- subset(x = labeled, subset = percent.mito < 25 & nFeature_RNA < 3500)

mito.genes.use <- setdiff(mito.genes,setdiff(mito.genes,rownames(labeled[["RNA"]])))

labeled[["percent.mito"]] <- PercentageFeatureSet(object = labeled, features = mito.genes.use)

labeled <- subset(x=labeled, subset=percent.mito < 25 & nFeature_RNA < 3500)

labeled <- NormalizeData(object = labeled, normalization.method = "LogNormalize", scale.factor = 1e4)

labeled <- FindVariableFeatures(object = labeled,selection.method = 'vst', nfeatures = 2000)

labeled <- ScaleData(labeled, verbose = TRUE, vars.to.regress = c("nFeature_RNA","nCount_RNA","percent.mito"))

labeled <- RunPCA(labeled, npcs = 100, verbose = TRUE)

labeled <- RunUMAP(labeled, reduction = "pca", dims = 1:20)
labeled <- FindNeighbors(labeled, reduction = "pca", dims = 1:20)
labeled <- FindClusters(labeled, resolution = 0.5)

labeled.markers <- FindAllMarkers(object = labeled, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

## prepare labeling efficiency per cell plot ####
labeling_eff <- read.table("/your/file/path/mutation_rate_per_cell.tsv",
                           sep = "\t", header = TRUE)
labeling_eff$CB <- sapply(labeling_eff$CB, function(x) substr(x, 1, 16))
rownames(labeling_eff) <- labeling_eff$CB
labeling_eff$TtoCr <- sapply(labeling_eff$TtoC, function(x) round(x, 3))
labeled@meta.data <- merge(labeled@meta.data,labeling_eff,by="row.names",all.x=T)

rownames(labeled@meta.data) <- labeled@meta.data$Row.names
labeled@meta.data$Row.names <- NULL

FeaturePlot(labeled, features = "TtoCr", cols = viridis(4)) #the labeling efficiency goes up to 4 percent, so we chose 4 colors from the viridis palette


#define cluster identities as identified by marker gene expression labeled.markers and unlabeled.markers
labeled_renamed <- labeled
labeled_renamed <- RenameIdents(object = labeled, `0` = "undetermined", `1` = "ectoderm (neuroectoderm)", `2` = "non-axial mesoderm",
                                `3` = "ectoderm (neuroectoderm)", `4` = "non-axial mesoderm, ventral mesoderm",
                                `5` = "axial mesoderm, dorsal organizer", `6` = "undetermined", `7` = "undetermined",
                                `8` = "non-axial mesoderm", `9` = "prospective neural plate", `10` = "axial mesoderm, dorsal organizer",
                                `11` = "enveloping layer", `12` = "undetermined", `13` = "prospective primordial germ cells")
##stash new identities
labeled_renamed[["named_clusters"]] <- Idents(object = labeled_renamed)

# add new identities to old RNA clusters
labeled_identites <- labeled_renamed@meta.data
labeled_identites[,1:6] <- NULL #remove all clumns from metadata despite the named cluster identity

unlabeled@meta.data <- merge(unlabeled@meta.data, labeled_identites, by="row.names",all.x=T)
#bring back row.names, as they get mogmled by merge
rownames(unlabeled@meta.data) <- unlabeled@meta.data$Row.names
unlabeled@meta.data$Row.names <- NULL

#Check clustering with DimPlots
karos.col <- brewer.pal(9, "Paired") #number of colors depends on number of clusters you expect
DimPlot(unlabeled, group.by = 'named_clusters', cols=karos.col) +
  ggtitle("clustering on unlabeled RNA") + 
  theme(text = element_text(size=12), axis.text = element_text(size=11), legend.position = "none") +
  xlab("UMAP1") +
  ylab("UMAP2")

DimPlot(unlabeled, group.by = 'seurat_clusters', cols=karos.col) +
  ggtitle("clustering on unlabeled RNA") + 
  theme(text = element_text(size=12), axis.text = element_text(size=11), legend.position = "none") +
  xlab("UMAP1") +
  ylab("UMAP2")

#
unlabeled.markers.lc <- FindAllMarkers(object = unlabeled, only.pos = TRUE)
#lowered cutoff: the parameter min.pct can be set to zero here because the fraction of recovered PGCs is very low - also known marker genes are expressed
#in just a small fraction
                             
#plot heatmaps for the datasets respective marker genes (Figure 3e)
#downsample to the number of prospective PGCs in the labeled_renamed dataset
length(WhichCells(labeled_renamed, idents = "prospective primordial germ cells"))
#54
DoHeatmap(subset(labeled_renamed, downsample=54), features = labeled.markers$gene, size=3,
          group.colors = karos.col) +
  theme(axis.text.x=element_blank(), axis.text.y=element_blank()) +
  labs(title = "marker gene expression in labeled data")

#cluster 7 marker genes as a dataframe
cluster7 <- unlabeled.markers.lc[unlabeled.markers.lc$cluster==7,]
                            
#Compare cluster 7 marker genes with vegetally localized genes in the one-cell stage embryo
#for that, load partial results of a tomo-seq dataset
veg.tomoseq <- read.csv("/your/file/path/vegetal.genes.dR.csv", header=T)

g=intersect(veg.tomoseq$Gene, cluster7$gene)
length(g) #28 - the overlap between vegetally localized genes and cluster 7 marker genes ('prospective PGCs')
g

###what is the average expression of the vegetally localised genes in each cluster?
#1. how many and which genes are expressed in the single cell dataset?
h=as.vector(veg.tomoseq$Gene)
length(h)

h %in% unlabeled@assays$RNA@counts@Dimnames[[1]]
sum(h %in% unlabeled@assays$RNA@counts@Dimnames[[1]])
#91 of the 97 vegetally localized genes are found back in the single cell dataset

Idents(unlabeled) = unlabeled@meta.data$named_clusters
k=h[h %in% unlabeled@assays$RNA@counts@Dimnames[[1]]]

test <- AverageExpression(unlabeled, features = k)
expr.veg2 <- test$RNA
expr.veg2$gene <- row.names(expr.veg2)

#calculating the mean expression over all cells of the dataset: 1. set identity to a unique value 2. calculate average
# 3. reset identity to seurat_clustering
unlabeled$ident1 <- "1"
Idents(unlabeled) <- "ident1"
av <- AverageExpression(unlabeled, features = expr.veg2$gene)
Idents(unlabeled) <- "named_clusters"
expr.veg2$mean.SLAM <- av$RNA[,1]

expr.veg3 <- expr.veg2[expr.veg2$mean.SLAM>0.1,] #filters out very lowly expressed genes (fold changes can get very noisy and are not meaningful to interpret. cutoff depends on your biological sample 
unlab.veg.FC <- as.data.frame(matrix(ncol=1, nrow = 47))
unlab.veg.FC$gene <- expr.veg3$gene

unlab.veg.FC$`undetermined` <- expr.veg3$undetermined/expr.veg3$mean.SLAM
unlab.veg.FC$`ectoderm (neuroectoderm)` <- expr.veg3$`ectoderm (neuroectoderm)`/expr.veg3$mean.SLAM
unlab.veg.FC$`non-axial mesoderm` <- expr.veg3$`non-axial mesoderm`/expr.veg3$mean.SLAM
unlab.veg.FC$`non-axial mesoderm, ventral mesoderm` <- expr.veg3$`non-axial mesoderm, ventral mesoderm`/expr.veg3$mean.SLAM
unlab.veg.FC$`axial mesoderm, dorsal organizer` <- expr.veg3$`axial mesoderm, dorsal organizer`/expr.veg3$mean.SLAM
unlab.veg.FC$`prospective neural plate` <- expr.veg3$`prospective neural plate`/expr.veg3$mean.SLAM
unlab.veg.FC$`enveloping layer` <- expr.veg3$`enveloping layer`/expr.veg3$mean.SLAM
unlab.veg.FC$`prospective primordial germ cells` <- expr.veg3$`prospective primordial germ cells`/expr.veg3$mean.SLAM

#plot the FC of filtered vegetal genes in different cluster idents:
melt.unlab.veg.FC <- melt(unlab.veg.FC[, -1], id="gene")  #melt the dataframe for plotting
colnames(melt.unlab.veg.FC)[2] <- "cluster"

ggplot(melt.unlab.veg.FC, aes(x=cluster, y = value))+
  theme_linedraw() +
  geom_jitter(alpha= 0.5, width = 0.2) +
  stat_summary(fun.y= mean, fun.ymin=mean, fun.ymax=mean, geom="crossbar", width=0.5, color="#CE3A39") +
  geom_hline(yintercept = 1, linetype= "dashed") + 
  #geom_text(label = melt.unlab.veg.FC$gene, check_overlap = T) +
  theme(text = element_text(size=15), axis.text.x=element_text(angle=60, hjust=1)) +
  labs(title = 'clusterwise enrichment of vegetally localized genes at 6 hpf', y= 'fold change', x= "")

#plot histogram: compare the log2 mean FC of the vegetally localised genes to randomly selected genes. For that,
#select the same number of random genes and calculate their log(meanFC). Iterate 1000 times to calculate a p-value for the hypothesis
# that the log2FC of vegetally localized genes is significantly higher than one expects for random genes on average. (not shown in the manuscript)

#filters genes with an average expression > 0.1 #apply the same filters for randomly sampled genes than for vegetally localized genes
Idents(unlabeled) <- "ident1"
test <- AverageExpression(unlabeled)
unlabeled.filt <- test$RNA
unlabeled.filt$gene <- rownames(unlabeled.filt)
colnames(unlabeled.filt)[1] <- "mean.SLAM"
unlabeled.filt <- unlabeled.filt[ unlabeled.filt$mean.SLAM>0.1,]
Idents(unlabeled) <- "named_clusters"

#filters out genes that are not expressed in PGCs
test <- AverageExpression(unlabeled, features = unlabeled.filt$gene)
av.random <- test$RNA
av.random$gene <- rownames(av.random)
av.random <- av.random[!av.random$`prospective primordial germ cells` == 0,]
unlab.filt <- merge(unlabeled.filt, av.random, all.y = T)

#initiate objects for the loop
random.av <- NULL
m <- as.data.frame(matrix(nrow = 47))
mean.distr <- as.vector(NULL)

for (i in 1:1000){
  random.genes <- unlab.filt[sample(nrow(unlab.filt)-1, 47),]
  m$FC = random.genes$`prospective primordial germ cells`/random.genes$mean.SLAM
  m$gene <- random.genes$gene
  d = log2(mean(m$FC))
  #print(m)
  mean.distr <- c(mean.distr, d)
  random.av <- rbind(random.av, m)
}

mean(random.av$FC) #2.1   
log2(mean(random.av$FC))  #= 1.07
max(random.av$FC)
log2(mean(unlab.veg.FC$`prospective primordial germ cells`)) #calculates the log2 mean FC of filtered vegetally localized genes in the PGC cluster
#1.58

#shows the distribution of log2 FC from 1000 iterations as a histogram
hist(mean.distr, breaks=50, main="distribution of log2 mean FC of 1000 iterations")
sum(mean.distr>1.58) #34
#p-value of log2 mean FC of vegetally localised genes (1.58) is 34/1000=0.034

random.av$log2FC <- log2(random.av$FC)
unlab.veg.FC$log2FC <- log2(unlab.veg.FC$`prospective primordial germ cells`)

#plots ALL log2FC of the 1000 iterations in comparison to the distribution of log2FC of the vegetally localised genes
ggplot() +
  theme_linedraw() +
  geom_density(data= random.av, aes(x=log2FC), fill="darkgray", alpha = 0.6) +
  geom_density(data= unlab.veg.FC, aes(x = log2FC), fill= "#FF7F00", alpha = 0.6) +
  geom_vline(xintercept = 0.92, colour= "darkgray", linetype ="dashed") +
  geom_vline(xintercept = 1.58, colour= "#FF7F00", linetype ="dashed") +
  scale_x_continuous(limits = c(- 5, 8)) +
  labs(title = "log2FC of 47 sampled genes in PGCs", x = "log2FC")


#fit a mixed model to the bimodal distribution if the PGC cluster:

veg.distribution <-  log2(unlab.veg.FC$`prospective primordial germ cells`)
mixmdl <- normalmixEM(veg.distribution)
plot(mixmdl,which=2, breaks =15)
lines(density(veg.distribution), lty=2, lwd=2, col="#2B3584")
curve(mixfunction, -2, 4, type ="l")
summary(mixmdl)
#          comp 1   comp 2
#lambda  0.281843 0.718157
#mu     -0.561875 1.521403
#sigma   0.390803 1.123669
#we can manually adjust the values with the result of two better fitted normal distributions:
mixfunction <- function(x) {
  0.27*dnorm(x, mean=-0.5, sd=0.65) + 0.73*dnorm(x, mean=1.52, sd=1.3)
}

norm1.function <- function(x){
  0.27*dnorm(x, mean = -0.5, sd = 0.65)
}
norm2.function <- function(x) {
  0.73*dnorm(x, mean=1.52, sd=1.3)
  }

plot(mixfunction, -2, 5, type ="l", ylab = "")
lines(density(veg.distribution), lty=2, lwd=2, col ="#2B3584")
curve(norm2.function, -2, 5, type ="l", col = "blue", ylab ="", add =T)

#add the fitted normal distribution of germ cell factor candidates to the random distribution and the log2FC of the vegetal genes in PGCs
ggplot(data= unlab.veg.FC) +
  theme_linedraw() +
  stat_density(data= random.av, aes(x=log2FC), position= "stack", fill="darkgray", alpha = 0.5, col="darkgray") +
  stat_function(fun = norm1.function, col = "#00A6DD") +
  stat_function(fun = norm2.function, col = "#2B3584") +
  stat_density(aes(x = log2FC), position= "stack", geom="line", col= "black", linetype= "dashed") +
  geom_vline(xintercept = 0.4, colour= "darkgray", linetype ="dotted") +
  geom_vline(xintercept = 1.52, colour= "#2B3584", linetype ="dotted") +
  scale_x_continuous(limits = c(- 5, 8)) +
  labs(title = "vegetal genes in PGCs fall into 2 populations", x = "log2FC")

#calculate the p-value of the norm2.function to the random (background) distribution
#1. calculate Welchs t-statistic and degree of freedom (by hand)
#t-value= 3.765, degree of freedom df= 73.79
#2. use pt() function in R to get p-value
pt(3.765, df=72.79, lower.tail=F)
#0.000167467 = 1.7*10^-4
                             
#plot the Average Expression for in PGCs enriched genes as a barplot (Figure 4c):
test <- AverageExpression(unlabeled, features = c("ppp1r3b", "sh2d5", "anln", "itpkca", "ndel1b", "krtcap2"))
plot.data <- test$RNA
plot.data$gene <- row.names(plot.data)
melt.plot.data <- melt(plot.data, id.vars = "gene")
colnames(melt.plot.data)[2:3] <- c("cell.type", "Expression")
                             
melt.plot.data$gene <- as.factor(melt.plot.data$gene)
for(i in 1:6){
  gene_col <- levels(melt.plot.data$gene)
  p11 <- ggplot(data=melt.plot.data[melt.plot.data$gene == gene_col[i],])+
    theme_linedraw() + 
    geom_bar(aes(x= cell.type, y= Expression, fill= gene), stat= "identity") +
    theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
    scale_fill_manual(values = c('black')) 
  print(p11)
}
                             
#plot the average expression of classical germ cell factor genes on a UMAP:
FeaturePlot(unlabeled, features = c("h1m", "dnd1", "buc", "nanos3"), cols = c('lightgrey', 'blue4'))
                             
#calculate the fold change of expression in PGC cells between labeled and unlabeled reads for maternally vegetally loclaized
#genes that are above the expression cutoff of the dataset (Figure S4c):
q= expr.veg3$gene
                             
#dataset-wise, set the identifier to an equal value
labeled.ren.q <-AverageExpression(labeled_renamed, features = q, slot = "counts")
unlabeled.q <-AverageExpression(unlabeled, features = q, slot = "counts")

plotq <- labeled.ren.q$RNA$`prospective primordial germ cells`/unlabeled.q$RNA$`prospective primordial germ cells`
                             
hist(plotq, breaks = 30, col= "darkgrey", main = "detection of vegetally localized genes in PGCs", xlab = "fold change labeled/unlabeled counts")
