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

#mitochondrial genes in with annotation of dR11 ####
#that file sits on the server in my local/genomes
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

#define cluster identities
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
# row.names get mingled by merge, therefore we have to bring them back
rownames(unlabeled@meta.data) <- unlabeled@meta.data$Row.names
unlabeled@meta.data$Row.names <- NULL

#Check clustering with DimPlots
karos.col <- brewer.pal(9, "Paired")
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
#the parameter min.pct can be set to zero here because the fraction of PGCs is very low - also known marker genes are expressed
#in just a small fraction

#cluster 7 marker genes as a dataframe
cluster7 <- unlabeled.markers.lc[unlabeled.markers.lc$cluster==7,]
veg.tomoseq <- read.csv("/your/file/path/vegetal.genes.dR.csv", header=T)

g=intersect(veg.tomoseq$Gene, cluster7$gene)
length(g) #28 - the overlap between vegetally localized genes and cluster 7 marker genes ('prospective PGCs')
g

###what is the average expression of the vegetally localised genes in each cluster?
#1. which genes are expressed in the single cell dataset?
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
#expr.veg <- merge(plot.av.veg, veg.tomoseq, by= 'gene', all.x=T)
#expr.veg$mean.tomo <- rowMeans(expr.veg[,c(11,13,15)])

#calculating the mean expression over all cells: 1. set identity to a unique value 2. calculate average
# 3. reset identity to seurat_clustering
unlabeled$ident1 <- "1"
Idents(unlabeled) <- "ident1"
av <- AverageExpression(unlabeled, features = expr.veg2$gene)
Idents(unlabeled) <- "named_clusters"
expr.veg2$mean.SLAM <- av$RNA[,1]

expr.veg3 <- expr.veg2[expr.veg2$mean.SLAM>0.1,]
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

#plot the FC of filteres vegetal genes in different cluster idents:
melt.unlab.veg.FC <- melt(unlab.veg.FC[, -1], id="gene")
colnames(melt.unlab.veg.FC)[2] <- "cluster"

ggplot(melt.unlab.veg.FC, aes(x=cluster, y = value))+
  theme_linedraw() +
  geom_jitter(alpha= 0.5, width = 0.2) +
  stat_summary(fun.y= mean, fun.ymin=mean, fun.ymax=mean, geom="crossbar", width=0.5, color="#CE3A39") +
  geom_hline(yintercept = 1, linetype= "dashed") + 
  #geom_text(label = melt.unlab.veg.FC$gene, check_overlap = T) +
  theme(text = element_text(size=15), axis.text.x=element_text(angle=60, hjust=1)) +
  labs(title = 'clusterwise enrichment of vegetally localised genes at 6 hpf', y= 'fold change', x= "")

#histogram plot: compare the log2 mean FC of the vegetally localised genes to randomly selected genes. For that,
#select the same number of random genes and calculate their log(meanFC)

#filters genes with an average expression > 0.1
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
log2(mean(unlab.veg.FC$`prospective primordial germ cells`)) #calculates the log2 mean FC of filtered veg genes in the PGC cluster
#1.58

#shows the distribution of log2 FC from 1000 iterations as a histogram
hist(mean.distr, breaks=50, main="distribution of log2 mean FC of 1000 iterations")
sum(mean.distr>1.58) #34
#p-value of log2 mean FC of vegetally localised genes (1.58) is 34/1000=0.034
#why is this not centred around 0?

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


