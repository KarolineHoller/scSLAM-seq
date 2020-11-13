#mitochondrial genes in dr11 names ####
#that file sits on the server in my local/genomes
mito.genes <- read.table("/local/users/aneusch/genomes/mito.genes.vs.txt",sep = ",")
mito.genes <- mito.genes$V3
mito.genes <- as.character(mito.genes)

#prepare seurat-object for the unlabeled data
unlabeled <- Read10X(data.dir = paste0("/local/users/aneusch/experiments/10x/20190306_AB/workdir/20190820_cellranger3_labeled_unlabeled_seurat/remap_unlabeled_Q20/Solo.out_/"))
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
labeled <- Read10X(data.dir = paste0("/local/users/aneusch/experiments/10x/20190306_AB/workdir/20190820_cellranger3_labeled_unlabeled_seurat/remap_labeled_Q20/Solo.out_/"))
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
veg.tomoseq <- read.csv("/local/Karo/scSLAM/6h/vegetal.genes.dR.csv", header=T)

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


