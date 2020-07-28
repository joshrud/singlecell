############# scRNAseqAnalysis.R #######################
#### processes aligned scRNA-seq files using Seurat

##### NOTES ####
# This file should be supplemented by a list of functions used by the processes in this script (something like "scNRAseqFunctions.R")
# Not sure how much control over the pipeline we want in this script -- alignment? 

################### Packages ###########################
library(Seurat)
library(scatterpie)
library(limma)
library(ggrepel)
library(Seurat)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(viridis)
library(reshape2)
library(RColorBrewer)
library(cowplot)
library(future)
library(SingleR)
library(matrixStats)
library(monocle3)

################### Sources ############################
scriptPath <- rstudioapi::getSourceEditorContext()$path #only works from rstudio (haven't tested on server)
source(paste0(scriptPath, "/scRNAseqFunctions.R")) #should just be in the same folder as this script 

################### Settings ###########################
options(future.globals.maxSize = 8000 * 1024^2)
runQC = TRUE
DEMULTIPLEX = TRUE #if [de/free]muxlet was run to demultiplex samples
multiLane = FALSE #if 2 lanes were run for a sample (or all samples)
dataDir = "" #wherever the folder of sample folders of matrix files are 
SPECIES = "UNICORN"


################### Constants ###########################
DARKCOLORS <- c(RColorBrewer::brewer.pal(8,"Set1"), RColorBrewer::brewer.pal(8,"Dark2"), "#FFFFFF")
DEFAULT_RES <- 0.5 #subject to change
LIGHTGRAY = "#D4D3D2"
MAX_PERCENT_MT = 5
MAX_PERCENT_RIBO = 5
MAX_PERCENT_HEMO = 5
MIN_CELLS_PER_FEATURE = 3
MIN_FEATURES_PER_CELL = 100
MIN_COUNTS_PER_CELL = 100
MAX_FEATURES_PER_CELL = 2500
HEMOGENES.MMU <- c("Hbb-bt", "Hbb-bs", "Hba-a1", "Hba-a2", "Hbq1a", "Hbq1b")
HEMOGENES.HSA <- c("HBB", "HBA1", "HBA2")
VAR.FEATURES <- 3000 #num variable genes

################### Analysis set up #####################
if (SPECIES == "HUMAN") {
  HEMOGENES = HEMOGENES.HSA
} else if (SPECIES == "MOUSE") {
  HEMOGENES = HEMOGENES.MMU 
} else {
  print("only mmu and hsa supported for species")
}

# importing matrix.mtx, setting up Experiment and Creating objects
# get sample names from matrix path x/y/sampleName/matrix.mtx
## import all the data
seurat_project <- list()
for (f in list.files(dataDir)) { #changed to .jfgc
  print(paste0(f, "/", list.files(paste0(dataDir, "/", f, "/Solo.out")))) #assuming we use STAR solo
  curfolder <- paste0(dataDir, "/", f,"/Solo.out")
  if (length(list.files(curfolder)) > 0) {
    print(paste0("importing: ", f))
    temp.data <- Read10X(curfolder)
    seurat_project[[f]] <- CreateSeuratObject(counts = temp.data, project = f, min.cells = MIN_CELLS_PER_FEATURE, min.features = MIN_FEATURES_PER_CELL)
  }
}

if (DEMULTIPLEX) {
  demult.dir <- "demultiplexing" #directory with clust1.samples.gz files
  demultfiles <- list.files(demult.dir)
  demultfiles <- paste0(demult.dir, "/", demultfiles[grep("clust1.samples.gz", demultfiles)])
  for (i in 1:length(demultfiles)) {
    curfile <- read.table(gzfile(demultfiles[i]), sep = "\t")
    colnames(curfile) <- as.character(as.matrix(curfile[1,]))
    curfile <- curfile[-1,]
    if (i == 1) {
      demult_cell_table <- curfile
    } else {
      demult_cell_table <- rbind(demult_cell_table, curfile)
    }
  }
}

lapply(seurat_project, pp_plots)
seurat_project <- lapply(seurat_project, pp, subset=FALSE) 


# saveRDS(seurat_project, "saves/seurat_project_loaded_beforeIntegration.rds")
seurat_project.anchors <- FindIntegrationAnchors(object.list = seurat_project, dims = 1:30)
seurat_project.int <- IntegrateData(anchorset = seurat_project.anchors, dims = 1:30)
seurat_project.int <- ScaleData(seurat_project.int) #might not want to include nUMI in this because some celltypes don't have as many rna molecules
seurat_project.int <- RunPCA(object = seurat_project.int, features = VariableFeatures(object = seurat_project.int))
seurat_project.int <- RunUMAP(object = seurat_project.int, dims=1:20)
seurat_project.int <- FindNeighbors(object = seurat_project.int, dims=1:20)
seurat_project.int <- FindClusters(object = seurat_project.int, dims=1:20, resolution = 0.4)
saveRDS(seurat_project.int, "seurat_project_integrated_processed.rds")
# med12pos_myo_int_final@meta.data$condition <- "med12pos"
# med12pos_myo_int_final@meta.data$condition[grep("myometrium", med12pos_myo_int_final@meta.data$orig.ident)] <- "myometrium"
# med12pos_myo_int_final
seurat_project.int.markers <- FindAllMarkers(seurat_project.int, min.pct = 0.25, logfc.threshold = 0.25, random.seed = 777)
write.table(seurat_project.int.markers, "seurat_project_int_allMarkers.txt", sep = "\t", quote = F, col.names = NA)
saveRDS(seurat_project.int.markers, "seurat_project_int_allMarkers.rds")
top10 <- seurat_project.int.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(seurat_project.int, features = as.vector(top10$gene)) 
height80dpi <- (length(as.vector(top10$gene))*12)*(1/80)+1
ggsave("heatmap_allMarkers.txt", height= height80dpi, width = 8)

#remove all hemo-containing cells (maybe use some sort of quantile cut-off here?)
seurat_project.int.noHBB <- subset(seurat_project.int, subset=percent.hbb == 0)
seurat_project.int.noHBB <- ScaleData(seurat_project.int.noHBB)
seurat_project.int.noHBB <- RunPCA(object = seurat_project.int.noHBB, features = VariableFeatures(object = seurat_project.int.noHBB))
seurat_project.int.noHBB <- RunUMAP(object = seurat_project.int.noHBB, dims=1:20)
seurat_project.int.noHBB <- FindNeighbors(object = seurat_project.int.noHBB, dims=1:20)
seurat_project.int.noHBB <- FindClusters(object = seurat_project.int.noHBB, dims=1:20, resolution = 0.4)
saveRDS(seurat_project.int.noHBB, "seurat_project_integrated_noBBcells_processed.rds")
Idents(seurat_project.int.noHBB) <- seurat_project.int.noHBB@meta.data$seurat_clusters
seurat_project.int.noHBB.markers <- FindAllMarkers(seurat_project.int.noHBB, min.pct = 0.25, logfc.threshold = 0.25, only.pos = TRUE, random.seed = 777)
write.table(seurat_project.int.noHBB.markers, "seurat_project_int_noHBBcells_allMarkers_onlyUp.txt", sep = "\t", quote = F, col.names = NA)
saveRDS(seurat_project.int.noHBB.markers, "seurat_project_int_noHBBcells_allMarkers_onlyUp.rds")
top10 <- seurat_project.int.noHBB.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(seurat_project.int.noHBB, features = as.vector(top10$gene)) 
height80dpi <- (length(as.vector(top10$gene))*12)*(1/80)+1
ggsave("heatmap_allMarkers_onlyUp.pdf", height= height80dpi, width = 8)

dir.create("Features")
for (p in bigList) {
  l=get(p)
  curdir = paste0("Features/", p)
  dir.create(curdir)
  for (gene in top10$gene) {
    print(gene)
    if (gene %in% rownames(seurat_project.int.noHBB@assays$RNA@counts)) {
      FeaturePlot(seurat_project.int.noHBB, features = gene)
      ggsave(paste0(curdir, "/", gene, ".png"), device = "png")
    } else {
      print("not found")
    }
  }
}

################### QC  & Filtration ####################


########## Merging On Experimental Conditions ###########


############### Dim.Red., UMAP, Clustering ##############

############### FindMarkers ############################# 
# & include singleR (or some derivative for finding celltypes)


############### DE per celltype or cluster ##############



############### Misc. Plots #############################
# see scRNAseqFunctions.R




















