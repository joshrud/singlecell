# These are various useful functions for analysis of single cell data 
# if using Rstudio, the "outline" feature is VERY useful for finding each function#prevents a sourcefile() run of the whole thing

library(ggrepel)
library(Seurat) #using v3.1.1
library(monocle) #2.14.0
library(ggplot2)
library(dplyr)
library(pheatmap)
library(viridis)
library(reshape2)
library(RColorBrewer)
library(cowplot)
library(future)
library(ggthemes)
library(roxygen2)
library(data.table)
library(reshape2)
library(ks) #this is for kde, don't need this for (most) functions
library(jmuOutlier)
library(scales)
library(fgsea)
library(biomaRt)
library(ggrepel)

####### CONSTANTS ######
GSEA_REF_DIR = "/data/annotations/gsea"
LIGHTGRAY = "#EBEBEB"
GRAY <- "#e6e6e6"
PURPLE <- "#8533ff"
GREEN <- "#00cc00"
RED_ORANGE <- "#FF531A"
ENERGETIC_RED <- "#FF004C"
PURPLE_HSIEH <- "#7F559B"
BRIGHT_PURPLE <- "#9445FF"
PIONYR_COLS <- c("#323B94", "#2D9646", "#3F3F3F", "#000000")
  
########################


###### Preprocess #########
HUMAN_GENES <- useMart(biomart="ensembl",
                       host = "https://dec2021.archive.ensembl.org",
                       dataset="hsapiens_gene_ensembl")
MOUSE_GENES <- useMart(biomart="ensembl",
                       host = "https://dec2021.archive.ensembl.org",
                       dataset="mmusculus_gene_ensembl")
MMU_HSA_HOMOLOGY_ <- read.table(file.path("/data/annotations",
                                          "HOM_MouseHumanSequence.rpt.txt"),
                                sep = "\t", header=T)
MMU_HSA_HOMOLOGY_.mouse <- MMU_HSA_HOMOLOGY_[MMU_HSA_HOMOLOGY_$Common.Organism.Name == "mouse, laboratory",] 
MMU_HSA_HOMOLOGY_.human <- MMU_HSA_HOMOLOGY_[MMU_HSA_HOMOLOGY_$Common.Organism.Name == "human",] 
MMU_HSA_HOMOLOGY <-  merge(MMU_HSA_HOMOLOGY_.mouse, MMU_HSA_HOMOLOGY_.human,
                           by="DB.Class.Key", suffixes=c(".mmu", ".hsa"))
###########################


# roxygen2::roxygenize(package.dir = "/Users/jrudolph/Repos/pionyrtx_codebase/scrna")

source("~/Repos/old/singlecell/clusterGSEA.R") #add directory path, R doesn't know where this is

timeFmt <- function() {
  format(Sys.time(), "%y%m%d_%H%M%S")
}


# this doesn't work exactly how I want it, but is decent for now
#   colors get duplicated because doing a 180º offset is tricky:
#   when we do the second color angle vector, we need an offset,
#   the offset, unfortunately, causes duplications anyway
#   can't figure out how to do both color wheels simultaneously 
#   without duplicating a color; maybe I'm missing something critical? 
barplot_cols <- function(n) {
  sub360 <- function(x) {
    if (x > 360) {
      return((x-(360 * (floor(x/360)))))
    } else return(x)
  }
  anglesize = round(360 / n)
  angleoffset = round(anglesize / 2)
  angles1 <- seq(0,360,anglesize)
  angles2 <- seq((0+180),(180+360),anglesize)
  angles2 <- unlist(sapply(angles2, sub360))
  angles2 <- sapply(angles2, `+`, angleoffset)
  cols1 <- sapply(angles1, hcl, l=70, c=100)
  cols2 <- sapply(angles2, hcl, l=70, c=100)
  cols <- cols1 %zip% cols2
  return(cols)
}


#' `%zip`
#' @description Works just like the zip() definition in python, but is used like so: x %zip% y 
#' @param x A vector 
#' @param y another vector
#'t
#' @return An interleaved vector of x and y
#'
#' @examples c(1,2,3) %zip% c(4,5,6) = c(1,4,2,5,3,6)
`%zip%` <- function(x,y) {
  if (length(x) != length(y)) {
    stop("unequal lengths of vectors to zip...")
  }
  final <- c()
  for (i in 1:min(length(x), length(y))) {
    final <- c(final, x[1], y[1])
    x <- x[-1]
    y <- y[-1]
  }
  return(final)
}

#' convertHumanGeneList
#' @author https://www.r-bloggers.com/2016/10/converting-mouse-to-human-gene-names-with-biomart-package/
#' @author https://www.biostars.org/p/315520/
#' @description For converting gene lists from human to mouse
#' @param x Human gene list
#'
#' @return a mouse gene list
convertHumanGeneList <- function(x){
  genesV2 = getLDS(attributes = c("hgnc_symbol"), 
                   filters = "hgnc_symbol", 
                   values = ANNOT_HUMAN$`Gene name` , 
                   mart = HUMAN_GENES, 
                   attributesL = c("mgi_symbol"), 
                   martL = MOUSE_GENES, 
                   uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  return(humanx)
}

#' convertMouseGeneList
#' @note DEPRECATION: this function is now deprectated in favor of using convertMouseRanks2Hum()
#' 
#' 
#' 
#' 
#' @author https://www.r-bloggers.com/2016/10/converting-mouse-to-human-gene-names-with-biomart-package/
#' @author https://www.biostars.org/p/315520/
#' @description For converting gene lists from mouse to human
#' @param x Mouse gene list
#'
#' @return a huan gene list
# convertMouseGeneList <- function(x){
#   genesV2 = getLDS(attributes = c("mgi_symbol"), 
#                    filters = "mgi_symbol", 
#                    values = x , 
#                    mart = MOUSE_GENES, 
#                    attributesL = c("hgnc_symbol"), 
#                    martL = HUMAN_GENES, 
#                    uniqueRows=T)
#   humanx <- unique(genesV2[, 2])
#   return(humanx)
# }


#' hsa2mmu
#' @description Converts human genes to mouse genes for GSEA analysis 
#' @author JR
#' @param x the array of human gene names 
#' @param opt the option for choosing the matching mouse id; can be "first" or "all" 
#' @return an array of mouse gene names matching the order of the original human gene names
hsa2mmu <- function(x, opt = "all") {
  y <- c()
  count=0
  for (i in seq_along(x)) {
    temp <- MMU_HSA_HOMOLOGY$Symbol.mmu[which(MMU_HSA_HOMOLOGY$Symbol.hsa == x[i])]
    if (length(temp) > 1) {
      message(paste0("Symbol ", x[i], " has ", length(temp), " features: ", paste(temp, collapse=", ")))
      count <- count+1
      if (opt=="first") {
        y <- c(y, temp[1]) # include only the first matching one, can miss some
        
        ### should be improved so that each match has a score based on how mismatched it is, least mismatched is chosen instead of first ....
        
      } else if (opt=="all") {
        for (j in seq_along(temp)) { # include all of the matching ones
          y <- c(y, temp[j])
        }
      } else {
        message("opt needs to be either \"first\" or \"all\" ")
        return(character(0)) 
      }
    } else {
      y <- c(y, temp)
    }
  }
  message("there were ", count, " genes in human list that had more than 1 mouse gene")
  message("there are ", length(y[duplicated(y)]), " duplicated genes in final list, removing..." )
  y <- y[!duplicated(y)]
  message("there are now ", length(y), " genes in final list." )
  return(y)
}



# input is a rank file as if it were going into fgsea, but\
#  instead we'll output a similar rank file with human names
convertMouseRanks2Hum <- function(x) {
  # some mmu genes have multiple homologous genes in human,\
  #  Cxcl3, for example, is homologous to CXCL1, CXCL2,and CXCL3
  # solution: add ALL to the new list and return that 
  y <- c()
  if (exists("MMU_HSA_HOMOLOGY", sys.frame())) {
    for (i in seq_along(x)) {
      x.human <- MMU_HSA_HOMOLOGY$Symbol.hsa[which(MMU_HSA_HOMOLOGY$Symbol.mmu == names(x)[i])]
      if (identical(x.human, character(0))) {
        next
      }
      if (length(x.human) > 1) {
        for (j in seq_along(x.human)) {
          if (length(y) == 0) {
            y <- x[i]
            names(y) <- x.human[j]
          } else {
            y <- c(y, x[i])
            names(y)[length(y)] <- x.human[j]
          }
        }
      } else {
        if (length(y) == 0) {
          y <- x[i]
          names(y) <- x.human
        } else {
          y <- c(y, x[i])
          names(y)[length(y)] <- x.human
        }
      }
    }
  } else {
    stop("must load MMU_HSA_HOMOLOGY...")
  }
  
  # now get unique names
  y.u <- unique(names(y))
  y.final <- c()
  for (i in seq_along(y.u)) {
    if (length(y[names(y) == y.u[i]]) > 1) {
      if (length(y.final) == 0) {
        y.final <- median(y[which(names(y) == y.u[i])])
        names(y.final)[i] <- y.u[i]
      } else {
        y.final <- c(y.final, median(y[which(names(y) == y.u[i])]))
        names(y.final)[i] <- y.u[i]
      }
    } else {
      if (length(y.final) == 0) {
        y.final <- y[which(names(y) == y.u[i])]
      } else {
        y.final <- c(y.final, y[which(names(y) == y.u[i])])
      }
    }
  }
  return(y.final)
}

#' Assign mitochondrial percent, hbb percent, ribo percent
#'
#' @param seurat_object The seurat object to assign meta data to 
#' @param mt.genes Inputted mitochondrial features, can be left blank
#' @param ribo.genes Inputted ribosomal features, can be left blank
#' @param hbb.genes Inputted hemoglobin features, can be left blank
#'
#' @return The updated seurat object
#' @examples iso <- add_meta_percents(iso)
add_meta_percents <- function(seurat_object, mt.genes = NULL, ribo.genes = NULL, hbb.genes = NULL) {
  allgenes = rownames(seurat_object@assays$RNA@counts)
  if (is.null(mt.genes)) {
    mt.genes <- allgenes[grep("^[Mm][Tt]-", allgenes)]
  }
  if (is.null(ribo.genes)) {
    ribo.genes <- allgenes[grep("^[Rr][Pp][Ss]|[Rr][Pp][Ll]", allgenes)]
  }
  if (is.null(hbb.genes)) {
    hbb.genes <- allgenes[grep("^H[Bb][BbAaQq]", allgenes)]
  }
  print(paste0("mt.genes: ", paste(mt.genes, collapse = " ")))
  print(paste0("ribo.genes: ", paste(ribo.genes, collapse = " ")))
  print(paste0("hbb.genes: ", paste(hbb.genes, collapse = " ")))
  
  #assign and return
  seurat_object[["percent.mt"]] <- PercentageFeatureSet(object = seurat_object, features = mt.genes)
  seurat_object[["percent.ribo"]] <- PercentageFeatureSet(object = seurat_object, features = ribo.genes)
  seurat_object[["percent.hbb"]] <- PercentageFeatureSet(object = seurat_object, features = hbb.genes)
  return(seurat_object)
}

pp_plots <- function(seurat_object, dir=".") {
  seurat_object[["percent.mt"]] <- PercentageFeatureSet(object = seurat_object, pattern = "^[Mm][Tt]-")
  seurat_object[["percent.ribo"]] <- PercentageFeatureSet(object = seurat_object, pattern = "^[Rr][Pp][Ss]|[Rr][Pp][Ll]")
  seurat_object[["percent.hbb"]] <- PercentageFeatureSet(object = seurat_object, pattern = "^H[Bb][BbAaQq]")
  featureinfo <- seurat_object@meta.data[,c("nCount_RNA", "nFeature_RNA", "orig.ident", "percent.mt", "percent.ribo", "percent.hbb")]
  VlnPlot(seurat_object, features=c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo", "percent.hbb"), pt.size = 0)
  ggsave(paste0(dir, "/", seurat_object@project.name, "_violin.pdf"), width=10,height=10)
  FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 1)
  ggsave(paste0(dir, "/", seurat_object@project.name, "_featuresCounts.pdf"), width=6,height=4)
}

pp <- function(seurat_object, HEMOGENES=NULL, subset=TRUE, SCT = FALSE, VAR.FEATURES=3000, MAX_FEATURES_PER_CELL=2500, MIN_FEATURES_PER_CELL=100, MAX_PERCENT_MT=5, DIMS=40) {
  genes <- rownames(seurat_object@assays$RNA@counts)
  genes.mt <- genes[grep("^M[Tt]-", genes)]
  print(paste0("found ", length(genes.mt), " genes: ", paste(genes.mt, collapse = " ")))
  seurat_object[["percent.mt"]] <- PercentageFeatureSet(object = seurat_object, pattern = "^M[Tt]-")
  
  genes.ribo <- genes[grep("^R[Pp][Ss]|R[Pp][Ll]", genes)]
  print(paste0("found ", length(genes.ribo), " genes: ", paste(genes.ribo, collapse = " ")))
  seurat_object[["percent.ribo"]] <- PercentageFeatureSet(object = seurat_object, pattern = "^R[Pp][Ss]|R[Pp][Ll]")
  
  if (is.null(HEMOGENES)) {
    genes.hemo <- genes[grep("^H[Bb][Bb]|^H[Bb][Aa][12]", genes)]
  } else genes.hemo <- HEMOGENES
  print(paste0("found ", length(genes.hemo), " genes: ", paste(genes.hemo, collapse = " ")))
  seurat_object[["percent.hbb"]] <- PercentageFeatureSet(object = seurat_object, features = genes.hemo)
  featureinfo <- seurat_object@meta.data[,c("nCount_RNA", "nFeature_RNA", "orig.ident", "percent.mt", "percent.ribo", "percent.hbb")]
  #subset on filtered cells and continue
  if (subset) {
    seurat_object <- subset(x = seurat_object, subset = nFeature_RNA > MIN_FEATURES_PER_CELL & nFeature_RNA < MAX_FEATURES_PER_CELL & percent.mt < MAX_PERCENT_MT)
  }
  if (SCT) {
    seurat_object <- SC
    Transform(seurat_object, vars.to.regress = "percent.mt") #use this to replace scaledata if it's not good enough
  } else {
    seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 1e4)
    seurat_object <- FindVariableFeatures(object = seurat_object, selection.method = 'vst', nfeatures = VAR.FEATURES)
    seurat_object <- ScaleData(object = seurat_object, vars.to.regress = c('percent.mt') )
  }
  seurat_object <- RunPCA(object = seurat_object, features = VariableFeatures(object = seurat_object))
  seurat_object <- RunUMAP(object = seurat_object, dims=1:DIMS) #VariableFeatures(seurat_object) makes it error -> need dims for some reason 
  seurat_object <- FindNeighbors(object = seurat_object, dims=1:DIMS)
  seurat_object <- FindClusters(object = seurat_object, dims=1:DIMS, resolution = 0.4)
  message("saving...")
  saveRDS(seurat_object, paste0(seurat_object@project.name, "_pp.rds"))
  return(seurat_object)
}

allUMAPs <- function(seurat_object) {
  UMAPPlot(seurat_object, group.by="seurat_clusters", label=TRUE)
  ggsave("umap_byCluster.png", device = "png", units = "in", width = 7, height = 6)
  UMAPPlot(seurat_object, group.by="condition")
  ggsave("umap_byCondition.png", device = "png", units = "in", width = 7, height = 6)
  UMAPPlot(seurat_object, group.by="orig.ident")
  ggsave("umap_bySample.png", device = "png", units = "in", width = 7, height = 6)
}

clusterDE <- function(seurat_object, name) {
  Idents(seurat_object) <- seurat_object@meta.data$seurat_clusters
  seurat_object.markers <- FindAllMarkers(seurat_object, min.pct = 0.25, logfc.threshold = 0.25, random.seed = 777)
  write.table(seurat_object.markers, paste0("markers_", name, ".txt"), sep = "\t", quote = F, col.names = NA)
  top10 <- seurat_object.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  DoHeatmap(seurat_object, features = as.vector(top10$gene)) 
  height80dpi <- (length(as.vector(top10$gene))*12)*(1/80)+1
  ggsave(paste0("heatmap_", name, ".png"), height= height80dpi, width = 8)
  
  ### print a umap with top n markers annotated for each cluster
  
}

toClipboard <- function(x) {
  clip <- pipe("pbcopy", "w")
  write.table(x, file=clip)
  close(clip)
  print("saved to clipboard")
}

zipTable <- function(x, separator = ":") {
  zipped <- c()
  for (i in 1:length(x)) {
    zipped <- c(zipped, paste0(names(x)[i], separator, " ", (round(x[i],3) * 100))) 
  }
  return(zipped)
}

#basically python's zip function ex: zip((x1,x2), (y1,y2)) -> [(x1,y1), (x2,y2)]
zip <- function(x, separator = "") {
  zipped <- c()
  for (i in 1:length(x)) {
    zipped <- c(zipped, paste0(names(x)[i], separator, " ", x[i]))
  }
  return(zipped)
}

wrapper <- function(x, ...) 
{
  paste(strwrap(x, ...), collapse = "\n")
}

#used to sort cluster numbers, especially when they contain whole-digit floats (like 1.0)
sortCharArrayByNum <- function(x) {
  if (!is.numeric(x)) {
    return(x)
  }
  positions <- length(x)
  for (i in 1:(positions-1)) {
    for (j in 1:(positions-i)) {
      if (as.numeric(as.character(x[j]))  > as.numeric(as.character(x[j+1]))) { #flip
        temp <- x[j]
        x[j] <- x[j+1]
        x[j+1] <- temp
      } 
    }  
  }
  return(x)
}

#makes a table of the number of cells per celltype_col, split by sample_col, and makes a barplot from the data (this part should probably just be taken out)
#requires ggplot
getCellNumbers <- function(seurat_object, celltype_col = "seurat_clusters", sample_col = "orig.ident", write_output = F, outdir = NULL) {
  
  sampleNames <- unique(seurat_object@meta.data[,sample_col])
  meta <- seurat_object@meta.data
  celltypes <- unique(meta[,celltype_col])
  if (!is.na(as.numeric(celltypes[1]))) {
    celltypes <- sortCharArrayByNum(celltypes)
  } else print("chosen a non-default column, non-numeric")
  
  clh.incelltype <- c()
  cur.vec <- c()
  for (c in celltypes) {
    cur.vec <- c()
    for (i in 1:length(sampleNames)) {
      if (i==1) {
        cur.vec <- c(length(meta[which(meta[,sample_col] == sampleNames[i] & meta[,celltype_col] == c), celltype_col]))
      }
      cur.vec <- c(cur.vec,  length(meta[which(meta[,sample_col] == sampleNames[i] & meta[,celltype_col] == c), celltype_col]))
    }
    clh.incelltype <- rbind(clh.incelltype, cur.vec)
  }
  clh.incelltype <- clh.incelltype[,-1]
  colnames(clh.incelltype) <- sampleNames
  rownames(clh.incelltype) <- celltypes
  clh.incelltype.melt <- melt(clh.incelltype)
  colnames(clh.incelltype.melt) <- c("celltype", "sample", "number_cells")
  clh.incelltype.melt$celltype <- as.factor(clh.incelltype.melt$celltype)
  if (write_output) {
    if (is.null(outdir)) {
      write.table(clh.incelltype, paste0("cellNumbers_", seurat_object@project.name, ".txt"), quote = F, col.names = NA, sep = "\t")
    } else {
      write.table(clh.incelltype, paste0(outdir, "/cellNumbers_", seurat_object@project.name, ".txt"), quote = F, col.names = NA, sep = "\t")
    }
    
    ## barplot1 
    ggplot(data=clh.incelltype.melt,aes(x=celltype, y=number_cells, group=sample, color=sample)) + geom_bar(stat="identity", position = "dodge", width=0.7) +
      labs(y="Number of cells", x="celltype") + 
      guides(fill=guide_legend(title="Sample")) +
      scale_y_continuous(expand = c(0,0)) + #removes the annoying space between x axis label and the bottom of plot
      theme(axis.text.x = element_text(angle = 45, hjust=1,vjust=1), panel.background = element_rect(fill="#E5E5E5"), 
            panel.grid.major = element_line(size=0.5,color="white"), panel.grid.minor = element_line(size=0.2,color="white"),
            panel.border = element_blank() )
    ggsave("cells_per_sample_per_celltype.pdf")
    
    ### barplot2
    ggplot(data=clh.incelltype.melt,aes(x=sample, y=number_cells, group=celltype, color=celltype)) + geom_bar(stat="identity", position = "dodge", width=0.7) +
      labs(y="Number of cells", x="celltype") + 
      guides(fill=guide_legend(title="Sample")) +
      scale_y_continuous(expand = c(0,0), limits=c(0,3000)) + #removes the annoying space between x axis label and the bottom of plot
      theme(axis.text.x = element_text(angle = 45, hjust=1,vjust=1), panel.background = element_rect(fill="#E5E5E5"), 
            panel.grid.major = element_line(size=0.5,color="white"), panel.grid.minor = element_line(size=0.2,color="white"),
            panel.border = element_blank() )
    ggsave("cells_per_celltype_per_sample.pdf")
  }

  return(clh.incelltype)
  print("done")
}


#makes umaps with pie charts (needs to be tested, but I know it works for me)
#requires ggplot, scatterpie
# to add: check variability of location of cells for a cluster, if it's past a threshold, choose the largest "chunk" of cells for center
#' umap_withPies
#' @description Prints a ggplot object and saves in getwd() with table of pie-chart data
#' @param seurat_object The seurat object with umap already calculated
#' @param pie_cond The column in meta.data that will be used for proportioning pies
#' @param pie_loc The column in meta.data that will be used for positioning pies 
#' @param normalize Whether to normalize (logical)
#'
#' @return NULL
umap_withPies <- function(seurat_object, pie_cond = "orig.ident", pie_loc = "seurat_clusters", normalize = TRUE) {
  
  if (normalize) {
    message("WARNING: normalize function normalizes by proportion of PIE_LOC, not PIE_COND")
  }
  clusteringRes <- colnames(seurat_object@meta.data)[grep("snn_res", colnames(seurat_object@meta.data))]
  clusteringRes <- unlist(strsplit(clusteringRes, "_"))[length(unlist(strsplit(clusteringRes, "_")))]
  totalCells <- nrow(seurat_object@meta.data)
  samplePercentages <- table(seurat_object@meta.data[,pie_cond]) / totalCells
  clusters <- sort(as.character(unique(seurat_object@meta.data[,pie_loc])))
  
  if (length(intersect(clusters, unique(seurat_object@meta.data[,pie_loc]))) != length(clusters)) {
    print("pie_loc not the same")
    print(clusters)
    print(unique(seurat_object@meta.data[,pie_loc]))
    return(NULL)
  }
  
  umap <- seurat_object@reductions$umap@cell.embeddings #SOlist[[x]]@dr$tsne@cell.embeddings
  umap.withmeta <- merge(umap, seurat_object@meta.data, by='row.names')
  rownames(umap.withmeta) <- umap.withmeta[,1]
  umap.withmeta <- umap.withmeta[,-1]
  umap.withmeta$pie.x <- NA
  umap.withmeta$pie.y <- NA
  umap.withmeta$pie.xmed <- NA
  umap.withmeta$pie.ymed <- NA
  for (i in clusters) { #take average of x,y for each cluster for location of pie chart
    cur.x <- umap.withmeta[which(umap.withmeta[,pie_loc] == i),'UMAP_1']
    cur.y <- umap.withmeta[which(umap.withmeta[,pie_loc] == i),'UMAP_2']
    cur.x.mean <- mean(cur.x)
    cur.y.mean <- mean(cur.y)
    cur.x.median <- median(cur.x)
    cur.y.median <- median(cur.y)
    umap.withmeta$pie.x[which(umap.withmeta[,pie_loc] == i)] <- cur.x.mean #mean
    umap.withmeta$pie.y[which(umap.withmeta[,pie_loc] == i)] <- cur.y.mean #mean
    umap.withmeta$pie.xmed[which(umap.withmeta[,pie_loc] == i)] <- cur.x.median #median
    umap.withmeta$pie.ymed[which(umap.withmeta[,pie_loc] == i)] <- cur.y.median #median
  }
  pies <- data.frame("cluster"=clusters, "pie.x"=rep(0,length(clusters)), "pie.y"=rep(0,length(clusters)),"pie.xmed"=rep(0,length(clusters)), "pie.ymed"=rep(0,length(clusters)))
  for (c in unique(seurat_object@meta.data[,pie_cond]) ) {
    pies[,c] <- rep(0, length(clusters))
  }
  
  umap.withmeta[,pie_cond] <- as.character(umap.withmeta[,pie_cond])
  for (j in clusters) { #pie data: number of cells from each clh in each cluster
    for (i in unique(umap.withmeta[,pie_cond])) {
      pies[which(pies$cluster==j),i] <- length(umap.withmeta[,pie_cond][which(umap.withmeta[,pie_cond]==i & umap.withmeta[,pie_loc]==j)])
    }
  }
  
  for (i in 1:length(clusters)) { #pie posiiton: median or mean
    pies$pie.x[i] <- unique(umap.withmeta$pie.x[which(umap.withmeta[,pie_loc]==clusters[i])])
    pies$pie.y[i] <- unique(umap.withmeta$pie.y[which(umap.withmeta[,pie_loc]==clusters[i])])
    pies$pie.xmed[i] <- unique(umap.withmeta$pie.xmed[which(umap.withmeta[,pie_loc]==clusters[i])])
    pies$pie.ymed[i] <- unique(umap.withmeta$pie.ymed[which(umap.withmeta[,pie_loc]==clusters[i])])
  }
  pies$totals <- rowSums(pies[,unique(seurat_object@meta.data[,pie_cond])]) 
  pies.norm <- pies[,unique(seurat_object@meta.data[,pie_cond])]
  print(pies$totals)
  print(pies.norm)
  cond.coefs <- sum(pies.norm)/colSums(pies.norm)
  pies.norm.f <- t(apply(pies.norm , 1, function(x) return(cond.coefs * x)))
  colnames(pies.norm.f) <- paste0(colnames(pies.norm.f), ".norm")
  print(pies.norm.f)
  pies <- cbind(pies, pies.norm.f)
  pies.norm.percentages <- pies.norm.f / rowSums(pies.norm.f)
  colnames(pies.norm.percentages) <- paste0(colnames(pies.norm.f), ".pct")
  pies <- cbind(pies, pies.norm.percentages)
  if (normalize) {
    #finally plot the data with the pie charts
    ggplot(data=umap.withmeta, aes_string(x="UMAP_1", y="UMAP_2", group=pie_loc, color=pie_loc)) + geom_point(size=0.3) +
      geom_scatterpie(data=pies, aes(x=pie.xmed, y=pie.ymed), cols=paste0(unique(seurat_object@meta.data[,pie_cond]), ".norm"), pie_scale = 4) + 
      theme_classic() + 
      theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.background = element_rect(fill = "white", color="black"), panel.grid = element_blank()) +
      guides(color = guide_legend(override.aes = list(size = 10))) + 
      labs(title = paste0(seurat_object@project.name, " UMAP with ", clusteringRes), 
           subtitle = wrapper(paste0("total cells: ", nrow(seurat_object@meta.data), " and percents: ", paste(zipTable(samplePercentages), collapse="% "), "%" ), width=140) ) 
    #, r=log2(totals)/log2(sum(totals)))
  } else {
    #finally plot the data with the pie charts
    ggplot(data=umap.withmeta, aes_string(x="UMAP_1", y="UMAP_2", group=pie_loc, color=pie_loc)) + geom_point(size=0.3) +
      geom_scatterpie(data=pies, aes(x=pie.xmed, y=pie.ymed), cols=unique(seurat_object@meta.data[,pie_cond]), pie_scale = 4) + 
      theme_classic() + #, r=log2(totals)/log2(sum(totals)))
      theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.background = element_rect(fill = "white", color="black"), panel.grid = element_blank()) +
      guides(color = guide_legend(override.aes = list(size = 10))) + 
      labs(title = paste0(seurat_object@project.name, " UMAP with ", clusteringRes), 
           subtitle = wrapper(paste0("total cells: ", nrow(seurat_object@meta.data), " and percents: ", paste(zipTable(samplePercentages), collapse="% "), "%"), width=140) ) 
  }
  ggsave("umap_with_pie.pdf")
  write.table(pies, "umap_pies_info.txt", sep = "\t", quote = F, col.names = NA)
  print(paste0("figures and table are in: ", getwd()))
  return(NULL)
} 

#found online, useful for adding stats to a dataframe of celltype proportions (like median bar, SEM, etc..)
##http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/#Helper%20functions
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}


normDataWithin <- function(data=NULL, idvar, measurevar, betweenvars=NULL,
                           na.rm=FALSE, .drop=TRUE) {
  library(plyr)
  
  # Measure var on left, idvar + between vars on right of formula.
  data.subjMean <- ddply(data, c(idvar, betweenvars), .drop=.drop,
                         .fun = function(xx, col, na.rm) {
                           c(subjMean = mean(xx[,col], na.rm=na.rm))
                         },
                         measurevar,
                         na.rm
  )
  
  # Put the subject means with original data
  data <- merge(data, data.subjMean)
  
  # Get the normalized data in a new column
  measureNormedVar <- paste(measurevar, "_norm", sep="")
  data[,measureNormedVar] <- data[,measurevar] - data[,"subjMean"] +
    mean(data[,measurevar], na.rm=na.rm)
  
  # Remove this subject mean column
  data$subjMean <- NULL
  
  return(data)
}

#http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/#Helper%20functions
summarySEwithin <- function(data=NULL, measurevar, betweenvars=NULL, withinvars=NULL,
                            idvar=NULL, na.rm=FALSE, conf.interval=.95, .drop=TRUE) {
  
  # Ensure that the betweenvars and withinvars are factors
  factorvars <- vapply(data[, c(betweenvars, withinvars), drop=FALSE],
                       FUN=is.factor, FUN.VALUE=logical(1))
  
  if (!all(factorvars)) {
    nonfactorvars <- names(factorvars)[!factorvars]
    message("Automatically converting the following non-factors to factors: ",
            paste(nonfactorvars, collapse = ", "))
    data[nonfactorvars] <- lapply(data[nonfactorvars], factor)
  }
  
  # Get the means from the un-normed data
  datac <- summarySE(data, measurevar, groupvars=c(betweenvars, withinvars),
                     na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)
  
  # Drop all the unused columns (these will be calculated with normed data)
  datac$sd <- NULL
  datac$se <- NULL
  datac$ci <- NULL
  
  # Norm each subject's data
  ndata <- normDataWithin(data, idvar, measurevar, betweenvars, na.rm, .drop=.drop)
  
  # This is the name of the new column
  measurevar_n <- paste(measurevar, "_norm", sep="")
  
  # Collapse the normed data - now we can treat between and within vars the same
  ndatac <- summarySE(ndata, measurevar_n, groupvars=c(betweenvars, withinvars),
                      na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)
  
  # Apply correction from Morey (2008) to the standard error and confidence interval
  #  Get the product of the number of conditions of within-S variables
  nWithinGroups    <- prod(vapply(ndatac[,withinvars, drop=FALSE], FUN=nlevels,
                                  FUN.VALUE=numeric(1)))
  correctionFactor <- sqrt( nWithinGroups / (nWithinGroups-1) )
  
  # Apply the correction factor
  ndatac$sd <- ndatac$sd * correctionFactor
  ndatac$se <- ndatac$se * correctionFactor
  ndatac$ci <- ndatac$ci * correctionFactor
  
  # Combine the un-normed means with the normed results
  merge(datac, ndatac)
}



#adds statistics to a barplot
addStats <- function(x, d) {
  install.packages("ggsignif")
  
  # Or for the latest development version
  devtools::install_github("const-ae/ggsignif")
  
  
  
  x <- x + 
  geom_signif(comparisons = list(c("compact", "midsize"), c("minivan", "suv")),
              map_signif_level = TRUE, textsize=6) 
    
    
}


#' Performs Differential Expression testing for all idents of a seurat object, \n
#' outputs volcano plots, de marker table, and go analysis for all comparisons \n
#' in each cluster
#'
#' @param seurat_object The seurat object.
#' @param output The directory we want to output DE results directory to 
#' @param comps The comparisons to perform in the form "cond1vscond2", comparing cond1 vs. cond2; both should be in "conditions" column in meta.data
#' @param coi The 'column of interest' in the meta.data that we want to split object by to perform DE testing, ex: "seurat_clusters" 
#' @param with_labels Whether we want labels in the volcano plot 
#' @param padj_thresh Adjusted pvalue threshold
#' @param org_db The database used for performing GO analysis, either \
#' "org.Mm.eg.db" or "org.Hs.eg.db", or "auto" for detecting from gene names which one s
#' @return DE results between comps in each coi, with volcanos and GO results for each
#'
easyDE <- function(seurat_object, 
                   output, 
                   comps, 
                   coi="seurat_clusters", 
                   with_labels = T,
                   padj_thresh=0.05, 
                   orgdb = "auto") {
  if (orgdb == "auto") {
    # check for beta-actin, a HK gene 
    if ("ACTB" %in% rownames(seurat_object)) { # Actb is the mouse gene
      orgdb <- "org.Hs.eg.db" 
    } else {
      orgdb <- "org.Mm.eg.db" 
    }
  } 
  
  # initialize path.df
  path.df <- c()
  
  # notify 
  timestamp()
  print(paste0("Detected ", orgdb, " as database."))
  
  req.packages <- c("ggthemes",
                    "ggplot2",
                    orgdb)
  pack.mat <- sapply(req.packages, FUN = function(x) {
             if (paste0("package:",x) %in% search()) {
               return(TRUE)
             } else return(FALSE)
           } )
  if (any(!pack.mat)) {
    print(paste0("packages: ", paste0(req.packages[!pack.mat], collapse = " "), 
                 " need to be loaded or installed..."))
    return(NULL)
  } else {
    print("good to go :)")
  }
  
  if (!exists("scRNA_pathways", where=sys.frame())) {
    print("GO function not loaded, breaking...")
    return(NULL)
  }
  
  clusters = unique(seurat_object@meta.data[,coi])
  for (i in seq_along(clusters)) {
    clusterfolder <- file.path(output, paste0("cluster_", clusters[i], "_DE"))
    dir.create(clusterfolder)
    
    cluster.cells <- rownames(seurat_object@meta.data[which(seurat_object@meta.data[,coi] == clusters[i]),])
    cluster.set <- subset(seurat_object, cells = cluster.cells) #only cells in this cluster
    for (comp in comps) { #
      
      timestamp()
      print(paste0(" working on ", comp, " in cluster ", 
                   clusters[i], " of object ", seurat_object@project.name))
      compfolder <- file.path(clusterfolder, comp)
      dir.create(compfolder)
      splitted <- unlist(strsplit(comp, "vs"))
      cond1 <- splitted[1]
      cond2 <- splitted[2]
      
      cond1.cells <- rownames(cluster.set@meta.data[which(cluster.set@meta.data$condition == cond1),])
      cond2.cells <- rownames(cluster.set@meta.data[which(cluster.set@meta.data$condition == cond2),])
      if (length(cond1.cells) <= 3 | length(cond2.cells) <=3 ) {
        print(paste0("comp ", comp, " in cluster ", clusters[i],
                     " has too few cells. Skipping"))
        next
      }
      
      cluster.set <- SetIdent(cluster.set, cells = cond1.cells, value=cond1)
      cluster.set <- SetIdent(cluster.set, cells = cond2.cells, value=cond2)
      notEnoughGenes <- tryCatch({
        cur.de <- FindMarkers(cluster.set, ident.1 = cond1, ident.2 = cond2, 
                              min.pct = 0.1, logfc.threshold = 0.1) 
      }, error = function(err) {
        notEnoughGenes <- TRUE
        print(paste("Caught an error: ", err, " ...continuing on...", sep = ""))
        return(TRUE)
      })
      if (typeof(notEnoughGenes)!="logical") { # this means an error wasn't thrown in the above tryCatch 
        notEnoughGenes <- FALSE
      }
      if (notEnoughGenes) next #this is where not enough genes is actually checked for errors
      
      # write both of the tables
      cur.de <- cur.de[order(cur.de$avg_log2FC,decreasing=T),]
      write.table(cur.de, file.path(compfolder, paste0("DEmarkers_", 
                                                       clusters[i], "_",cond1, 
                                    "_all", ".tsv")), 
                  quote = F, sep = "\t", col.names = NA)
      cur.de.sig <- cur.de[cur.de$p_val_adj < padj_thresh,]
      write.table(cur.de.sig, file.path(compfolder, paste0("DEmarkers_", 
                                                      clusters[i], "_",cond1, 
                                                       "_sig", ".tsv")), 
                  quote = F, sep = "\t", col.names = NA)
      
      # pathway analysis 
      go_folder <- file.path(compfolder, "GO")
      hall_folder <- file.path(compfolder, "Hallmark")
      dir.create(go_folder)
      dir.create(hall_folder)
        
      # submit the marker matrix to the function 
      gsea_res.go <- scRNA_pathways(cur.de,  
                                  outfolder =  go_folder, 
                                  type = "go",
                                  paths = NULL,
                                  pathname = paste0(clusters[i],"_",comp, "_GO_paths"),
                                  orgdb = orgdb)
      gsea_res.h <- scRNA_pathways(cur.de, 
                                 outfolder =  hall_folder, 
                                 type = "hallmark", 
                                 paths = NULL,
                                 pathname = paste0(clusters[i],"_",comp, "_Hall_paths"),
                                 orgdb = orgdb)
      if (nrow(gsea_res.h) > 0) {
        temp.h <- cbind.data.frame("pathway"=gsea_res.h$pathway, 
                                    "pval"=gsea_res.h$pval,
                                    "padj"=gsea_res.h$padj, 
                                    "NES"=gsea_res.h$NES,
                                    "celltype"=clusters[i],
                                    "comp"=comp)
        if (is.null(path.df) ) {
          path.df <- temp.h
        } else {
          path.df <- rbind(path.df, temp.h)
        }
      }
      
      # volcano plots 
      if (with_labels) {
        easyDE.volcano(cur.de = cur.de, clust = clusters[i], 
                       comp = comp, with_labels = T, #with labels
                       compfolder = compfolder)
      } else {
        easyDE.volcano(cur.de = cur.de, clust = clusters[i], 
                       comp = comp, with_labels = F, #without labels
                       compfolder = compfolder)
      }
    }  
  }
  # write all of the combined hallmark paths so that we won't \
  #  have to collect them later down the road 
  write.table(path.df, file=file.path(output, "all_hallmark_paths.tsv"),
              quote=F,sep = "\t", col.names = NA)
}


#' Volcano Plot code for easyDE
#' @description Prints the volcano plot for easyDE so that easyDE is easier to follow
#' @param cur.de Requires same columns to a markers text file, with "logPadjf" = -log(pvalue) and gene names in rownames
#' @param clust current cluster from the easyDE function
#' @param comp current comparison string
#' @param with_labels Whether we want labels in volcano
#' @param compfolder The folder we're outputting to, already written in easyDE
#' @return NULL 
easyDE.volcano <- function(cur.de, 
                           clust = NULL, 
                           comp = NULL, 
                           with_labels = TRUE, 
                           compfolder) {
  
  if (is.null(comp)) {
    comp <- "placeholder1vsplaceholder2"
  }
  splitted <- unlist(strsplit(comp, "vs"))
  cond1 <- splitted[1]
  cond2 <- splitted[2]
  
  # volcano plot ----
  cur.de$logPadjf <- -log10(cur.de$p_val_adj + 1e-300)
  volcdat <- cur.de[,grep("avg_log2FC|logPadjf", colnames(cur.de))]
  cur.de$labels <- NA
  if (with_labels) {
    
    #new strategy: label top 5 genes by significance, top 5 by FC, and lowest 5 by FC
    top5.p <- head(rownames(cur.de[order(cur.de$logPadjf, decreasing = T),]), n=5)
    top5.l2fc.up <- head(rownames(cur.de[order(cur.de$avg_log2FC, decreasing = T),]), n=5)
    top5.l2fc.down <- head(rownames(cur.de[order(cur.de$avg_log2FC, decreasing = F),]), n=5)
    
    if (nrow(cur.de) < 20) {
      cur.de$labels <- rownames(cur.de)
    } else {
      cur.de$labels[which(rownames(cur.de) %in% c(top5.p, 
                                                  top5.l2fc.up, 
                                                  top5.l2fc.down))] <- 
        rownames(cur.de[which(rownames(cur.de) %in% c(top5.p, 
                                                      top5.l2fc.up, 
                                                      top5.l2fc.down)),])
    }
  }
  
  cur.de$coloring <- NA
  cur.de$coloring[which(cur.de$p_val_adj < 0.05 & 
                          cur.de$avg_log2FC > 1)] <- "de"
  cur.de$coloring[which(is.na(cur.de$coloring))] <- "not_de"
  de.colors <- c("de"="red", 
                 "not_de"="lightgray")
  if (with_labels) {
    cur.de$coloring[which(!is.na(cur.de$labels))] <- "labeled"
    de.colors <- c(de.colors,
                   "labeled"="blue")
  } 
  
  cur.de$cumulative_diff_btw_clusters <- abs(cur.de$pct.2 - cur.de$pct.1)
  r <- max(abs(cur.de$avg_log2FC)) + 1

  # make the ggplot 
  volcano <- ggplot2::ggplot(data=cur.de, aes(x = avg_log2FC, y = logPadjf, 
                                              size=cumulative_diff_btw_clusters)) + 
    geom_point(aes(alpha=0.5, color=coloring)) + 
    geom_text_repel(aes(label = labels), point.padding = 0.25, size = 3) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", 
               size=0.5, alpha = 0.4) +
    geom_vline(xintercept=c(-1,1), linetype="dashed", 
               size=0.5, alpha = 0.4) +
    scale_color_manual(name="Diff. Expressed", values=de.colors) +
    scale_x_continuous(limits = c(-r, r)) +
    xlab(paste0("log2 Fold change (", cond1, "/", cond2, ")")) + 
    ylab("-log10 adjusted p-value") +
    ggtitle(paste0(comp)) +
    guides(alpha=F) + 
    theme_bw() + 
    theme(legend.key.size = unit(2, 'lines'),
          text = element_text(size=10),
          legend.text=element_text(size=8), legend.key.height = unit(2, 'lines'),
          legend.position = "right")
  
  if (with_labels) {
    pdf(file.path(compfolder, 
                  paste0("volcano_", clust, "_", comp, "_withLabels.pdf")))
  } else {
    pdf(file.path(compfolder, 
                  paste0("volcano_", clust, "_", comp, "_noLabels.pdf")))
  }
  print(volcano)
  dev.off()
}


#' scRNA_pathways
#' @description Runs the pathway analysis for a markers DE file for\
#'  a single cell object 
#'
#' @param markers the markers object from performing FindMarkers() \
#' in memory, NOT a file path
#' @param orgdb org.Hs.eg.db or org.Mm.eg.db, for determining if we \
#'  need to switch gene names to human or not
#' @param outfolder the file path we'll output the results to 
#' @param type either "hallmark" or "go" pathways
#' @param paths the pathways we want to run 
#' @param pathname what to call the path file output once we make a table 
#'
#' @return NULL
scRNA_pathways <- function(markers,
                           orgdb,
                           type = "hallmark",
                           paths = NULL,
                           pathname,
                           outfolder) {
  
  markers <- markers[order(markers$avg_log2FC, decreasing = T),]
  rank <- markers[,"avg_log2FC"]
  names(rank) <- rownames(markers)
  
  # get paths loaded
  if (is.null(paths)) {
    if (type == "hallmark") {
      paths <- gmtPathways(file.path(GSEA_REF_DIR, 
                                     "h.all.v7.4.symbols.gmt.txt"))
    } else if (type == "go") {
      paths <- gmtPathways(file.path(GSEA_REF_DIR, 
                                     "c5.all.v7.4.symbols.gmt.txt"))
    }
  }
  
  # convert rank file to human genes if it's mouse
  if (orgdb == "org.Mm.eg.db") {
    rank <- convertMouseRanks2Hum(rank)
  }

  # run fgsea 
  fgsea_res <- fgsea(pathways = paths,
                     stats = rank,
                     minSize = 15,
                     maxSize = 500,
                     nperm= 1000)
  
  # by default take the top 60 pathway rows 
  if (nrow(fgsea_res) > 60) {
    fgsea_res <- fgsea_res[order(abs(fgsea_res$NES),decreasing=T),]
    fgsea_res <- head(fgsea_res, n=60)
  }
  fgsea_res$leadingEdge <- do.call(rbind, 
                                       lapply(as.list(fgsea_res$leadingEdge), 
                                                 paste, collapse=", "))
  fgsea_res <- as.data.frame(fgsea_res)
  fgsea_res$pathway <- factor(fgsea_res$pathway, 
            levels=fgsea_res$pathway[order(fgsea_res$NES,decreasing=F)])
  
  # write table 
  write.table(fgsea_res, file.path(outfolder, paste0(pathname, ".tsv")),
              sep = "\t", quote = F, col.names = NA)
  fgsea_res$significant <- ifelse(fgsea_res$padj < 0.05, "sig", "no")
  barp_cols <- c(sig="red", no="royalblue")
  
  # plot 
  pp <- ggplot(data=fgsea_res, aes(x=NES, y=pathway)) + 
    geom_bar(aes(fill=significant),
             stat="identity", orientation="y") + 
    scale_fill_manual(values=barp_cols) + 
    theme_bw() 
  
  pdf(file.path(outfolder,paste0(pathname, "_plot.pdf")))
  print(pp)
  dev.off()
  
  return(fgsea_res)
}



#single highlight umap subclusters for microglia, mono_dc, nk_t
#' Highlights cells of Umaps by condition
#' @description Loops through all unique factors of a condition and outputs highlighted umaps for each of those factors
#' @param seurat_object The seurat object you're making a UMAP for
#' @param coi The column name that we'll find in meta.data where we want to display the highlighted cells for
#' @param name What we're going to call the background sections of the umap
#' @param outdir The directory we're outputting to
#' @return NULL
umap_highlight <- function(seurat_object, coi, name = "background", outdir, foc_col = "#ff3300") {
  if (!(coi %in% colnames(seurat_object@meta.data))) {
    print("column not found in seurat object, check column names")
    return(NULL)
  }
  cois <- sort(unique(seurat_object@meta.data[,coi]))
  for (i in 1:length(cois)) {
    seurat_object@meta.data$highlight <- name
    seurat_object@meta.data$highlight[which(seurat_object@meta.data[,coi] == cois[i])] <- paste0(coi, " ", cois[i])
    seurat_object@meta.data$highlight <- 
      factor(seurat_object@meta.data$highlight, levels = c(paste0(coi, " ", cois[i]), "background"))
    #print figure
    UMAPPlot(seurat_object, group.by = "highlight", cols = c(foc_col, LIGHTGRAY))
    ggsave(paste0(outdir, "/umap_highlight_", coi, "_", cois[i], ".pdf"))
  }
  return(NULL)
}

#' @author https://stackoverflow.com/questions/15660829/how-to-add-a-number-of-observations-per-group-and-use-group-mean-in-ggplot2-boxp
#' @param vshift multiplier for shifting the annotation vertically in the figure
give.n <- function(x, vshift=1.05){
  return(c(y = median(x)*vshift, label = length(x))) 
}

# function for mean labels
#' @author https://stackoverflow.com/questions/15660829/how-to-add-a-number-of-observations-per-group-and-use-group-mean-in-ggplot2-boxp
#' @param vshift multiplier for shifting the annotation vertically in the figure
mean.n <- function(x,vshift=0.97){
  return(c(y = median(x)*vshift, label = round(mean(x),2))) 
}

# wrapper for violin plot function in seurat, which adds the number of cells for each violin plot group 
VlnPlot_withnums <- function(object, 
                             features,
                             group.by,
                             ...
) {
  order = sort(as.character(unique(object@meta.data[,group.by])))
  object@meta.data[,group.by] <- factor(object@meta.data[,group.by],
                                        levels=order)
  p1 <- VlnPlot(object = object, 
                features = features,
                group.by = group.by,
                assay="RNA",
                ...)
  p1 <- p1 + stat_summary(fun.data = give.n, geom = "text", fun = median)
  return(p1)
}







#' Runs a standard slingshot analysis with a seurat object as input
#' @param seurat_object 
#' @return pseudotime gene heatmap
#' 
slingshotAnalysis <- function(seurat_object) {
  
}

#' UMAP with KDE overlay
#' @description Makes a umap with kernel density estimation path overlay of 50th,85th,and 95th quantiles
#' @param seurat_object The seurat object with umap already calculated
#' @param gene The gene to be used for feature plotting
#' @param contour.input The contours to be made in the figure
#' @param reduction The dimensionality reduction in the object ("tsne" or "umap")
#' @note Taken from here: https://stackoverflow.com/questions/23437000/how-to-plot-a-contour-line-showing-where-95-of-values-fall-within-in-r-and-in
#' @return No return
dimred_withKDE <- function(seurat_object,genes,reduction,
                         contour.input = c(50,85,95)) {
  
  dimred.withmeta <- seurat_object@reductions[[reduction]]@cell.embeddings
  if (length(genes) > 1) {
    message("more than 1 gene, creating module score...")
    seurat_object <- AddModuleScore(seurat_object, features = list(genes), name = "genes",
                                    assay = "RNA")
    genes = "genes"
    colnames(seurat_object@meta.data)[ncol(seurat_object@meta.data)] <- genes
    dimred.withmeta <- merge(dimred.withmeta, seurat_object@meta.data, by = "row.names")
    rownames(dimred.withmeta) <- dimred.withmeta[,1]
    dimred.withmeta <- dimred.withmeta[,-1]
    dimred.withmeta <- dimred.withmeta[order(dimred.withmeta[,genes],decreasing = F),]
    
  } else {
    message("only 1 gene, fetching data...")
    dimred.withmeta <- merge(dimred.withmeta, seurat_object@meta.data, by = "row.names")
    rownames(dimred.withmeta) <- dimred.withmeta[,1]
    dimred.withmeta <- dimred.withmeta[,-1]
    genes <- genes[1]
    dimred.withmeta <- merge(dimred.withmeta, FetchData(seurat_object,genes,slot="data"), by = "row.names")
    rownames(dimred.withmeta) <- dimred.withmeta[,1]
    dimred.withmeta <- dimred.withmeta[,-1]
    dimred.withmeta <- dimred.withmeta[order(dimred.withmeta[,genes],decreasing = F),]
  }
  #weighted KDE because I don't care about where points are as much as where expression is
  set.seed(777)
  if (reduction == "umap") {
    colname = "UMAP"
  } else if (reduction == "tsne") {
    colname = "tSNE"
  }
  kd <- suppressWarnings(ks::kde(dimred.withmeta[,grep(colname, colnames(dimred.withmeta))], w=dimred.withmeta[,genes]) )
  contours <- list()
  if (length(contour.input) == 1) {
    col.grad <- "#000000" #just make it black if we only want 1 contour
  } else {
    col.grad <- scales::seq_gradient_pal("#888888", "#000000")(seq(0,1,length.out=length(contour.input)))
  }
  for (i in seq(1,length(contour.input))) {
    contours[[as.character(contour.input[i])]] <- with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                                         z=estimate, levels=cont[paste0(contour.input[i],"%")])[[1]]) #cont is a part of the kd object
    contours[[as.character(contour.input[i])]] <- data.frame(contours[[as.character(contour.input[i])]])
    contours[[as.character(contour.input[i])]] <- geom_path(data=contours[[as.character(contour.input[i])]], aes(x, y),
                                                            color = col.grad[i], size=0.3) 
  }
  
  ggplot(data=dimred.withmeta, aes_string(x=paste0(colname, "_1"),y=paste0(colname, "_2"),color=genes)) +
    geom_point() +
    contours + 
    scale_color_gradient(low = "#E5E5E5", high = "#F95738") + #"#FF006E";  #F95738 is orange
    theme(panel.background = element_rect(fill="white",color="black"))
}


#' Feature Plot with Boundaries
#' @description Draws boundaries around all the clusters in a dim.red., with expression pattern overlayed, doesn't work great...
#' 
#' @param seurat_object The seurat object with umap already calculated
#' @param gene The gene to be used for feature plotting
#' @param reduction The dimensionality reduction in the object ("tsne" or "umap")
#' @param coi The cluster condition of interest, i.e. where the splines are drawn
#'
#' @return
feature_with_boundaries <- function(seurat_object,genes,reduction,coi) {
  
  dimred.withmeta <- seurat_object@reductions[[reduction]]@cell.embeddings
  if (length(genes) > 1) {
    message("more than 1 gene, creating module score...")
    seurat_object <- AddModuleScore(seurat_object, features = list(genes), name = "genes",
                                    assay = "RNA")
    genes = "genes"
    colnames(seurat_object@meta.data)[ncol(seurat_object@meta.data)] <- genes
    dimred.withmeta <- merge(dimred.withmeta, seurat_object@meta.data, by = "row.names")
    rownames(dimred.withmeta) <- dimred.withmeta[,1]
    dimred.withmeta <- dimred.withmeta[,-1]
    dimred.withmeta <- dimred.withmeta[order(dimred.withmeta[,genes],decreasing = F),]
    
  } else {
    message("only 1 gene, fetching data...")
    dimred.withmeta <- merge(dimred.withmeta, seurat_object@meta.data, by = "row.names")
    rownames(dimred.withmeta) <- dimred.withmeta[,1]
    dimred.withmeta <- dimred.withmeta[,-1]
    genes <- genes[1]
    dimred.withmeta <- merge(dimred.withmeta, FetchData(seurat_object,genes,slot="data"), by = "row.names")
    rownames(dimred.withmeta) <- dimred.withmeta[,1]
    dimred.withmeta <- dimred.withmeta[,-1]
    dimred.withmeta <- dimred.withmeta[order(dimred.withmeta[,genes],decreasing = F),]
  }
  #weighted KDE because I don't care about where points are as much as where expression is
  set.seed(777)
  if (reduction == "umap") {
    colname = "UMAP"
  } else if (reduction == "tsne") {
    colname = "tSNE"
  }
  clusters = sort(unique(dimred.withmeta[,coi]))
  col.grad <- scales::hue_pal()(length(clusters))
  contours <- list()
  for (i in seq(1,length(clusters))) {
    dimred.withmeta[,coi] <- as.numeric(as.character(dimred.withmeta[,coi]))
    dimred.withmeta.subset <- dimred.withmeta[which(dimred.withmeta[,coi] == clusters[i]),]
    kd <- suppressWarnings(ks::kde(dimred.withmeta.subset[, 
                  grep(colname, colnames(dimred.withmeta))]) )
    contours[[paste0("cluster_", clusters[i])]] <- with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                                        z=estimate, levels=cont["5%"])[[1]]) #5% looks the best
    contours[[paste0("cluster_", clusters[i])]] <- data.frame(contours[[paste0("cluster_", clusters[i])]])
    contours[[paste0("cluster_", clusters[i])]] <- geom_path(data=contours[[paste0("cluster_", clusters[i])]],
                                                            aes(x, y), color = col.grad[i], size=0.5) 
  }
  
  ggplot(data=dimred.withmeta, aes_string(x=paste0(colname, "_1"),y=paste0(colname, "_2"),color=genes)) +
    geom_point() +
    contours + 
    scale_color_gradient(low = "#E5E5E5", high = "black") + #"#FF006E";  #F95738 is orange
    theme(panel.background = element_rect(fill="white",color="black"))
}

#' two_object_sample_correlation
#'
#' @description Prints a correlation heatmap with Pearson R values \
#' for two objects' samples split by each given COI
#' @note Using spearman correlation uses raw cell numbers instead of percentages of sample
#' @note Version 3: returns permutation test results included in pheatmap 
#' @param SOs A list of seurat objects to sample so1 and so2 from 
#' @param COIs A condition of interest for each object to split heatmap by EX: COIs <- list(so1_coi, so2_coi)
#' @param statmethod The correlation statistic, either pearson or spearman 
#' @param cellpcts Boolean, whether we want percentage of orig.ident or raw cell numbers
#' @param num_sim The number of simulations to run for the correlation permutation test for p-values
#' @param noplot Whether to output a figure or not
#' @param with_permtest Whether to add the p-values to the figure if noplot = F, or to the returned lists if noplot=T
#' @param norm_vector An external vector to normalize sample counts by, for instance, num CD45+ cells from larger object \
#' (only works with cellpcts=F)
#' @return NULL (prints the figure instead)
two_object_sample_correlation <- function(SOs, COIs, statmethod="pearson", 
                                          cellpcts=T, num_sim = 10000,
                                          noplot = F, with_permtest = T, 
                                          norm_vector = NULL
                                          # sample_col=NULL
                                          ) {
  so1 <- SOs[[1]]
  so2 <- SOs[[2]]
  coi1 <- COIs[[1]]
  coi2 <- COIs[[2]]
  
  if (statmethod != "pearson" && statmethod != "spearman") {
    print("statmethod must either be pearson or spearman, breaking...")
    return(NULL)
  }
  # if (is.null(sample_col)) {
  #   sample_col = "orig.isdent"
  # } else {
  #   if ((!sample_col %in% colnames(so1)) || (!sample_col %in% colnames(so2))) {
  #     stop("check sample_col, one object doesn't have this column...")
  #   }
  # }
  if (!(all(sort(unique(so1@meta.data$orig.ident)) == sort(unique(so2@meta.data$orig.ident))))) {
    print("object orig.ident's do not match...")
    return(NULL)
  } else {
    name.ordering = sort(unique(so1@meta.data$orig.ident))
  }
  if (so1@project.name == "SeuratProject") {
    name1 = "object1"
  } else {
    name1 = so1@project.name
  }
  if (so2@project.name == "SeuratProject") {
    name2 = "object2"
  } else {
    name2 = so2@project.name
  }
  cor.df <- data.frame(matrix(nrow = length(unique(so1@meta.data[,coi1])),
                              ncol = length(unique(so2@meta.data[,coi2]))))
  clusters1 <- sort(unique(so1@meta.data[,coi1]))
  clusters2 <- sort(unique(so2@meta.data[,coi2]))
  rownames(cor.df) <- paste0(name1, "_", clusters1)
  colnames(cor.df) <- paste0(name2, "_", clusters2)
  if (with_permtest) {
    cor.df.p <- cor.df #make a copy to use for p-values from permutation test
  }
  for (i in clusters1) {
    for (j in clusters2) {
      if (cellpcts) {
        cur.so1 <- table(factor(so1@meta.data$orig.ident[which(so1@meta.data[,coi1]==i)], levels=name.ordering)) / 
          table(so1@meta.data$orig.ident)
        cur.so2 <- table(factor(so2@meta.data$orig.ident[which(so2@meta.data[,coi2]==j)], levels=name.ordering)) / 
          table(so2@meta.data$orig.ident)
      } else {
        if (is.null(norm_vector)) {
          cur.so1 <- table(factor(so1@meta.data$orig.ident[which(so1@meta.data[,coi1]==i)], levels=name.ordering))
          cur.so2 <- table(factor(so2@meta.data$orig.ident[which(so2@meta.data[,coi2]==j)], levels=name.ordering))
        } else {
          if (length(intersect(names(norm_vector), name.ordering)) == length(name.ordering)) { #check if the same elements are in both
            norm_vector <- norm_vector[match(name.ordering, names(norm_vector))]
          } else {
            stop("norm_vector names don't match names provided by objects")
          }
          cur.so1 <- table(factor(so1@meta.data$orig.ident[which(so1@meta.data[,coi1]==i)], levels=name.ordering)) / 
            norm_vector
          cur.so2 <- table(factor(so2@meta.data$orig.ident[which(so2@meta.data[,coi2]==j)], levels=name.ordering)) / 
            norm_vector
        }
      }
      cur.so1 <- cur.so1[match(name.ordering, names(cur.so1))]
      cur.so2 <- cur.so2[match(name.ordering, names(cur.so2))]
      if (with_permtest) {
        set.seed(415)
        cur.permtest <- perm.cor.test(as.numeric(cur.so1), as.numeric(cur.so2), alternative = "two.sided", method = statmethod, num.sim = num_sim)
        cor.df.p[paste0(name1, "_", i),paste0(name2, "_", j)] <- cur.permtest$p
      }
      cor.df[paste0(name1, "_", i),paste0(name2, "_", j)] <- cor(cur.so1, cur.so2, method = statmethod)
    }
  }
  if (noplot) {
    if (with_permtest) {
      return(list("correlations"=cor.df, "permutation_test_pvalues"=cor.df.p))
    } else {
      return(list("correlations"=cor.df, "permutation_test_pvalues"=c()))
    }
  } else {
    if (with_permtest) {
      return(pheatmap(cor.df, display_numbers=cor.df.p, cluster_cols = F, cluster_rows = F, angle_col = 315))
    } else {
      return(pheatmap(cor.df, cluster_cols = F, cluster_rows = F, angle_col = 315))
    }
  }
}

#' get_sample_correlation_list
#'
#' @description VERY similar to two_object_sample_correlation() except this outputs a list 
#' of sample %s to use in scatterplot
#' @note Using spearman correlation uses raw cell numbers instead of percentages of sample
#' 
#' @param SOs A list of seurat objects to sample so1 and so2 from 
#' @param COIs A condition of interest for each object to split heatmap by EX: COIs <- list(so1_coi, so2_coi)
#' @param cellpcts Boolean, whether we want percentage of orig.ident or raw cell numbers
#' @return NULL (prints the figure instead)
get_sample_correlation_list <- function(SOs, COIs, cellpcts=T) {
  so1 <- SOs[[1]]
  so2 <- SOs[[2]]
  coi1 <- COIs[[1]]
  coi2 <- COIs[[2]]
  
  if (!(all(sort(unique(so1@meta.data$orig.ident)) == sort(unique(so2@meta.data$orig.ident))))) {
    print("object orig.ident's do not match...")
    return(NULL)
  } else {
    name.ordering = sort(unique(so1@meta.data$orig.ident))
  }
  if (so1@project.name == "SeuratProject") {
    name1 = "object1"
  } else {
    name1 = so1@project.name
  }
  if (so2@project.name == "SeuratProject") {
    name2 = "object2"
  } else {
    name2 = so2@project.name
  }
  cor.df <- data.frame(matrix(nrow = length(unique(so1@meta.data[,coi1])),
                              ncol = length(unique(so2@meta.data[,coi2]))))
  clusters1 <- sort(unique(so1@meta.data[,coi1]))
  clusters2 <- sort(unique(so2@meta.data[,coi2]))
  rownames(cor.df) <- paste0(name1, "_", clusters1)
  colnames(cor.df) <- paste0(name2, "_", clusters2)
  cor.list <- list()
  for (i in clusters1) {
    for (j in clusters2) {
      if (cellpcts) {
        cur.so1 <- table(factor(so1@meta.data$orig.ident[which(so1@meta.data[,coi1]==i)], levels=name.ordering)) / 
          table(so1@meta.data$orig.ident)
        cur.so2 <- table(factor(so2@meta.data$orig.ident[which(so2@meta.data[,coi2]==j)], levels=name.ordering)) / 
          table(so2@meta.data$orig.ident)
      } else {
        cur.so1 <- table(factor(so1@meta.data$orig.ident[which(so1@meta.data[,coi1]==i)], levels=name.ordering))
        cur.so2 <- table(factor(so2@meta.data$orig.ident[which(so2@meta.data[,coi2]==j)], levels=name.ordering))
      }
      return(NULL)
      cur.so1 <- cur.so1[match(name.ordering, names(cur.so1))]
      cur.so2 <- cur.so2[match(name.ordering, names(cur.so2))]
      
      cor.temp <- cbind.data.frame(as.vector(cur.so1), as.vector(cur.so2))
      colnames(cor.temp) <- c(paste0(name1, "_", i),paste0(name2, "_", j))
      rownames(cor.temp) <- name.ordering
      fullname <- paste0(paste0(name1, "_", i), "__", paste0(name2, "_", j))
      cor.list <- append(cor.list, list(cor.temp))
      names(cor.list)[length(cor.list)] <- fullname
    }
  }
  return(cor.list)
}

#' get_cell_percents
#' @description Outputs a table of percentages of samples that come from cells in each coi
#'
#' @param so 
#' @param COIs 
#' @param raw Whether to get actual cell numbers instead of percentages of orig.ident
#'
#' @return
#' @export
#'
#' @examples
get_cell_percents <- function(so, coi, raw=F) {
  
  name.ordering = sort(unique(so@meta.data$orig.ident))
  cor.df <- data.frame(matrix(nrow = length(unique(so@meta.data[,coi])),
                              ncol = length(unique(so@meta.data$orig.ident))))
  clusters.in <- as.character(sort(unique(so@meta.data[,coi])))
  clusters.in.names <- paste0(coi, "_", clusters.in)
  clusters.orig <- as.character(sort(unique(so@meta.data$orig.ident)))
  rownames(cor.df) <- clusters.in
  colnames(cor.df) <- clusters.orig
  
  if (raw) {
    for (i in clusters.in) {
      cur <- table(factor(so@meta.data$orig.ident[which(so@meta.data[,coi]==i)], levels=name.ordering))
      cur <- cur[match(name.ordering, names(cur))]
      cor.df[i,] <- cur
    }
  } else {
    for (i in clusters.in) {
      cur <- table(factor(so@meta.data$orig.ident[which(so@meta.data[,coi]==i)], levels=name.ordering)) / 
        table(so@meta.data$orig.ident)
      cur <- cur[match(name.ordering, names(cur))]
      cor.df[i,] <- cur
    }
  }
  
  return(cor.df)
}


#' get_cell_percents_mod
#' @description Modification of get_cell_percents from scRNAseqFunctions.R where you can choose divisor AND dividend
#' @param so Seurat Object
#' @param coi1 Condition of interest 1, will be used as numerator
#' @param coi2 Condition of interest 2, will be used as denominator
#' @param raw Whether to output raw cell numbers or use percents (F=percents)
#'
#' @return
#' @export
#'
#' @examples
get_cell_percents_mod <- function(so, coi1, coi2="orig.ident", raw=F) {
  
  name.ordering = sort(unique(so@meta.data[,coi2]))
  cor.df <- data.frame(matrix(nrow = length(unique(so@meta.data[,coi1])),
                              ncol = length(unique(so@meta.data[,coi2]))))
  clusters.in <- as.character(sort(unique(so@meta.data[,coi1])))
  clusters.in.names <- paste0(coi2, "_", clusters.in)
  clusters.orig <- as.character(sort(unique(so@meta.data[,coi2])))
  rownames(cor.df) <- clusters.in
  colnames(cor.df) <- clusters.orig
  
  if (raw) {
    for (i in clusters.in) {
      cur <- table(factor(so@meta.data[,coi2][which(so@meta.data[,coi1]==i)], levels=name.ordering))
      cur <- cur[match(name.ordering, names(cur))]
      cor.df[i,] <- cur
    }
  } else {
    for (i in clusters.in) {
      cur <- table(factor(so@meta.data[,coi2][which(so@meta.data[,coi1]==i)], levels=name.ordering)) / 
        table(so@meta.data[,coi2])
      cur <- cur[match(name.ordering, names(cur))]
      cor.df[i,] <- cur
    }
  }
  
  return(cor.df)
}

# makes a barplot from the output of get_cell_percents_mod()
# 20220623 | slight update: changed to print the plot an return the dataframe of pcts instead of vice versa
cell_percents_barplot <- function(so, coi1, coi2="orig.ident", withBarplotCols = F) {
  x <- get_cell_percents_mod(so=so, coi1=coi1, coi2=coi2)
  x.m <- melt(x)
  x.m$coi2 <- rep(rownames(x), length(unique(coi2)))
  colnames(x.m) <- c(coi2, "percent_of_total_cells", coi1)
  # print(x.m)
  p <- ggplot(data=x.m, aes_string(x=coi2,y="percent_of_total_cells",
                            fill=coi1)) + 
    geom_bar(stat="identity") 
  if (withBarplotCols) {
    p <- p = scale_fill_manual(values=barplot_cols(length(unique(x.m[,coi1]))))
  } 
  print(p)
  return(x.m)
}


# correlation code ----


#' corr_gene_sig
#' correlates a gene with a signature in a selected seurat object and returns correlation scatter
#' @param so the seurat object to find correlation in
#' @param gene the gene used to compare correlation against signature "sig"
#' @param sig the vector of genes used to compare against "gene"
#' @param signame name of signature which will be used for figure at the end
#' @param plot boolean determining whether to plot a correlation scatter or not
#' @return  either a correlation scatter plot or the matrix used to make such plot
corr_gene_sig <- function(so, gene, sig, signame = "test_signature", plot=T) {
  
  so.temp <- AddModuleScore(so, features=list(sig), seed = 415, name="temp")
  m <- aggregate(so.temp@meta.data$temp1 ~ so.temp@meta.data$orig.ident, 
                 data=so.temp@meta.data, FUN="mean")
  m.names <- m[,1]
  m <- as.vector(m[,2])
  names(m) <- m.names

  m <- m[match(names(trem2_markers),m.names)]
  comp <- cbind.data.frame(trem2_markers, m)
  estimate <- summary(lm(comp$trem2_markers ~ comp$m))$coefficients[[2]]
  rsq <- summary(lm(comp$trem2_markers ~ comp$m))$r.squared
  if (estimate < 0) {
    pearsonr <- -(sqrt(rsq))
  } else {
    pearsonr <- sqrt(rsq)
  }
  p <- ggplot(data=comp, aes(x=trem2_markers, y=m, label=rownames(comp))) +
    geom_point() +
    geom_text(hjust=-0.2,vjust=-0.2) +
    geom_smooth(method="lm", se=F) +
    geom_abline(intercept = 0, slope=1) +
    labs(title=paste0("Average expression of TREM2 per tumor sample vs. \naverage expression of signature: ", name, " per tumor sample"),
         subtitle=paste0("Pearson r = ", round(pearsonr, 4)),
         x = "TREM2 Average Expression",
         y = name)
  return(p)
}

#' gene_correlation_heatmap
#' @description This creates a heatmap of the top n genes correlated with selected gene of interest, in a single cell object\
#' we can also split by high and low to get genes that are correlated with high expression of the selected gene 
#' @param so 
#' @param gene the gene we're correlating with 
#' @param splitby whether we want to split the object by a metadata column's  factors
#' @param cur_assay the assay ("RNA" or "integrated") we want to gather correlation expression data from
#' @param topp a percentage representing how many genes we want to check from each object tested
#' @param topn a number of top genes to take from each object used to make the plot
#' @return a pheatmap object
gene_correlation_heatmap <- function(so, gene, splitby=NULL, cur_assay="RNA", 
                                     topn = NULL, topp = NULL) {
  if (!is.null(splitby)) {
    if (!(splitby %in% colnames(so@meta.data))) {
      stop("can't find splitby column in seurat object meta data")
    } else {
      message(paste0("splitting objects by column: ", splitby))
      so.list <- SplitObject(so, split.by = splitby)
    }
  } 
  if (!is.null(topp) && !is.null(topn)) {
    stop("can't choose both topp and topn")
  }
  
  gcorr <- function(cur.so) {
    message(paste0("working on: ", unique(cur.so@meta.data[,splitby])))
    x.mat <- as.matrix(cur.so@assays[[cur_assay]]@data)
    g.mat <- as.numeric(x.mat[gene,])
    x.mat <- x.mat[grep(paste0("^", gene,"$"), rownames(x.mat), invert = T),] #get rid of the target gene
    cor.rna.mat <- apply(x.mat,1,function(y){cor(g.mat,y)})
    cor.rna.mat <- cor.rna.mat[order(cor.rna.mat, decreasing = T)]
    return(cor.rna.mat)
  }
  
  if (is.null(splitby)) {
    
    so.mat <- gcorr(so)
    so.mat <- so.mat[!is.na(so.mat)]
    if (!is.null(topn)) {
      so.mat <- head(so.mat, n=topn)
    }
    if (!is.null(topp)) {
      so.mat <- so.mat[which(so.mat > quantile(so.mat, probs=1-topp/100))]
    } 
    
    so.mat <- sort(so.mat, decreasing = T)
    cor_mat <- as.data.frame(so.mat)
    rownames(cor_mat) <- names(so.mat)
    
    # returns a heatmap  
    p.heat <- pheatmap(cor_mat, cluster_rows = F, cluster_cols=F, cellwidth = 48,cellheight = 12)
    return(p.heat)
    
  } else {
    so.list.mat <- lapply(so.list, gcorr)
    so.list.mat <- lapply(so.list.mat, function(x){return(x[!is.na(x)])})
    so.list.meds <- lapply(so.list.mat, median)
    largest <- which.max(so.list.meds)
    if (!is.null(topn)) {
      so.list.mat <- lapply(so.list.mat, head, n=topn)
    }
    if (!is.null(topp)) {
      so.list.mat <- lapply(so.list.mat, function(xy) {
        return(xy[which(xy > quantile(xy, probs=1-topp/100))])
      })
    } 
    
    so.list.mat.gnames <- unique(Reduce(union, lapply(so.list.mat,names))) # gets  all of the possible genes from the correlation vecs
    cor_mat <- setNames(data.frame(matrix(ncol=length(so.list.mat),
                                 nrow=length(so.list.mat.gnames)), 
                                 row.names = so.list.mat.gnames), names(so.list.mat))
    for (i in 1:length(so.list.mat)) {
      cor_mat[,i] <- so.list.mat[[i]][match(so.list.mat.gnames, names(so.list.mat[[i]]))]
    }
    # cor_mat[is.na(cor_mat)] <- 0
    scheme <- apply(cor_mat, 1, function(x){return(length(x[!is.na(x)]))}) #number of columns with a non NA correlation
    cor_mat <- cor_mat[which(scheme > (ncol(cor_mat)/2)),] # at least half of the celltypes have each gene 
    scheme <- apply(cor_mat, 1, function(x){return(length(x[!is.na(x)]))}) #again because the last one was based off the unfiltered cor_mat
    cor_mat <- cor_mat[order(scheme, decreasing = T),]
    
    # returns a heatmap  
    p.heat <- pheatmap(cor_mat, cluster_rows = F, cluster_cols = T,
                       cellwidth = 48, cellheight = 12)
    return(p.heat)
  }
  
}



####### PATHWAY FUNCTIONS ######
# some of these are copied over from the functionsUsed.R file

# for inputting plotenrichment easier
plot.enrichment <- function (geneset, pathway, ranked_list) {
  plotEnrichment(geneset[[pathway]], ranked_list) + labs (title = pathway)
}

#ranks <- tibble::deframe(ranks)
fgsea_object <- function(mypathway, rank, write_table = F) { #creates the fgsea object
  fgsea_obj <- fgsea(pathways = mypathway,
                     stats = rank,
                     minSize = 15,
                     maxSize = 500,
                     nperm= 1000)
  
}

processed_ranks <- function(rank) {
  if (class(rank$log2FoldChange) != "numeric") { #if the log2FoldChange is no numeric- makes it numeric
    rank$log2FoldChange <- as.numeric(rank$log2FoldChange)
    
  }
  if(dim(rank[duplicated(rank$Gene.name),]) [1] > 0) { #checks if there are duplicates
    dupe = rank[,"Gene.name"]
    dupes <- rank[duplicated(dupe) | duplicated(dupe, fromLast=TRUE),] #finds the gene name
    print("Data has duplicates, average of duplicate names was taken.")
    show(dupes) #shows the duplicates
    rank <- aggregate(.~Gene.name, FUN = mean, data = rank) #takes the mean of the duplicated rows
  }
  
  if (any(is.na(rank)) == T) { #if there are any NA's in the ranks this will remove them
    rank <- rank[complete.cases(rank), ]
  }
  rank <- tibble::deframe(rank) #makes ranks into a vector
  rank
}

#' signif_barplot
#'
#' @param fgseas_obj fgsea object 
#' @param pattern a pattern to grep for in all of the fgsea object's pathways
#' @param mytitle what to name the plot
#' @param stat which column in the fgsea object to threshold for
#' @param statcutoff the threshold value for the stat 
#' @param topn prevents the output from containing a thousand pathways by taking the top n pathways
#' @param discrete whether to color the bars if passing/notpassing or by continuous value
#' @return a ggplot object 
signif_barplot <- function(fgseas_obj, 
                           pattern = "", 
                           mytitle, 
                           stat="padj", 
                           statcutoff=0.05,
                           topn = 50,
                           discrete=T) {
  require(viridis)
  if (!(stat%in%c("padj", "pval"))) {
    stop("stat must be either \"padj\" or \"pval\"")
  }
  if (pattern != "") {
    fgseas_obj <- fgseas_obj[grep(pattern, fgseas_obj$pathway),]
    if (nrow(fgseas_obj) == 0) {
      warning(paste0("no pathways match this pattern (", 
                     pattern ,") printing an empty figure..."))
    }
  }
  fgseas_obj$pathway <- gsub("^[^_]*_*", "\\1", fgseas_obj$pathway) #removes up to and including the first underscore of each pathway
  fgseaRes_tidy <- fgseas_obj %>% as_tibble() %>% arrange(desc(NES))
  if (nrow(fgseas_obj) > topn) { 
    fgseaRes_tidy <- fgseaRes_tidy %>% slice_head(n=topn)
  }
  if (discrete) {
    ggplot(fgseaRes_tidy, aes(reorder(pathway, NES), NES)) +
      geom_col(aes(fill = padj < statcutoff)) + coord_flip() +
      labs(x = "Pathway", y = "Normalized Enrichment Score", title = mytitle, fill = "") +
      theme(axis.text.y = element_text(size = 5)) +
      scale_fill_discrete(name = paste0("padj <", statcutoff))
  } else {
    ggplot(fgseaRes_tidy, aes(reorder(pathway, NES), NES))+
      geom_col(aes(fill = padj)) + coord_flip() +
      viridis::scale_fill_viridis(option="A" , direction = -1, limits=c(0,1)) + 
      labs(x = "Pathway", y = "Normalized Enrichment Score", title = mytitle, fill = "") +
      theme(axis.text.y = element_text(size = 5)) 
  }
}

#' fgsea_barplot_hallmark
#'
#' @param rankfile_loc the location of the .rnk file 
#' @param hallmark.pathway the pathway object (default is hallmark, which needs to be loaded) 
#' @param GO.pathway the go pathway object
#' @param fromfile 
#'
#' @return
#' @export
#'
#' @examples
fgsea_barplot_hallmark <- function(rankfile_loc,
                                   hallmark.pathway=hallmark.pathway,
                                   GO.pathway=GO.pathway,
                                   fromfile=T
) {
  if (!fromfile) {
    ranks <- rankfile_loc
  } else {
    ranks <- read.table(rankfile_loc, header = T,
                        colClasses = c("character", "numeric"), stringsAsFactors = F)
  }
  rank <- processed_ranks(ranks)
  hallmark_object <- fgsea_object(hallmark.pathway, ranks, write_table = T)
  signif_barplot(hallmark_object, "hallmark genes")
}


####### MONOCLE FUNCTIONS ######

#' Start with either monocle or seurat object, process and return object with default monocle settings
#' 
#' @param sc_object The monocle/seurat object used to process
#' @param var_feats The variable features used to order the cells of the monocle object
#' 
#' 
#' @return The processed monocle object (after OrderCells())
monocle_proc <- function(sc_object, var_feats) {
  message("converting to monocle object...")
  if (as.character(class(sc_object)) == "Seurat") {
    print("identifed as a seurat object, converting to monocle...")
    monocle_object <- as.CellDataSet(sc_object)
  }
  message("performing normalization...")
  monocle_object <- estimateSizeFactors(monocle_object)
  monocle_object <- estimateDispersions(monocle_object)
  monocle_object <- setOrderingFilter(monocle_object, var_feats)
  message("calculating ddrtree...")
  monocle_object <- reduceDimension(monocle_object, max_components = 2,
                                     method = 'DDRTree')
  message(paste0("ordering cells with variable genes (",length(var_feats) ,")..."))
  monocle_object <- orderCells(monocle_object)
  return(monocle_object)
} 

#' Finding root node in monocle (amended from monocle website)
#' @description Finds the suggested root node in a monocle DDRTree object
#' @param cds The monocle object
#' @param root_node_name The name of the celltype/timepoint/treatment that would designate a 'root' e.g. undifferentiated cells for a celltype
#' @param coi Condition of Interest; the column name in phenotype data that was used to bifurcate trajectory
#' @return An integer: the state number that has the most cells in root_node_name of coi
GM_state <- function(cds, root_node_name, coi){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)[,coi])[,root_node_name]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}


#' For creating Pseudotime figures, represented in a percentage bar graph (by coi)
#'
#' @param monocle_obj The monocle object used to make the figure, must have already done orderCells()
#' @param coi The column in pData() that we're using to bin percentages of cells 
#' @param pseudotime_col The column in pData() that is used as the x-axis (pseudotime is default)
#' @param binsize The size of a bin used to calculate a percentage
#' @param overlap The amount of overlap there will be on each bin (by absolute pseudotime_col units)
#' 
#' @return A data frame to make a ggplot object representing the bar graph we're trying to make
#'
#' @note Freedman-Diaconis suggested: \
#' https://stats.stackexchange.com/questions/798/calculating-optimal-number-of-bins-in-a-histogram/862
psuedotime_percent <- function(monocle_obj, 
                               pseudotime_col = "Pseudotime",
                               coi = "orig.ident", 
                               binsize = NULL,
                               overlap = 0) {
  if (is.null(binsize)) { #calculate the Freedman-Diaconis suggested binwidth if missing
    binsize = 2 * IQR(pData(monocle_obj)[,pseudotime_col]) / 
      length(pData(monocle_obj)[,pseudotime_col])^(1/3)
  }
  if (overlap != 0) {
    if (binsize %% overlap == 0) {
      breaks = seq(0, max(pData(monocle_obj)[,pseudotime_col]), overlap)
      skip = binsize / overlap
    } else {
      print(paste0("binsize: ", binsize, 
                   " is not divisible by overlap: ", overlap , " using overlap = 0"))
      overlap = 0
      skip = 1
      breaks = seq(0, max(pData(monocle_obj)[,pseudotime_col]), binsize)
    }
  } else {
    skip = 1
    breaks = seq(0, max(pData(monocle_obj)[,pseudotime_col]), binsize)
  }
  nbins = length(breaks)
  all_coi <- unique(pData(monocle_obj)[,coi])
  pseudotime_dat <- c()
  pseudotime_col.max <- max(pData(monocle_obj)[,pseudotime_col])
  last = F
  for (i in seq(1,nbins-skip)) {
    break.last = breaks[i]
    break.cur = breaks[i+skip] 
    if (i == nbins) {
      break.cur = pseudotime_col.max
      last = T
    }
    pseudotime_dat <- rbind(pseudotime_dat, 
                             get_pct_melt(monocle_obj,  #probably should just change this so that we're feeding it the postiion we're at along with the entire breaks array
                                          pseudotime_col, 
                                          coi,
                                          break.last,
                                          break.cur,
                                          binsize,
                                          last)
    )
  }
  colnames(pseudotime_dat)[3] <- pseudotime_col #it's easier to just name it here...
  pseudotime_dat[,coi] <- as.factor(pseudotime_dat[,coi])
  # bp <- ggplot(data=pseudotime, aes_string(x = pseudotime_col, y="percent", fill=coi)) +
  #   geom_bar(stat = "identity", position = "stack")
  return(pseudotime_dat)
}

#' Utility function for pseudotime_percent(), creates a single bar in the final barchart
#'
#' @param monocle_obj The monocle object used to make the figure, must have already done orderCells()
#' @param pseudotime_col The column in pData() that is used as the x-axis (pseudotime is default)
#' @param coi The column in pData() that we're using to bin percentages of cells 
#' @param break.last The smaller value for the bin we're calculating
#' @param break.cur The larger value for the bin we're calculating
#' @param binsize The size of a bin used to calculate a percentage
#' @param last Whether this is the last barplot of the psuedotime_percent() loop
#' 
#' @return The melted data to be appended to the full table
get_pct_melt <- function(monocle_obj, 
                         pseudotime_col,
                         coi,
                         break.last,
                         break.cur,
                         binsize,
                         last) {
  cur.frame <- pData(monocle_obj)[
    which(pData(monocle_obj)[,pseudotime_col] > break.last & 
            pData(monocle_obj)[,pseudotime_col] < break.cur),]
  cur.tab <- table(cur.frame[,coi])
  all_coi <- unique(pData(monocle_obj)[,coi])
  if (length(cur.tab) < length(all_coi)) {
    tosettozero <- all_coi[!(all_coi %in% names(cur.tab))]
    for (x in tosettozero) {
      cur.tab <- c(cur.tab, x=0)
      names(cur.tab)[length(cur.tab)] <- x
    }
    cur.tab <- as.table(cur.tab)
  }
  if (sum(cur.tab) != 0) {
    cur.tab.pct <- round((cur.tab / sum(cur.tab)), 3)*100 
  } else {
    cur.tab.pct <- cur.tab
  }
  cur.tab.pct.melt <- suppressWarnings(melt(cur.tab.pct)) #large warning is annoying
  # cur.tab.pct.melt[,coi] <- rownames(cur.tab.pct.melt)
  if (all(unlist(lapply(cur.tab.pct.melt, is.numeric)) == c(TRUE, FALSE))) { 
    #this means percents are in 1st column and labels are in 2nd, need to flip  
    cur.tab.pct.melt <- cur.tab.pct.melt[,c(2,1)]
  }
  colnames(cur.tab.pct.melt) <- c(coi, "percent")
  if (last) {
    cur.tab.pct.melt <- cbind(cur.tab.pct.melt, "break"= (break.last + (binsize / 2))) #want to save: coi ident, %, (centerpoint between last and current break)     
  } else {
    cur.tab.pct.melt <- cbind(cur.tab.pct.melt, "break"=base::mean(c(break.last, break.cur))) #want to save: coi ident, %, (centerpoint between last and current break)
  }
  return(cur.tab.pct.melt)
}

#' Helpful function for plotting the percent of coi and State on top of each other
#'
#' @param monocle_obj The monocle object used to make the figure, must have already done orderCells()
#' @param binsize The size of a bin used to calculate a percentage
#' @param overlap The amount of overlap there will be on each bin (by absolute pseudotime_col units)
#' @param coi The column in pData() that we're using to bin percentages of cells 
#' @param lineplot Whether we want it plotted as a line plot or a barplot
#' 
#' @return A grid of 2 ggplot 2 objects stacked 
pseudotime_pct_plt <- function(monocle_obj,
                               coi = "orig.ident",
                               binsize = 0.5,
                               overlap = 0,
                               lineplot = F
                               ) {
  pseudotime.dat.coi <- psuedotime_percent(monocle_obj = monocle_obj, 
                                       pseudotime_col = "Pseudotime",
                                       coi = coi, 
                                       binsize = binsize,
                                       overlap = overlap)
  pseudotime.dat.state <- psuedotime_percent(monocle_obj = monocle_obj, 
                                       pseudotime_col = "Pseudotime",
                                       coi = "State", 
                                       binsize = binsize,
                                       overlap = overlap)
  pseudotime_col = "Pseudotime"
  # barplot of coi
  p.coi <- ggplot(data=pseudotime.dat.coi, aes_string(x = pseudotime_col, y="percent", fill=coi))
  p.state <- ggplot(data=pseudotime.dat.state, aes_string(x = pseudotime_col, y="percent", fill="State"))
  if (lineplot) {
    p.coi <- p.coi + geom_area(position = "stack") 
    p.state <- p.state + geom_area(position = "stack") 
  } else {
    p.coi <- p.coi + geom_bar(stat = "identity", position = "stack") 
    p.state <- p.state + geom_bar(stat = "identity", position = "stack") 
  }
  
  p.coi <- p.coi + theme(legend.position = "top")
  p.state <- p.state + theme(legend.position = "bottom")
  
  plot_grid(p.coi, p.state, nrow = 2)
}



# SPATIAL SEURAT -----

# use graph to identify relationships of cell-cell interaction
# load in all of the spot locations via csvs in raw data (<sample>/outs/spatial/tissue_positions_list.csv)
#' Adds an adjacency matrix to a seurat object using igraph package 
#' @description Returns a seurat object with an igraph added to seurat_obj@graphs$spatial
#' @param seurat_obj the seurat object that we're adding an adjacency list to 
#' @param file the file (usually called tissue_positions_list.csv) that 
#' @param meta colnames in seurat_obj@meta.data that we want to keep for spot meta information in the graph
#' @param debug prints messages that will help with debugging
#' @return the same seurat object, but with a graph object added to seurat_obj@graphs$spatial
graph_lattice_from_seurat <- function(seurat_obj, meta, debug=T) {
  requires(igraph)
  requires(dplyr)
  requires(Seurat)
  
  # lazy, just wait for it to error instead of checking if exists
  spot_loc <- cbind(rownames(seurat_obj@images$slice1@coordinates), 
                    seurat_obj@images$slice1@coordinates)
  
  # should have these columns every time, unless spaceranger wasn't run, or 10x changes the code output 
  colnames(spot_loc) <- c("barcode", 
                          "in_tissue", 
                          "array_row", 
                          "array_col", 
                          "pxl_row_in_fullres", 
                          "pxl_col_in_fullres")
  spot_loc <- spot_loc[which(spot_loc$in_tissue == 1),]
  
  # merge meta data with spot_loc
  spots_meta <- merge(spot_loc, 
                      seurat_obj@meta.data, 
                      by.x = "barcode",
                      by.y = "row.names")
  if (nrow(spots_meta) < nrow(spot_loc)) {
    stop("spot barcodes don't align, exiting...")
  }
  
  spots_meta <- spots_meta[,c("barcode", 
                              "in_tissue", 
                              "array_row", 
                              "array_col", 
                              "pxl_row_in_fullres", 
                              "pxl_col_in_fullres",
                              meta)]
  
  spots_meta <- spots_meta[order(spots_meta$array_row, 
                                 spots_meta$array_col, 
                                 decreasing = F),]
  newgraph <- make_empty_graph(directed=F) #initialize it the only way I know how...
  for (i in 1:nrow(spot_loc)) {
    v.attr <- as.list(spots_meta[i,]) # can just take all cols, already filtered
    newgraph <- add_vertices(graph=newgraph, nv=1,
                             attr = v.attr)
  }
  for (i in V(newgraph)) {
    cur <- c(V(newgraph)[i]$array_row,
             V(newgraph)[i]$array_col)
    
    # all necessary surrounding edges to add
    surrounding <- list(c(cur[1],cur[2]+2), # +2 because in the same row we skip a column in the array; same row, right
                        c(cur[1]+1,cur[2]+1), # below to the right
                        c(cur[1]+1,cur[2]-1)) # below to the left
    
    # need to check to make sure a candidate vertex exists before we make an edge with it 
    for (j in surrounding) {
      # if location is missing from the graph, skip this new vertex, it's probably a border edge
      if (!any(j[1] == V(newgraph)$array_row & 
               j[2] == V(newgraph)$array_col)){
        if (debug) message(paste0("surrounding vertex at position: (", j[1], ",", j[2], ") doesn't exist..."))
        next
      }
      
      # get the actual vertex in the graph (usually just 1 number, the position in the vertex array of newgraph)
      # can just get the order because that's the name of the vertex
      cur.sur <- which(V(newgraph)$array_row == j[1] & V(newgraph)$array_col == j[2])
      cur.edge <- c(i,cur.sur) # undirected: order doesn't matter...
      
      #####  shouldn't need to check for duplicate edges anymore because we're only using the 3 vertices in front of a vertex
      ## if this becomes a problem, just use overloaded function igraph::unique() on the graph at the end...
      
      # add the edge, we know it's fine to do so
      if (debug) message(paste0("adding edge: ", cur.edge[1], " - ", cur.edge[2]))
      newgraph <- add_edges(newgraph, cur.edge)
    }
  }
  seurat_obj@graphs$spatial <- newgraph
  message("done.")
  return(seurat_obj)
}


#' @title spat_from_seurat
#' @description  create a spatstat graph with a seurat object 
#' @param object the seurat object to be used
#' @param feats any features to be added to the meta data of the spatstat points
#' @param assay we want to use from the seurat object
#' @return a spatstat point pattern object
#' @importFrom spatstat.geom ppp 
spat_from_seurat <- function(object=NULL, 
                             feats=NULL, 
                             assay="spatial") {
  DefaultAssay(object) <- assay
  cf <- object@images$slice1@coordinates
  if (!is.null(feats)) {
    feat_df <- FetchData(object, feats)
  }
  # reverse the y values
  cf$imagerow <- range_y[2] - cf$imagerow
  
  range_x <- c(min(cf$imagecol), max(cf$imagecol))
  range_y <- c(min(cf$imagerow), max(cf$imagerow))
  
  # create the ppp
  spat_object <- ppp(x=cf$imagecol, 
                     y=cf$imagerow,
                     xrange=range_x,
                     yrange=range_y,
                     marks=feat_df)
  return(spat_object)
}






                        

# ggplot custom theme
theme_v1 <- function(){ 
  #font <- "Georgia"   #assign font family up front
  
  theme_bw() %+replace%    #replace elements we want to change
    
    theme(
      
      #grid elements
      # panel.grid.major = element_blank(),    #strip major gridlines
      # panel.grid.minor = element_blank(),    #strip minor gridlines
      # axis.ticks = element_blank(),          #strip axis ticks
      
      #since theme_minimal() already strips axis lines, 
      #we don't need to do that again
      
      #text elements
      # plot.title = element_text(             #title
      #   family = font,            #set font family
      #   size = 20,                #set font size
      #   face = 'bold',            #bold typeface
      #   hjust = 0,                #left align
      #   vjust = 2),               #raise slightly
      
      # plot.subtitle = element_text(          #subtitle
      #   family = font,            #font family
      #   size = 14),               #font size
      # 
      # plot.caption = element_text(           #caption
      #   family = font,            #font family
      #   size = 9,                 #font size
      #   hjust = 1),               #right align
      # 
      # axis.title = element_text(             #axis titles
      #   family = font,            #font family
      #   size = 10),               #font size
      # 
      # axis.text = element_text(              #axis text
      #   family = font,            #axis famuly
      #   size = 9),                #font size
      
      axis.text.x = element_text(            #margin for axis text
        margin=margin(5, b = 10),
        angle=45, 
        hjust=1,
        vjust=1)
      
      #since the legend often requires manual tweaking 
      #based on plot content, don't define it here
    )
}
