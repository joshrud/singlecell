library(ggrepel)
library(Seurat) #using v3.1.1
library(ggplot2)
# library(ggradar)
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
library(ks)

roxygen2::roxygenize(package.dir = "/Users/jrudolph/Repos/singlecell/")

source("~/Repos/singlecell/clusterGSEA.R") #add directory path, R doesn't know where this is

timeFmt <- function() {
  format(Sys.time(), "%y%m%d_%H%M%S")
}

LIGHTGRAY = "#D4D3D2"
GRAY <- "#e6e6e6"
PURPLE <- "#8533ff"
GREEN <- "#00cc00"
RED_ORANGE <- "#ff531a"
PURPLE_HSIEH <- "#7F559B"
BRIGHT_PURPLE <- "#9445FF"

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
    seurat_object <- SCTransform(seurat_object, vars.to.regress = "percent.mt") #use this to replace scaledata if it's not good enough
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


library(DirichletReg)
calcDirichlet <- function(counts, covariates) {
  counts$counts <- DR_data(counts)
  data <- cbind(counts, covariates)
  data$condition <- factor(data$condition, levels=c("myometrium", "med12pos", "med12neg"))
  fit = DirichReg(counts ~ condition, data, verbosity = 0)
  u = summary(fit)
  pvals = u$coef.mat[, 4] #grep('Intercept', rownames(u$coef.mat), invert=T)
  v = names(pvals)
  pvals = matrix(pvals, ncol=length(u$varnames))
  rownames(pvals) = gsub('condition', '', v[1:nrow(pvals)])
  # colnames(pvals) = u$varnames
  colnames(pvals) <- paste0("cluster_", sort(unique(seurat_object@meta.data$seurat_clusters)))
  rownames(pvals)[grep("Intercept", rownames(pvals))] <- "Intercept" #paste0(as.character(unique(data$condition)[which(!(unique(data$condition) %in% rownames(pvals)))]), " (Int)")
  return(pvals)
}

###
# barPlotsWmedian() {} compiles boxplot data from seurat object for printing in ggplot
### makes celltype proportion boxplots, x=cluster, y=% of sample, group=condition 
### condition & orig.ident REQUIRED for this function to work
# seurat_object: a seurat object with 'condition' column in meta.data, clustering performed, and samples named in orig.ident
# name: some name to use for output naming
# ordering: expected order of clustering factors
# conditional_ordering: listed conditions (1 for each sample) in expected order
## condition ordering ex: c(rep("myometrium", 5), rep("med12pos", 5), rep("med12neg", 3))
# celltype_col: column in meta.data used for grouping the data (usually 'celltype' or 'seurat_clusters')
# type: using median or mean style figures (median currently preferred)
# outdir: directory to output to 
# withDirichlet: an option to write a table with pvalues according to dirichlet multinomial regression (debatable, anova is better)
###
barPlotsWmedian <- function(seurat_object, name, ordering = NULL, #condition_ordering = NULL, # condition order assignment messes things up
                            celltype_col = NULL, type="median", stats = "anova",  outdir) { #colorscheme=NULL,
  # if (!is.null(colorscheme)) {
  #   if (!(length(colorscheme) == length(unique(seurat_object@meta.data[,celltype_col])))) {
  #     print("celltype col length and colorscheme don't match, breaking...")
  #     break
  #   }
  # }
  if (is.null(celltype_col)) {
    celltypeNums <- getCellNumbers(seurat_object, celltype_col = "seurat_clusters", sample_col = "orig.ident")
  } else {
    celltypeNums <- getCellNumbers(seurat_object, celltype_col = celltype_col, sample_col = "orig.ident")
  }
  if (is.null(ordering)) {
    ordering = sort(colnames(celltypeNums))
  }
  celltypeNums <- celltypeNums[,ordering] 
  celltypeNums.pctOfSmple <- round(t(celltypeNums) / colSums(celltypeNums), 3)*100 
  condition_ordering <- as.vector(sapply(ordering, function(x) {
    return(unique(seurat_object@meta.data$condition[which(seurat_object@meta.data$orig.ident == x)]))}))
  covariates <- data.frame("condition"=condition_ordering) #when celltypeNums uses orig.ident
  celltypeNums.wconds <- cbind(celltypeNums.pctOfSmple, covariates)
  
  #dirichlet regression for pvalues (probably should try using anova also)
  print(celltypeNums.pctOfSmple)
  print(covariates)
  if (stats=="dirichlet" || stats == "all") {
    covariates <- cbind("condition"=condition_ordering, ordering)
    colnames(celltypeNums) <- paste0("cluster", colnames(celltypeNums))
    counts <- data.frame(cbind(t(celltypeNums)))
    dirichletStats <- calcDirichlet(counts, covariates)
    write.table(dirichletStats, paste0(outdir, "/dirichlet_regression_pvals.txt"), quote = F, sep = "\t", col.names = NA)
  } 
  if (stats == "anova" || stats == "all") {
    counts <- data.frame(cbind(t(celltypeNums)))
    covariates <- data.frame("condition"=condition_ordering)
    counts.wconds <- cbind(counts[,order(as.numeric(as.matrix(colnames(counts))))], covariates)
    counts.wconds <- counts.wconds[order(counts.wconds$condition),]
    pval <- c()
    for (i in 1:(ncol(counts.wconds)-1)) { 
      cur.stats <- counts.wconds[,c(i,ncol(counts.wconds))]
      colnames(cur.stats) <- c("value", "condition")
      model=lm( cur.stats$value ~ cur.stats$condition )
      ANOVA=aov(model)
      
      # Tukey test to study each pair of condition :
      TUKEY <- TukeyHSD(x=ANOVA, 'cur.stats$condition', conf.level=0.95)
      tukey.stats <- TUKEY$`cur.stats$condition`
      tukey.stats <- cbind(tukey.stats[,ncol(tukey.stats)], colnames(counts.wconds)[i])
      pval <- rbind(pval, tukey.stats)
    }
    pval <- data.frame(pval)
    colnames(pval) <- c("adjusted_pvalue", "cluster")
    pval$comparison <- gsub("\\.[0-9]*$", "", rownames(pval))
    rownames(pval) <- 1:nrow(pval)
    write.table(pval, paste0(outdir, "/anova_tukey_pvals.txt"), quote = F, sep = "\t", col.names = NA)
  } 
  celltypeNums.wconds.melt <- melt(celltypeNums.wconds)
  colnames(celltypeNums.wconds.melt) <- c("condition", "cluster",  "percent_of_sample")
  celltypeNums.wconds.melt$cluster <- factor(celltypeNums.wconds.melt$cluster)  
  # if (length(unique(celltypeNums.wconds.melt$condition)) == length(unique(seurat_object@meta.data$condition)))  {
  #   celltypeNums.wconds.melt$condition <- factor(celltypeNums.wconds.melt$condition, levels = c("myometrium", "med12pos", "med12neg"))
  # }
  #add stats (mean and median)
  celltypeNums.wconds.melt$median <- NA
  celltypeNums.wconds.melt$mean <- NA
  celltypeNums.wconds.melt$total <- NA
  celltypeNums.wconds.melt$condition <- as.character(celltypeNums.wconds.melt$condition)
  celltypeNums.wconds.melt$cluster <- as.character(celltypeNums.wconds.melt$cluster)
  for (i in unique(celltypeNums.wconds.melt$cluster)) {
    for (j in unique(celltypeNums.wconds.melt$condition)) {
      cur.frame <- celltypeNums.wconds.melt[which(celltypeNums.wconds.melt$condition == j & celltypeNums.wconds.melt$cluster == i),]
      median <- median(cur.frame$percent_of_sample)
      mean <- mean(cur.frame$percent_of_sample)
      total <- sum(cur.frame$percent_of_sample)
      celltypeNums.wconds.melt[which(celltypeNums.wconds.melt$condition == j & celltypeNums.wconds.melt$cluster == i),"median"] <- median
      celltypeNums.wconds.melt[which(celltypeNums.wconds.melt$condition == j & celltypeNums.wconds.melt$cluster == i),"mean"] <- mean
      celltypeNums.wconds.melt[which(celltypeNums.wconds.melt$condition == j & celltypeNums.wconds.melt$cluster == i),"total"] <- total
    }
  }
  celltypeNums.stats <- summarySEwithin(celltypeNums.wconds.melt, measurevar="percent_of_sample", withinvars=c("condition","cluster"))
  celltypeNums.stats <- celltypeNums.stats[order(celltypeNums.stats$condition,celltypeNums.stats$cluster),]
  celltypeNums.wconds.melt$se <- NA
  for (i in unique(celltypeNums.wconds.melt$cluster)) {
    for (j in unique(celltypeNums.wconds.melt$condition)) {
      cur.frame <- celltypeNums.wconds.melt[which(celltypeNums.wconds.melt$condition == j & celltypeNums.wconds.melt$cluster == i),]
      se <- celltypeNums.stats[which(celltypeNums.stats$condition == j & celltypeNums.stats$cluster == i), "se"]
      celltypeNums.wconds.melt[which(celltypeNums.wconds.melt$condition == j & celltypeNums.wconds.melt$cluster == i),"se"] <- se
    }
  }
  
  condition_ordering.factor <- unique(as.vector(as.matrix(condition_ordering))) #this ordering should be able to be changed 
  celltypeNums.wconds.melt$condition <- factor(celltypeNums.wconds.melt$condition, levels = condition_ordering.factor)
  ggplotAndSaveBarplots(celltypeNums.wconds.melt, name, type, outdir)
  write.table(celltypeNums.wconds, paste0(outdir, "celltypePctsWconds_", name, ".txt"), quote=F, col.names = NA, sep = "\t")
}

###
# ggplotAndSaveBarplots() {} utilized by barPlotsWmedian() to print ggplots to dir 
# plotdata: input dataframe for plotting
# name: some name to use for output naming
# type: using median or mean style figures (median currently preferred)
# outdir: directory to output to 
###
ggplotAndSaveBarplots <- function(plotdata, name, type, outdir) {
  prevdir <- getwd()
  setwd(outdir)
  
  print(summary(plotdata))
  plotdata$cluster <- factor(plotdata$cluster, levels = sort(unique(as.numeric(as.character(plotdata$cluster))), decreasing = F))
  if (type == "median") {
    ggplot(plotdata, aes(x=cluster, y=median, fill = condition)) + 
      geom_bar(stat="identity", position = position_dodge(width=0.7),  width = 0.6) +
      geom_errorbar(aes(ymin=ifelse(median-se < 0, 0,median-se), ymax=median+se),
                    width=.2,                    # Width of the error bars
                    position=position_dodge(.7)) + 
      # scale_fill_manual(values=COLORSCHEME2[c(2,4,3)]) +
      ggtitle(paste0("Median ", name)) +
      xlab("Cluster") + 
      ylab("Percent of Sample") +
      theme(axis.text.x = element_text(vjust = 1, size = 14), 
            axis.title = element_text(face = "bold", size = 16, margin = margin(0.3,0.3,0.1,0.1)),
            panel.background = element_rect(fill="white"), 
            panel.grid = element_blank(), 
            axis.line.x = element_line(color = "black", size = 0.5),
            axis.line.y = element_line(color = "black", size = 0.5),
            plot.margin = unit(c(0.1,0.1,0.3,0.3), "in")) + 
      scale_y_continuous(expand = c(0,0))
    ggsave(paste0("barplot_median_", name, ".pdf"), width=12, height=7, dpi=300)
    
    #now make the median box plots
    ggplot(plotdata, aes(x=cluster, y=percent_of_sample, color = condition)) + 
      geom_boxplot() + 
      # scale_color_manual(values=COLORSCHEME2[c(2,4,3)]) +
      ggtitle(paste0("Median ", name)) +
      xlab("Cluster") + 
      ylab("Percent of Sample") +
      theme(axis.text.x = element_text(vjust = 1, size=14), 
            axis.title = element_text(face = "bold", size = 16, margin = margin(0.3,0.3,0.1,0.1)),
            panel.background = element_rect(fill="white", color = "black"), 
            panel.grid = element_blank(), 
            axis.line.x = element_line(color = "black", size = 0.5),
            axis.line.y = element_line(color = "black", size = 0.5),
            plot.margin = unit(c(0.1,0.1,0.3,0.3), "in")) +
      scale_y_continuous(expand = c(0,0))
    ggsave(paste0("boxplot_median_", name, ".pdf"), width=12, height=7, dpi=300)
    
  } else if (type == "mean") {
    #bar plot
    ggplot(plotdata, aes(x=cluster, y=mean, fill = condition)) + 
      geom_bar(stat="identity", position = position_dodge(width=0.7),  width = 0.6) +
      geom_errorbar(aes(ymin=ifelse(mean-se < 0, 0,mean-se), ymax=mean+se),
                    width=.2,                    # Width of the error bars
                    position=position_dodge(.7)) + 
      # scale_fill_manual(values=COLORSCHEME2[c(2,4,3)]) +
      ggtitle(paste0("Mean ", name)) +
      xlab("Cluster") + 
      ylab("Percent of Sample") +
      theme(axis.text.x = element_text(vjust = 1, size=14),
            axis.title = element_text(face = "bold", size = 16, margin = margin(0.3,0.3,0.1,0.1)),
            panel.background = element_rect(fill="white"), 
            panel.grid = element_blank(), 
            axis.line.x = element_line(color = "black", size = 0.5),
            axis.line.y = element_line(color = "black", size = 0.5), 
            plot.margin = unit(c(0.1,0.1,0.3,0.3), "in")) +
      scale_y_continuous(expand = c(0,0))
    ggsave(paste0("barplot_mean_", name, ".pdf"), width=12, height=7, dpi=300)
    
  } else if (type == "total") {
    ggplot(plotdata, aes(x=cluster, y=total, fill = condition)) + 
      geom_bar(stat="identity", position = position_dodge(width=0.7),  width = 0.6) +
      geom_errorbar(aes(ymin=ifelse(total-se < 0, 0,total-se), ymax=total+se),
                    width=.2,                    # Width of the error bars
                    position=position_dodge(.7)) + 
      # scale_fill_manual(values=COLORSCHEME2[c(2,4,3)]) +
      ggtitle(paste0("total ", name)) +
      xlab("Cluster") + 
      ylab("Percent of Sample") +
      theme(axis.text.x = element_text(vjust = 1, size=14), 
            panel.background = element_rect(fill="white"), 
            axis.title = element_text(face = "bold", size = 16, margin = margin(0.3,0.3,0.1,0.1)),
            panel.grid = element_blank(),
            axis.line.x = element_line(color = "black", size = 0.5),
            axis.line.y = element_line(color = "black", size = 0.5),
            plot.margin = unit(c(0.1,0.1,0.3,0.3), "in")) +
      scale_y_continuous(expand = c(0,0))
    ggsave(paste0("barplot_total_", name, ".pdf"), width=12, height=7, dpi=300)
  }
  setwd(prevdir)
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
#' @param plotinsig Whether we want to plot the genes above FDR threshold in volcano plot
#' @param withLabels Whether we want labels in the volcano plot 
#' @param padj_thresh Adjusted pvalue threshold
#' @param orgdb The database used for performing GO analysis, either "org.Mm.eg.db" or "org.Hs.eg.db"
#' @param split Whether to do an "UP" and "DOWN" table for markers and GO analysis (best to keep it on)
#' @return DE results between comps in each coi, with volcanos and GO results for each
#'
easyDE <- function(seurat_object, 
                   output, comps, 
                   coi="seurat_clusters", 
                   plotinsig = TRUE, 
                   withLabels = T,
                   padj_thresh=0.05, 
                   orgdb = "org.Mm.eg.db", 
                   split = TRUE) {
  req.packages <- c("clusterProfiler",
                    "org.Hs.eg.db",
                    "ggthemes",
                    "enrichplot",
                    "DOSE",
                    "ggplot2")
  pack.mat <- sapply(req.packages, FUN = function(x) {
             if (paste0("package:",x) %in% search()) {
               return(TRUE)
             } else return(FALSE)
           } )
  if (any(!pack.mat)) {
    print(paste0("packages: ", paste0(req.packages[pack.mat], collapse = " "), " need to be loaded or installed..."))
    return(NULL)
  } else {
    print("good to go :)")
  }
  
  if (!exists("go_analysis", where=sys.frame())) {
    print("GO function not loaded, breaking...")
    return(NULL)
  }
  clusters = unique(seurat_object@meta.data[,coi])
  for (c in clusters) {
    clusterfolder <- paste0(output, "/cluster_", c, "_DE")
    dir.create(clusterfolder)
    
    cluster.cells <- rownames(seurat_object@meta.data[which(seurat_object@meta.data[,coi] == c),])
    cluster.set <- subset(seurat_object, cells = cluster.cells) #only cells in this cluster
    for (comp in comps) { #
      
      print(paste0(timeFmt(), " ----- working on ", comp, " in cluster ", c, " of object ", seurat_object@project.name))
      compfolder <- paste0(clusterfolder, "/",comp)
      dir.create(compfolder)
      splitted <- unlist(strsplit(comp, "vs"))
      cond1 <- splitted[1]
      cond2 <- splitted[2]
      cond1.cells <- rownames(cluster.set@meta.data[which(cluster.set@meta.data$condition == cond1),])
      cond2.cells <- rownames(cluster.set@meta.data[which(cluster.set@meta.data$condition == cond2),])
      if (length(cond1.cells) <= 3| length(cond2.cells) <=3 ) {
        print(paste0("comp ", comp, " in cluster ", c, " has too few cells. Skipping"))
        next
      }
      
      cluster.set <- SetIdent(cluster.set, cells = cond1.cells, value=cond1)
      cluster.set <- SetIdent(cluster.set, cells = cond2.cells, value=cond2)
      notEnoughGenes <- tryCatch({
        cur.de <- FindMarkers(cluster.set, ident.1 = cond1, ident.2 = cond2, min.pct = 0.1, logfc.threshold = 0.1) #only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25
      }, error = function(err) {
        notEnoughGenes <- TRUE
        print(paste("Caught an error: ", err, " ...continuing on...", sep = ""))
        return(TRUE)
      })
      if (typeof(notEnoughGenes)!="logical") {
        notEnoughGenes <- FALSE
        write.table(cur.de, paste0(compfolder, "/DEmarkers_", c, "_", comp, ".txt"), quote = F, sep = "\t", col.names = NA)
      }
      if (notEnoughGenes) next #this is where not enough genes is actually checked for errors
      if (split) {
        cur.de.up <- cur.de[which(cur.de$avg_logFC > 0 & cur.de$p_val_adj < padj_thresh),]
        cur.de.down <- cur.de[which(cur.de$avg_logFC < 0 & cur.de$p_val_adj < padj_thresh),]
        write.table(cur.de.up, paste0(compfolder, "/DEmarkers_", c, "_",cond1, "_sigUP", ".txt"), quote = F, sep = "\t", col.names = NA)
        write.table(cur.de.down, paste0(compfolder, "/DEmarkers_", c, "_",cond1, "_sigDOWN", ".txt"), quote = F, sep = "\t", col.names = NA)
      }  
      
      # GO Analysis (from Andrew's script)
      gofolder <- paste0(compfolder, "/GO")
      dir.create(gofolder)
      if (split) {
        goUPfolder <- paste0(gofolder, "/UP")
        dir.create(goUPfolder)
        goDOWNfolder <- paste0(gofolder, "/DOWN")
        dir.create(goDOWNfolder)
        
        upres <- go_analysis(cur.de.up, 
                    outdir =  goUPfolder, 
                    orgdb = "org.Mm.eg.db",
                    singlecomp = T)
        downres <- go_analysis(cur.de.down, 
                    outdir =  goDOWNfolder, 
                    orgdb = "org.Mm.eg.db",
                    singlecomp = T)
        if (!upres) {
          print("no up GO results, deleting GO folder...")
          unlink(goUPfolder, recursive = T) #remove GO dir for cleanliness
        }
        if (!downres) {
          print("no up GO results, deleting GO folder...")
          unlink(goDOWNfolder, recursive = T) #remove GO dir for cleanliness
        }
        if (!upres && !downres) {
          print("no up/down GO results, deleting GO folder...")
          unlink(gofolder, recursive = T) #remove GO dir for cleanliness
        }
      } else {
        res <- go_analysis(cur.de, 
                    outdir =  gofolder, 
                    orgdb = "org.Mm.eg.db",
                    singlecomp = T)
        if (!res) {
          print("no GO results, deleting GO folder...")
          unlink(gofolder, recursive = T) #remove GO dir for cleanliness
        }
      }
      
      if (withLabels) {
        easyDE.volcano(cur.de = cur.de, c = c, comp = comp, withLabels = T, #with labels
                       plotinsig = F, compfolder = compfolder)
      } else {
        easyDE.volcano(cur.de = cur.de, c = c, comp = comp, withLabels = F, #without labels
                       plotinsig = F, compfolder = compfolder)
      }
    }  
  }
}


#' Volcano Plot code for easyDE
#' @description Prints the volcano plot for easyDE so that easyDE is easier to follow
#' @param cur.de Requires same columns to a markers text file, with "logPadjf" = -log(pvalue) and gene names in rownames
#' @param c current cluster from the easyDE function
#' @param comp current comparison string
#' @param withLabels Whether we want labels in volcano
#' @param plotinsig Whether we want to plot genes above significance threshold
#' @param compfolder The folder we're outputting to, already written in easyDE
#' @return NULL 
easyDE.volcano <- function(cur.de, c = c, comp, withLabels = TRUE, plotinsig = TRUE, compfolder) {
  
  splitted <- unlist(strsplit(comp, "vs"))
  cond1 <- splitted[1]
  cond2 <- splitted[2]
  
  # volcano plot ----
  cur.de$logPadjf <- -log10(cur.de$p_val_adj + 1e-300)
  volcdat <- cur.de[,grep("avg_logFC|logPadjf", colnames(cur.de))]
  cur.de$labels <- NA
  if (withLabels) {
    
    # topLeft <- c(-300,300)
    # topRight <- c(300,300)
    # for (i in 1:nrow(cur.de)) {
    #   lDist <- dist(rbind(topLeft, volcdat[i,]), method = "euclidean")
    #   rDist <- dist(rbind(topRight, volcdat[i,]), method = "euclidean")
    #   cur.de$lDist[i] <- lDist
    #   cur.de$rDist[i] <- rDist
    # }
    # closestRight = cur.de[order(cur.de$lDist, decreasing = F),]
    # closestLeft = cur.de[order(cur.de$lDist, decreasing = F),]
    # closestRight.genes <- rownames(head(closestRight, n=10))
    # closestLeft.genes <- rownames(head(closestLeft, n=10))
    
    #new strategy: label top 5 genes by significance, top 5 by FC, and lowest 5 by FC
    top5.p <- head(rownames(cur.de[order(cur.de$logPadjf, decreasing = T),]), n=5)
    top5.l2fc.up <- head(rownames(cur.de[order(cur.de$avg_logFC, decreasing = T),]), n=5)
    top5.l2fc.down <- head(rownames(cur.de[order(cur.de$avg_logFC, decreasing = F),]), n=5)
    
    if (nrow(cur.de) < 20) {
      print(paste0("cur.de for comp: ", comp, " of cluster: ",c, " in object: ", seurat_object@project.name, " has fewer than 15 DE genes"))
      cur.de$labels <- rownames(cur.de)
    } else {
      # cur.de$labels[which(rownames(cur.de) %in% c(closestRight.genes, closestLeft.genes))] <- rownames(cur.de[which(rownames(cur.de) %in% c(closestRight.genes, closestLeft.genes)),])
      cur.de$labels[which(rownames(cur.de) %in% c(top5.p, top5.l2fc.up, top5.l2fc.down))] <- 
        rownames(cur.de[which(rownames(cur.de) %in% c(top5.p, top5.l2fc.up, top5.l2fc.down)),])
    }
  }
  
  cur.de$coloring <- NA
  cur.de$coloring[which(cur.de$p_val_adj < 0.05 & cur.de$avg_logFC > 1)] <- "de"
  cur.de$coloring[which(cur.de$p_val_adj < 0.05 & cur.de$avg_logFC < 1)] <- "significant"
  cur.de$coloring[which(cur.de$p_val_adj > 0.05 & cur.de$avg_logFC > 1)] <- "variable"
  cur.de$coloring[which(cur.de$p_val_adj > 0.05 & cur.de$avg_logFC < 1)] <- "nothing"
  viridis.4 <- viridis(4)
  de.colors <- c("de"=viridis.4[1], "significant"=viridis.4[2], "variable"=viridis.4[3], "nothing"="#000000")
  
  cur.de$cumulative_diff_btw_clusters <- abs(cur.de$pct.2 - cur.de$pct.1)
  r <- max(abs(cur.de$avg_logFC)) + 1
  if (!plotinsig) {
    print(paste0("not plotting insignificant genes for comparison: ", cond1, " vs ", cond2))
    cur.de <- cur.de[which(cur.de$coloring != "nothing"),]
    if (nrow(cur.de) == 0) {
      print("no genes left after significance filtration, continuing...")
      return(NULL)
    }
  }
  volcano <- ggplot2::ggplot(data=cur.de, aes(x = avg_logFC, y = logPadjf, size=cumulative_diff_btw_clusters)) + #not yet sure how this should be colored
    geom_point(aes(alpha=0.5, color=coloring)) + 
    geom_text_repel(aes(label = labels), point.padding = 0.25, size = 3) +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", size=0.5, alpha = 0.4) +
    geom_vline(xintercept=c(-1,1), linetype="dashed", size=0.5, alpha = 0.4) +
    scale_color_manual(name="Diff. Expressed", values=de.colors) +
    xlab(paste0("log2 Fold change (", cond1, "/", cond2, ")")) + scale_x_continuous(limits = c(-r, r)) +
    ylab("-log10 adjusted p-value") +
    ggtitle(paste0(comp)) +
    guides(alpha=F) + 
    theme_bw() + theme(legend.key.size = unit(2, 'lines'),
                       text = element_text(size=10),
                       legend.text=element_text(size=8), legend.key.height = unit(2, 'lines'),
                       legend.position = "right")
  if (withLabels) {
    pdf(paste0(compfolder, "/volcano_", c, "_", comp, "_withLabels.pdf")) # (png)
  } else {
    pdf(paste0(compfolder, "/volcano_", c, "_", comp, "_noLabels.pdf")) # (png)
  }
  print(volcano)
  dev.off()
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
#' @param SOs A list of seurat objects to sample so1 and so2 from 
#' @param COIs A condition of interest for each object to split heatmap by EX: COIs <- list(so1_coi, so2_coi)
#'
#' @return NULL (prints the figure instead)
two_object_sample_correlation <- function(SOs, COIs) {
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
  for (i in clusters1) {
    for (j in clusters2) {
      cur.so1 <- table(so1@meta.data$orig.ident[which(so1@meta.data[,coi1]==i)]) / 
        table(so1@meta.data$orig.ident)
      cur.so2 <- table(so2@meta.data$orig.ident[which(so2@meta.data[,coi2]==j)]) / 
        table(so2@meta.data$orig.ident)
      cur.so1 <- cur.so1[match(name.ordering, names(cur.so1))]
      cur.so2 <- cur.so2[match(name.ordering, names(cur.so2))]
      cor.df[paste0(name1, "_", i),paste0(name2, "_", j)] <- cor(cur.so1, cur.so2)
    }
  }
  print(pheatmap(cor.df, display_numbers=T, cluster_cols = F, cluster_rows = F))
  return(NULL)
}




####### MONOCLE FUNCTIONS

#' Start with either monocle or seurat object, process and return object with default monocle settings
#' 
#' @param sc_object The monocle/seurat object used to process
#' @param var_feats The variable features used to order the cells of the monocle object
#' 
#' 
#' @return The processed monocle object (after OrderCells())
monocle_proc <- function(sc_object, var_feats) {
  
  if (as.character(class(sc_object)) == "Seurat") {
    print("identifed as a seurat object, converting to monocle...")
    monocle_object <- as.CellDataSet(sc_object)
  }
  monocle_object <- estimateSizeFactors(monocle_object)
  monocle_object <- estimateDispersions(monocle_object)
  monocle_object <- setOrderingFilter(monocle_object, var_feats)
  monocle_object <- reduceDimension(monocle_object, max_components = 2,
                                     method = 'DDRTree')
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




                        

