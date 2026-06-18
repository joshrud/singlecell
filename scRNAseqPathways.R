library(ggrepel)
library(Seurat) 
library(ggplot2)
library(dplyr)
library(reshape2)
library(ggthemes)
library(scales)
library(fgsea)
library(biomaRt)
library(ggrepel)
library("org.Mm.eg.db")
library("org.Hs.eg.db")

####### CONSTANTS ######
GSEA_REF_DIR = "/data/annotations/gsea"
LIGHTGRAY = "#D4D3D2"
GRAY <- "#e6e6e6"
PURPLE <- "#8533ff"
GREEN <- "#00cc00"
RED_ORANGE <- "#ff531a"
PURPLE_HSIEH <- "#7F559B"
BRIGHT_PURPLE <- "#9445FF"

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




#' @name easyDE 
#' @description Performs Differential Expression testing for all idents of a \
#' seurat object, outputs volcano plots, de marker table, and go analysis \
#' for all comparisons in each cluster
#'
#' @param seurat_object The seurat object.
#' @param output The directory we want to output DE results directory to 
#' @param comps The comparisons to perform in the form "cond1vscond2", \
#' comparing cond1 vs. cond2; both should be in "conditions" column in meta.data
#' @param coi The 'column of interest' in the meta.data that we want to split \
#' object by to perform DE testing, ex: "seurat_clusters" 
#' @param with_labels Whether we want labels in the volcano plot 
#' @param padj_thresh Adjusted pvalue threshold
#' @param org_db The database used for performing GO analysis, either \
#' "org.Mm.eg.db" or "org.Hs.eg.db", or "auto" for detecting from gene names \
#' which one s
#' @return DE results between comps in each coi, with volcanos and GO\
#'  results for each
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
  
  # notify 
  timestamp()
  print(paste0("Detected ", orgdb, ", as database."))
  
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
      scRNA_pathways(cur.de,  
                     outfolder =  go_folder, 
                     type = "go",
                     paths = NULL,
                     pathname = paste0(clusters[i],"_",comp, "_GO_paths"),
                     orgdb = orgdb)
      scRNA_pathways(cur.de, 
                     outfolder =  hall_folder, 
                     type = "hallmark", 
                     paths = NULL,
                     pathname = paste0(clusters[i],"_",comp, "_Hall_paths"),
                     orgdb = orgdb)
      
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
}



