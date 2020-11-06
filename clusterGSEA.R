###########GSEA######################
# originally written by Andrew Shroeder 

#convert names
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggthemes)
library(enrichplot)
library(DOSE)
library(ggplot2)

#org.Mm.eg.db, org.Hs.eg.db
go_analysis <- function(seuratobject.markers, outdir, orgdb ,
                        singlecomp = F){
  prevdir <- getwd()
  db <- GOSemSim::load_OrgDb(orgdb)
  
  #convert from UNIPROT to symbol, ensembl, entrez id, and genename
  gsea.markers<- subset(seuratobject.markers) 
  if (nrow(gsea.markers) == 0) {
    print("not enough genes above threshold")
    return(FALSE)
  }
  if (singlecomp) {
    gsea.markers<- bitr(rownames(gsea.markers), fromType="SYMBOL", toType="ENTREZID", OrgDb=orgdb)
    colnames(gsea.markers)[1]<-'gene'
    gsea.markers<-merge(gsea.markers, seuratobject.markers, by.x='gene', by.y="row.names", all.x=T)
    #filter
    gsea.markers<- subset(gsea.markers)
    #order dataset
    # print(gsea.markers)
    gsea.markers<- gsea.markers[with(gsea.markers, order(-avg_logFC)),]
    ## feature 1: numeric vector
    geneList <- gsea.markers$avg_logFC
    if (length(geneList) < 25) {
      print(paste0("not enough genes, skipping..."))
      return(FALSE)
    }
    ## feature 2: named vector
    names(geneList) <- unique(as.character(gsea.markers$ENTREZID))
    ## feature 3: decreasing order
    geneList <- sort(geneList, decreasing = TRUE)
    enrichGO <- enrichGO(gene     = names(geneList),
                         #universe= names(geneList),
                         OrgDb        = db, 
                         ont          = "BP", #can change to MF or CC as well
                         keyType="ENTREZID",
                         pvalueCutoff = 0.1,
                         qvalueCutoff = 0.05,
                         readable=TRUE,
                         pool=FALSE)
    #assign(paste( 'gseGO.',i, sep = ""), gseGO)
    write.table(enrichGO@result, paste0(outdir, "/GO.txt"), sep = "\t", quote = F, col.names = NA)
    png(paste0(outdir, "/GO.png"),width=15,height=10,units = 'in', res =300)  
    
    #the showCategory parameter is # of BPs to plot
    print( barplot(enrichGO,x='Count', showCategory=30, title=paste('Enriched Biological Processes'), color='qvalue')+
             theme_tufte()+ylab('Genes in Set')+
             theme(axis.text.y = element_text(size=12))+
             theme(plot.title = element_text(hjust = 0.5))
    )
    dev.off()
    
  } else {
    gsea.markers<- bitr(gsea.markers$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb=orgdb)
    colnames(gsea.markers)[1]<-'gene'
    gsea.markers<-merge(gsea.markers, seuratobject.markers, by='gene',all.x=T)
    #filter
    # gsea.markers<- subset(gsea.markers, avg_logFC >0.25)
    
    #order dataset
    # print(gsea.markers)
    gsea.markers<- gsea.markers[with(gsea.markers, order(cluster, -avg_logFC)),]
    #start pipeline to create figure for cluster specific gsea
    for(i in sort(unique(gsea.markers$cluster)) ){
      
      ## feature 1: numeric vector
      geneList <- gsea.markers$avg_logFC[which(gsea.markers$cluster==i)]
      if (length(geneList) < 25) {
        print(paste0("not enough genes for cluster: ", i, " skipping..."))
        next #because of this, we'll return true for each time we run through a FindAllMarkers textfile
      }
      ## feature 2: named vector
      names(geneList) <- unique(as.character(gsea.markers$ENTREZID[which(gsea.markers$cluster==i)]))
      ## feature 3: decreasing order
      geneList <- sort(geneList, decreasing = TRUE)
      
      enrichGO <- enrichGO(gene     = names(geneList),
                           #universe= names(geneList),
                           OrgDb        = db, 
                           ont          = "BP", #can change to MF or CC as well
                           keyType="ENTREZID",
                           pvalueCutoff = 0.1,
                           qvalueCutoff = 0.05,
                           readable=TRUE,
                           pool=FALSE)
      #assign(paste( 'gseGO.',i, sep = ""), gseGO)
      write.table(enrichGO@result, paste0(outdir, "/GO.txt"), sep = "\t", quote = F, col.names = NA)
      png(paste0(outdir, "/GO.png"),width=15,height=10,units = 'in', res =300)  
      
      #the showCategory parameter is # of BPs to plot
      print( barplot(enrichGO,x='Count', showCategory=30, title=paste('Cluster',i,'Enriched Biological Processes'), color='qvalue')+
               theme_tufte()+ylab('Genes in Set')+
               theme(axis.text.y = element_text(size=12))+
               theme(plot.title = element_text(hjust = 0.5))
      )
      dev.off()
    }
  }
  return(TRUE)
}

#output is a folder with a figure for each cluster
