# ==ReactomePA====================https://www.bioconductor.org/packages/release/bioc/vignettes/ReactomePA/inst/doc/ReactomePA.html#pathway-visualization
#http://yulab-smu.top/clusterProfiler-book/chapter12.html

#====INSTALL PACKAGES=====

# Lista de paquetes de funciones a instalar
.packages = c("dplyr","data.table","plyr","ggplot2", "tidyverse","tidyr","org.Hs.eg.db",
              "GOsummaries","ReactomePA","clusterProfiler","DOSE")

# Instala los paquetes sin? los tienes instalados
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) {install.packages(.packages[!.inst])
  BiocManager::install(c("GOsummaries","ReactomePA","clusterProfiler","DOSE","doseplot"))}

# Carga los paquetes sin? los tienes cargados
lapply(.packages, require, character.only=TRUE)
#==============================================================

# browseVignettes("ReactomePA") 

path <- "output/output_ingrid/"
gene_selected_Tabla <- read.delim(paste(path,"gene selected_in1.txt",sep = "")) # Open the data
rownames(gene_selected_Tabla) <- gene_selected_Tabla$gene_set


gene_selected_Tabla$FC[gene_selected_Tabla$coefficient_estimate < 0] = -1
gene_selected_Tabla$FC[gene_selected_Tabla$coefficient_estimate > 0] = 1

# Pathway enrichment analysis with the geneset "gene_set_263".

gene_selected_Tabla1 <- subset(gene_selected_Tabla, gene_selected_Tabla$Result=="gene_set_263")

x =  gene_selected_Tabla1$gene_set

eg = clusterProfiler::bitr(x, fromType="SYMBOL", 
                           toType="ENTREZID", 
                           OrgDb="org.Hs.eg.db") #  convert biological IDs to Entrez gene ID.
de <- eg$ENTREZID

x1 <- ReactomePA::enrichPathway(gene=de,pvalueCutoff=0.05, readable=T)
head(as.data.frame(x1))

barplot(x1, showCategory=25)
dotplot(x1, showCategory=25)

FC0 = as.data.frame(gene_selected_Tabla1[,c('gene_set','FC')])
names(FC0)[1] <- 'SYMBOL'

FC.df <- plyr::join_all(list(FC0,eg),by='SYMBOL')

FC <- FC.df$FC

names(FC)<-FC.df$ENTREZID

ReactomePA::cnetplot(x1, categorySize="pvalue", 
                     layout = "kk",legend_n=2,pie="count",
                     colorEdge = T, 
                     circular=T,showCategory=9,
                     line_scale = 0.5,pie_scale=1.5,
                     foldChange = FC)


# Pathway enrichment analysis with the geneset "Metabolism".

gene_selected_Tabla1 <- subset(gene_selected_Tabla, gene_selected_Tabla$Result=="Metabolism")

x =  gene_selected_Tabla1$gene_set

eg = clusterProfiler::bitr(x, fromType="SYMBOL", 
                           toType="ENTREZID", 
                           OrgDb="org.Hs.eg.db") #  convert biological IDs to Entrez gene ID.
de <- eg$ENTREZID

x1 <- ReactomePA::enrichPathway(gene=de,pvalueCutoff=0.05, readable=T)
head(as.data.frame(x1))

barplot(x1, showCategory=25)
dotplot(x1, showCategory=25)

FC0 = as.data.frame(gene_selected_Tabla1[,c('gene_set','FC')])
names(FC0)[1] <- 'SYMBOL'

FC.df <- plyr::join_all(list(FC0,eg),by='SYMBOL')

FC <- FC.df$FC

names(FC)<-FC.df$ENTREZID

ReactomePA::cnetplot(x1, categorySize="pvalue", 
                     layout = "kk",legend_n=2,pie="count",
                     colorEdge = T, 
                     circular=T,showCategory=9,
                     line_scale = 0.5,pie_scale=1.5,
                     foldChange = FC)