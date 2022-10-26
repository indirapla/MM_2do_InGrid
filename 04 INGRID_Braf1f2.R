#========MELANOMA 2do. B-raf mutation=========INGRID================================   
#  https://dongjunchung.github.io/INGRID/ 
#   vignette("INGRID-example")
#   INGRID: Integrative Genomics Robust iDentification of cancer subgroups ===
#   
#   INGRID is a statistical approach that integrates information from biological 
#   pathway databases with high-throughput genomic data to improve the robustness for
#   dentification and interpretation of molecularly-defined subgroups of cancer patients.

#====INSTALL PACKAGES=================================================================

# List of packages to install
.packages = c("devtools","BiocManager","tidyverse","ggplot2",
              "ggbiplot","doParallel","ReactomePA","clusterProfiler",
              "msigdbr","INGRID","reactome.db")

# Installing not installed packages
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) {install.packages(.packages[!.inst])
  BiocManager::install(c("ReactomePA","clusterProfiler","reactome.db"))
  devtools::install_github("vqv/ggbiplot",force = TRUE)
  devtools::install_github("dongjunchung/INGRID",force = T)
  #remotes::install_github("sunny-zq/INGRID")
 }

# Loading packages 
lapply(.packages, require, character.only=TRUE)
#=========================================================================================

#========Loading tables=====================================================

# protein expression table
Big_Tabla <- read.delim("data/data_ingrid/20191128_MM_TMT_all_isotope_correction_100threshold_normalized_ratios_with samples names_1_ipp2.txt") # Open the data
rownames(Big_Tabla) <- paste(Big_Tabla$Accession,".",Big_Tabla$Gene.name,sep = "")

Big_Tabla1 <- Big_Tabla %>% dplyr::select(contains("MM")) # selecting only columns that contain "MM" in the column name

head(Big_Tabla1)

# clinical data containing survival information
annotation_colum <- read.delim("data/data_ingrid/20200520ClinicalData_144samplesLundMM_REVISEDandUPDATEDtill2020_no_formulas.txt", header = TRUE,row.names = 1)
annotation_colum$sample <- row.names(annotation_colum)

annotation_colum$time.surv <- annotation_colum[,"days.from.sample.collection.to.censoring.death"]/365.35
annotation_colum$cat.surv[annotation_colum$time.surv < 3] = 1
annotation_colum$cat.surv[annotation_colum$time.surv >= 3] = 0

# Pathways table

string.sig.path <- read.csv("data/data_STRING-react/enrichment.RCTM.string.csv")
string.sig.path$term.ID2 <- paste("R-",string.sig.path$X.term.ID, sep = "")
ReactomePathwaysDb <- read.delim("data/data_reactome/ReactomePathwaysDB.txt", header = F) # Open the data


#======================================================

#===============================FUNCTIONS==============================================

# filtering ReactomeDB based on the list of pathways selected to be used in INGRID
   # x = list of pathways selected to be used in INGRID. e.g. "string.sig.path"
   # y = ReactomeDB        # Reactome database. 
   # x.col.ID = "term.ID2" # column name by which x table will be merged with y table
   # y.col.ID = "V2"       # column name by which y table will be merged with x table  
   
   # Example
   # x = string.sig.path
   # y = ReactomePathwaysDb
   # x.col.ID = "term.ID2"
   # y.col.ID = "V2"

short.ReactomeDB <- function(x,y, x.col.ID,y.col.ID){

  x.path <- x[,x.col.ID]
  
  y1 <- subset(y, y[,y.col.ID] %in% x.path)
  
  return(y1)
}
short.ReactomeDB <- short.ReactomeDB(x=string.sig.path,
                                     y=ReactomePathwaysDb,
                                     x.col.ID = "term.ID2",
                                     y.col.ID = "V2")

#===converting ReactomePathwaysDb format to the one required for INGRID
    # x = ReactomePathwaysDb unformated

Path.list.format.function <- function(x){
  rownames(x) <- x[,2]  # naming rows with the pathway code (e.g R-HSA-164843) which are in column 2 
  path.list <- list()
 
  for (i in 1:nrow(x)) {
    
    x1<-x[i,3:ncol(x)]
    path <- x1[1,]
    path <- t(path)
    path <- as.character(path[,1])
    path <- as.character(path[!path %in% c("")])
    path.list[[i]] <- path
    names(path.list)[i] <- as.character(x[i,1])

  } 
  return(path.list)
}
short.ReactomeDB.list <- Path.list.format.function(x = short.ReactomeDB)


# Selecting only Braf samples==========================
    # x = Protein expression table
    # y = Table with clinical information
    # colnames.x = "MM" ; columns that contain "MM" in the column name in the x table
    # y.var = "BRAF.cat" ; column in table y by with the samples in x will be selected
    # y.var.T = 1 ; criterias to be selected in the y.var. It could be also a vector
    # conditions.level = 1 ; can be only 1 or 2


select.samples <- function(x,y, colnames.x,y.var,
                           y.var.T,t.column, d.column,
                           conditions.level){
  
  y1 <- na.omit(y)
  x1 <- x %>% dplyr::select(contains(colnames.x)) # selecting only columns that contain "MM" in the column name
  
  x1.T <- as.data.frame(t(x1))
  colnames(x1.T) <- x$Gene.name
  x1.T$sample <- row.names(x1.T)
  
  x1.T1 <- plyr::join_all(list(y1,x1.T),by="sample",type="inner")
  
  row.names(x1.T1) <- x1.T1$sample
  
  if(conditions.level==1){
  x1.T1.n <-  subset(x1.T1, x1.T1[,y.var] %in% y.var.T)
  }else{
    x1.T1.n <-  subset(x1.T1, x1.T1[,y.var] %in% y.var.T)
    x1.T1.n <-  subset(x1.T1.n, BRAF.conc. %in% c("L","H"))
  }
  
  t <- as.numeric(x1.T1.n[,t.column]/365.35)   # tid in years
  d <- as.numeric(x1.T1.n[,d.column])
  
  x1.T1.n <- x1.T1.n[,(ncol(y1)+1):ncol(x1.T1.n)]
  
  prot.exp <- as.data.frame(t(x1.T1.n))
  
  prot.exp.m <- na.omit(prot.exp)
  
  data.list <- list(t,d,prot.exp.m)
  names(data.list) <- c("t","d","prot.exp.m")
  
  return(data.list)
  
}

y <- annotation_colum[,c("sample","cat.surv","Remarks.origil.num","BRAF.cat","BRAF.conc.",
                         "days.from.sample.collection.to.censoring.death")]

Braf.data.learn.list <- select.samples(x = Big_Tabla,
                                 y=y,
                                 colnames.x="MM",
                                 y.var = "BRAF.cat",
                                 y.var.T = 1,
                                 t.column = "days.from.sample.collection.to.censoring.death",
                                 d.column = "Remarks.origil.num",
                                 conditions.level = 1)

Braf.data.test.list <- select.samples(x = Big_Tabla,
                                       y=y,
                                       colnames.x="MM",
                                       y.var = "BRAF.cat",
                                       y.var.T = 1,
                                       t.column = "days.from.sample.collection.to.censoring.death",
                                       d.column = "Remarks.origil.num",
                                       conditions.level = 2)

#= === Z-SCORE function================

# x = Braf.data.list$prot.exp.m ; protein expression 
# dim.class = 2 ; dim.class can be only 1(rows) or 2(columns)

zscore.function <- function(x,dim.class){
  
z <- matrix(nrow = nrow(x),ncol = ncol(x),
                        dimnames = list(rownames(x),colnames(x)))
if(dim.class==1){
  mean = apply(x, 1, mean)
    sd = apply(x, 1, sd)
  
  for (j in 1:ncol(x)) {
     for (i in 1:nrow(x)){
      zscore <- (x[i,j]-mean[i]/sd[i]) 
      
      z[i,j] <- zscore
    }
    } 
      return(z)} 
else 
  if(dim.class==2){
    mean = apply(x, 2, mean)
      sd = apply(x, 2, sd)
    
    for (i in 1:nrow(x)) {
      for (j in 1:ncol(x)) {
          zscore <- (x[i,j]-mean[j]/sd[j]) 
    
      z[i,j] <- zscore
  }
  }
    return(z)}

else print("'dim.class' can only be  1(rows) or 2(columns)")
}

prot.exp.learn.m <- Braf.data.learn.list$prot.exp.m

prot.exp.test.m <- Braf.data.test.list$prot.exp.m

zscore.prot.exp.learn.m <- zscore.function(x=prot.exp.learn.m,dim.class = 1)
zscore.prot.exp.learn.m <- as.data.frame(zscore.prot.exp.learn.m)

zscore.prot.exp.test.m <- zscore.function(x=prot.exp.test.m,dim.class = 1)
zscore.prot.exp.test.m <- as.data.frame(zscore.prot.exp.test.m)
#

#====selecting proteins belonging to selected pathways=

# x= zscore.prot.exp.learn.m ; matrix with protein expression
# selected.path = short.ReactomeDB.list ; 

gene.paths <- function(x,selected.path){

    gene.list.to.analysis <- character()
    p=1
    for (k in 1:length(selected.path)) {
      kk<-selected.path[[k]]
      for (r in 1:nrow(x)) {
         ifelse(rownames(x)[r] %in% kk,
                gene.list.to.analysis[p]<-rownames(x)[r],
                gene.list.to.analysis[p]<-"")
          p=p+1
      }
    }
    gene.list.to.analysis <- as.data.frame(gene.list.to.analysis)
    gene.list.to.analysis1 <- gene.list.to.analysis[!gene.list.to.analysis[,1] %in% c(""),]
    
    gene.list.to.analysis2<-as.data.frame(gene.list.to.analysis1[!duplicated(gene.list.to.analysis1)])
    colnames(gene.list.to.analysis2)<-"gene"
    
    x$gene <- rownames(x)
    
    x.T2<- plyr::join_all(list(gene.list.to.analysis2,x),by="gene", type = "inner")
    rownames(x.T2) <- x.T2$gene
    x.T2 <- dplyr::select(x.T2,-gene)
    
    return(x.T2)

}

zscore.gene.paths.learn <- gene.paths(x=zscore.prot.exp.learn.m, 
                                selected.path = short.ReactomeDB.list)

# zscore.gene.paths.test <- gene.paths(x=zscore.prot.exp.test.m, 
#                                  selected.path = short.ReactomeDB.list)

#=====Selecting from 'zscore.gene.paths' matrix only the 15 patients included in the B-raf paper (https://doi.org/10.3390/cancers11121981)

zscore.gene.paths.T <- as.data.frame(t(zscore.gene.paths.learn))
zscore.gene.paths.T$sample <- rownames(zscore.gene.paths.T)

zscore.gene.paths.T1 <- plyr::join_all(list(zscore.gene.paths.T,
                                           as.data.frame(annotation_colum[,c("sample","Braf_16pat")])),
                                      by="sample")

zscore.gene.paths.T1a <- subset(zscore.gene.paths.T1, zscore.gene.paths.T1$Braf_16pat %in% 1)
rownames(zscore.gene.paths.T1a) <- zscore.gene.paths.T1a$sample
zscore.gene.paths.T1a <- dplyr::select(zscore.gene.paths.T1a, -c("sample", "Braf_16pat"))

zscore.gene.paths.test <- as.data.frame(t(zscore.gene.paths.T1a))

##=============INGRID============================================================


# (1) gene expression z-scores in the form of a either data frame or matrix;
 
zscore.gene.paths.learn.1 <- t(zscore.gene.paths.learn)
#write.csv(zscore.gene.paths.learn.1, "zscore.gene.paths.learn.1.csv")

# (2) survival time and censoring indicator in the form of vectors;

t.learn = Braf.data.learn.list$t
d.learn = Braf.data.learn.list$d

# (3) pathway information in the form of a list, where each element is a vector of the names of gene belonging to the pathway.

Pathways.list <- short.ReactomeDB.list


#===
file.dir="output/output_ingrid/"
gene.exp=zscore.gene.paths.learn.1
path = Pathways.list
t=t.learn
d=d.learn
newx="not"  # yes or not
newx=newx  # yes or not

ingrid.funct <- function(gene.exp, path, t,d,file.dir,newx){
      
      # Gene Regrouping_small cohort
        geneRegroup.results <- geneRegroup(path)
      
   
      # # pre-filter. To refine the candidate set of genes. 
      #               It first eliminate the most unlikely gene predictors by conducting 
      #               a supervised pre-filtering using a Cox model applied
      #               to each gene separately and including only the genes with p-value < 0.5 for further analysis
      
      prefilter.results <- prefilter(data=gene.exp, 
                                     time=t, p.cut = 0.5, 
                                     status=d, plist=geneRegroup.results@gset)
      
      # Gene selection
      gene.results <- selectGene(prefilter.results, K = 1, etas = seq(0.1,0.9,0.1), 
                                 fold = 5, se1 = TRUE, method = "plik", par = FALSE,
                                 foldid = NULL, seed = 123)
      
      
      #The list of the SPLS regression coeficients of cancer-related genes can be generated using the function coef().
      
        # head(coef(gene.results)[[2]])
      
      coef.gene.results <- coef(gene.results)
      
      erer::write.list(z = coef.gene.results, 
                       file = paste(file.dir,"gene selected_in.csv",sep = ""))
       
      #==Pathway Selection
      
      path.results <- selectPath(gene.results)
      
      #=LASSO regression coeficients of cancer-related pathways can be generated using the function coef().
      head(coef(path.results))
      
      path.results.df <- coef(path.results)
      
      write.csv(path.results.df,paste(file.dir,"path.results.list.csv",sep = ""))
      
      png(filename = paste(file.dir,"HR_pathways.png",sep = ""), width = 500, height = 400, units = "px", pointsize = 16,
          bg = "white", res = NA, family = "", restoreConsole = TRUE,type = c("windows", "cairo", "cairo-png"))
      
      HR <- plot(path.results, type="HR")
      
      print(HR)
      dev.off()
      
      # Risk Group Prediction
      if(newx=="not"){
      predicted <- predict(path.results)
      erer::write.list(z = predicted, 
                       file = paste(file.dir,"predicted_in.csv",sep = ""))
      } else {
      predicted <- predict(path.results, newx=newx)
      erer::write.list(z = predicted, 
                       file = paste(file.dir,"predicted_in.csv",sep = ""))
      }
      
      list.patients <- as.data.frame(predicted$riskcat)

      write.csv(list.patients,paste(file.dir,"list of patients.csv",sep = ""))
      
      list.patients1 <- cbind(row.names(gene.exp),list.patients)
      
      write.csv(list.patients1,paste(file.dir,"list of patients1.csv",sep = ""))
      
      # Ploting results
      png(filename = paste(file.dir,"KM.png",sep = ""), width = 500, height = 400, units = "px", pointsize = 16,
          bg = "white", res = NA, family = "", restoreConsole = TRUE,type = c("windows", "cairo", "cairo-png"))
      
      #KM
      KM <- plot(path.results, type="KM")
      
      print(KM)
      dev.off()
      
      # Survival ROC
      png(filename = paste(file.dir,"ROC.png",sep = ""), width = 500, height = 400, units = "px", pointsize = 16,
          bg = "white", res = NA, family = "", restoreConsole = TRUE,type = c("windows", "cairo", "cairo-png"))
      ROC <- plot(path.results, type="ROC")
      
      print(ROC)
      
      dev.off()
      
     
      pdf(paste(file.dir,"KM.pdf",sep = "")) 
      
      #KM
      KM <- plot(path.results, type="KM")
      
      print(KM)
      dev.off()
      
      
      ingrid.results <- list(geneRegroup.results,
                             prefilter.results,
                             gene.results,
                             coef.gene.results,
                             path.results,
                             path.results.df,
                             predicted,
                             list.patients)
      names(ingrid.results) <- c("geneRegroup.results",
                                 "prefilter.results",
                                 "gene.results",
                                 "coef.gene.results",
                                 "path.results",
                                 "path.results.df",
                                 "predicted",
                                 "list.patients")
      return(ingrid.results)
}    

ingrid.results.learn <- ingrid.funct(file.dir="output/output_ingrid/",
                               gene.exp=zscore.gene.paths.learn.1,
                               path = Pathways.list,
                               t=t.learn,
                               d=d.learn,
                               newx="not")


# Testing

# (1) gene expression z-scores in the form of a either data frame or matrix;

zscore.gene.paths.test.1 <- t(zscore.gene.paths.test)

# (2) survival time and censoring indicator in the form of vectors;

t.test = Braf.data.test.list$t
d.test = Braf.data.test.list$d

# (3) pathway information in the form of a list, where each element is a vector of the names of gene belonging to the pathway.

ingrid.path <- ingrid.results.learn$path.results@pathSelected
names(ingrid.path) <- ingrid.results.learn$path.results@pathSelected

Pathways.list2 <- Pathways.list[ingrid.path]
Pathways.list2[[1]] <- ingrid.results.learn$gene.results@geneSelected$gene_set_263
names(Pathways.list2)[1] <- names(ingrid.results.learn$gene.results@geneSelected)[1]

Pathways.list2[[2]] <- ingrid.results.learn$gene.results@geneSelected$gene_set_264
names(Pathways.list2)[2] <- names(ingrid.results.learn$gene.results@geneSelected)[2]
  
# (4) InGRID

ingrid.results.test <- ingrid.funct(file.dir="output/output_ingrid/test/",
                                     gene.exp=zscore.gene.paths.test.1,
                                     path = Pathways.list2,
                                     t=t.test,
                                     d=d.test,
                                     newx="not")


#======
