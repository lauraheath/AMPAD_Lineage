#Loading data
synapser::synLogin()

#load rosmap filtered counts logCPM:
dlpfcCPMObj <- synapser::synGet('syn8456638')
Dat <- read.delim(dlpfcCPMObj$path,stringsAsFactors = F)

#synapse id of dat2 file (rosmap covariates): syn8466814
dlpfcCovObj <- synapser::synGet('syn11024258')
Dat2 <- read.delim(dlpfcCovObj$path,stringsAsFactors = F)

foobar <- synapser::synGet('syn11695124')
foobar2 <- data.table::fread(foobar$path,data.table=F)

#load ROSMAP_DLPFC_DiffExpression.tsv
de_file <- synapser::synGet('syn8456721')
de1 <- data.table::fread(de_file$path,data.table=F)
de3 <- dplyr::filter(de1,Model=='Diagnosis',Comparison=='AD-CONTROL')
AMP_mods <- data.frame(GeneID=de3$ensembl_gene_id,logPV= - log(de3$adj.P.Val)/log(10), stringsAsFactors=F)
AMP_mods <- dplyr::filter(AMP_mods,GeneID%in% foobar2$GeneID)
#subsetting genes based on differential expression 
#AMP_mods <-  read.csv('Data/DLPFC_DE.csv')
In <- which(AMP_mods$logPV >= 1)
AMP_mods <- AMP_mods[In,]


#Dat column names need to match the row names of gene_metadata (ensemble IDs)
#Dat 2 row names need to match Dat column names (patient IDs)
#Need to create gene_metadata file with row names=ensemble ids, matched to one column (gene_short_name)
#gene_metadata <- subset(de3, select=c("ensembl_gene_id", "hgnc_symbol"))

#Normalize all columns 
GeneNames <- Dat$ensembl_gene_id
GeneNamesAD <- AMP_mods$GeneID

Names <- colnames(Dat)

for (i in 1:length(Names)){
  
  Names[i] <- substring(Names[i],2)
  
}


colnames(Dat) <- Names
cNames <- Dat2$SampleID
l <- length(Names)

#deleting columns not in the covariate list
temp <- rep(T,l)
for (i in 1:l){
  if (!(Names[i] %in% cNames)){
    temp[i] <- F
  }
}

In <- which(temp)
#print(temp)
Dat <- Dat[,In]

#deleting extra rows in covariate list
Names <- Names[In]
l <- length(cNames)
temp <- rep(T,l)
for (i in 1:l){
  if (!(cNames[i] %in% Names)){
    temp[i] <- F
  }
}
In <- which(temp)
Dat2 <- Dat2[In,]

ColNorm <- function(Dat){
  
  M = max(colSums(Dat))
  l <- length(colnames(Dat))
  
  for( i in 1:l){
    
    Dat[,i] = Dat[,i]*(M/sum(Dat[,i]))
    
  }
  
  return(Dat)
}

DatNorm <- ColNorm(Dat)
In_genes <- which(GeneNames %in% GeneNamesAD)
DatNorm2 <- DatNorm[In_genes,]
GeneNamesAD <- GeneNames[In_genes]

#removing bad batches
DatNorm2 <- DatNorm2[,Dat2$Batch<7]
Dat2 <- Dat2[Dat2$Batch<7,] 



DatNorm3 <- DatNorm2
Dat3 <- Dat2

#Keeping only female data 
#Sex <- 'FEMALE'
In_S <- which(Dat3$msex == 0)
DatNorm4 <- DatNorm3[,In_S]
Dat4 <- Dat3[In_S,]


temp <- DatNorm4
temp2 <- Dat4

#converting ENSG to gene symbols
#convert to gene symbol
convertEnsemblToHgnc <- function(ensemblIds){
  
  ensembl=biomaRt::useMart('ENSEMBL_MART_ENSEMBL',
                           dataset = 'hsapiens_gene_ensembl',
                           host='useast.ensembl.org')
  
  genes<-getBM(attributes = c('ensembl_gene_id','external_gene_name'),
               filters='ensembl_gene_id',
               values=ensemblIds,
               mart=ensembl)
  return(genes)
}
Make.Gene.Symb <- function(GeneENSG){
  
  #source('convertEnsemblToHgnc.R')
  GeneConv <- convertEnsemblToHgnc(GeneENSG)
  Symb <- as.character(c(1:length(GeneENSG)))
  
  for (i in 1:length(GeneENSG)){
    In <- which(GeneConv$ensembl_gene_id == GeneENSG[i])
    if (length(In)>0){
      Symb[i] <- GeneConv$external_gene_name[In]
    }
  }
  
  return(Symb)
  
}
gene_short_name <- Make.Gene.Symb(GeneNamesAD)

rosmapObj <- synapser::synGet('syn3191087')
rosmap <- data.table::fread(rosmapObj$path,data.table=F)

#add in braak score & cerad score
rosmapIdObj <- synapser::synGet('syn3382527')
rosmapId <- data.table::fread(rosmapIdObj$path,data.table=F)
rosmapId <- dplyr::select(rosmapId,projid,rnaseq_id)
rosmapRNAid<-dplyr::left_join(rosmapId,rosmap)
#remove duplicate rows
rosmapRNAid <- unique(rosmapRNAid)
rosmapRNAid2 <- subset(rosmapRNAid, select=c(rnaseq_id,braaksc,ceradsc))
names(rosmapRNAid2)[names(rosmapRNAid2) == "rnaseq_id"] <- "SampleID"

temp2<-dplyr::left_join(temp2,rosmapRNAid2, by="SampleID")

temp2$braaksc <- factor(temp2$braaksc,levels = c(0:6))
temp2$ceradsc <- factor(temp2$ceradsc,levels = c(1:4))
temp2$cogdx <- factor(temp2$cogdx,levels = c(1:6))
temp2$cogdxNew <- ifelse(temp2$cogdx==1, 'NCI',
                         ifelse(temp2$cogdx==2, 'MCI',
                                ifelse(temp2$cogdx==4, 'LOAD', NA)))
temp2$cogdxNew<- factor(temp2$cogdxNew,levels = c('NCI','MCI','LOAD'))


GeneNamesAD <- as.data.frame(GeneNamesAD)
gene_short_name <- as.data.frame(gene_short_name)
gene_metadata <- cbind(GeneNamesAD, gene_short_name)
DatNorm5<-DatNorm4
rownames(DatNorm5) <- gene_metadata$GeneNamesAD
counts <- Matrix(as.matrix(DatNorm5), sparse=TRUE)
rownames(gene_metadata) <- gene_metadata$GeneNamesAD
rownames(temp2) <- temp2$SampleID

#to save the components of the monocle object for later use:
saveRDS(counts, file="~/AMPAD_Lineage/LH_reviewer_response/DLPFC_countmatrix_Mono3.rds")
saveRDS(gene_metadata, file="~/AMPAD_Lineage/LH_reviewer_response/DLPFC_gene_metadata_Mono3.rds")
saveRDS(temp2, file="~/AMPAD_Lineage/LH_reviewer_response/DLPFC_cell_metadata_Mono3.rds")

counts <- readRDS(file="~/AMPAD_Lineage/LH_reviewer_response/DLPFC_processed_countmatrix_Mono3.rds")
gene_metadata <- readRDS(file="~/AMPAD_Lineage/LH_reviewer_response/DLPFC_processed_gene_metadata_Mono3.rds")
temp2 <- readRDS(file="~/AMPAD_Lineage/LH_reviewer_response/DLPFC_processed_cell_metadata_Mono3.rds")


#for monocle: expression matrix=counts, cell_metadata=temp2, gene_metadata=gene_metadata

# must first unload synapser because causes multiple definitions of S4Vectors
detach("package:synapser", unload=TRUE)
unloadNamespace("PythonEmbedInR") 
#run monocle and get Monocle Object: cds_tcx
library(monocle3)

cds_dlpfc <- new_cell_data_set(counts, cell_metadata=temp2, gene_metadata=gene_metadata)

cds_dlpfc <- preprocess_cds(cds_dlpfc, method="PCA", norm_method="none", num_dim=30)
#plot_pc_variance_explained(cds_dlpfc)
#Monocle3 uses UMAP by default to reduce dimension (which is a manifold learning method)
cds_dlpfc <- reduce_dimension(cds_dlpfc, reduction_method="UMAP")
cds_dlpfc <- cluster_cells(cds_dlpfc, resolution=0.3)
cds_dlpfc <- learn_graph(cds_dlpfc)
p <- plot_cells(cds_dlpfc, 
                color_cells_by="Diagnosis", 
                cell_size=2) + theme(legend.text=element_text(size=6)) + theme(legend.position="right")
p


#a shiny app interface will allow you to select a root node; i selected the starting node in the cluster with fewer AD cases
cds_dlpfc <- order_cells(cds_dlpfc, reduction_method="UMAP") 
plot_cells(cds_dlpfc, 
           color_cells_by="pseudotime", 
           cell_size=2)

pstime_Mono3 <- pseudotime(cds_dlpfc, reduction_method="UMAP")
pstime_Mono3 <- as.data.frame(pstime_Mono3)
pstime_Mono3$SampleID <- rownames(pstime_Mono3)
write.csv(pstime_Mono3, file="~/AMPAD_Lineage/LH_reviewer_response/DLPFC_pstimes_Mono3.csv")
pstime_Mono3<-read.csv(file="~/AMPAD_Lineage/LH_reviewer_response/DLPFC_pstimes_Mono3.csv")
pstime_Mono3 <- pstime_Mono3[,-c(1)]

tiff(file='~/AMPAD_Lineage/LH_reviewer_response/DLPFC_celltraj_Mono3.tiff',height=85,width=100,units='mm',res=300)
p <- plot_cells(cds_dlpfc, 
                color_cells_by="Diagnosis", 
                cell_size=1.5,
                label_cell_groups=FALSE,
                label_leaves=FALSE,
                label_branch_points=FALSE) + theme(legend.text=element_text(size=6)) + theme(legend.position="right")
p
dev.off()

dlpfc <- dplyr::left_join(temp2, pstime_Mono3, by='SampleID')
#braak score plots:
tiff(file='~/AMPAD_Lineage/LH_reviewer_response/DLPFC_celltraj_braak_mono3.tiff',height=85,width=100,units='mm',res=300)
g<- plot_cells(cds_dlpfc, 
               color_cells_by="braaksc", 
               cell_size=1.5,
               label_cell_groups=FALSE,
               label_leaves=FALSE,
               label_branch_points=FALSE) + theme(legend.text=element_text(size=6)) + theme(legend.position="right")
#g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="Braak Score")
g
dev.off()
tiff(file='~/AMPAD_Lineage/LH_reviewer_response/DLPFC_boxplot_braak_mono3.tiff',height=85,width=100,units='mm',res=300)
g <- ggplot2::ggplot(dlpfc, aes(x=braaksc, y=scale(pstime_Mono3,center=F),fill=braaksc)) 
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::stat_summary(fun=mean, geom="point", shape=23, size=2)
g <- g + ggplot2::theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
                        legend.text=element_text(size=15)) 
#g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(fill="Braak\nScore",y="Pseudotime",x="Braak Score")
g
dev.off()

###CERAD score plots
tiff(file='~/AMPAD_Lineage/LH_reviewer_response/DLPFC_celltraj_cerad_mono3.tiff',height=85,width=100,units='mm',res=300)
g<- plot_cells(cds_dlpfc, 
               color_cells_by="ceradsc", 
               cell_size=1.5,
               label_cell_groups=FALSE,
               label_leaves=FALSE,
               label_branch_points=FALSE) + theme(legend.text=element_text(size=6)) + theme(legend.position="right")
#g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="CERAD Score")
g
dev.off()
tiff(file='~/AMPAD_Lineage/LH_reviewer_response/DLPFC_boxplot_cerad_mono3.tiff',height=85,width=100,units='mm',res=300)
g <- ggplot2::ggplot(dlpfc, aes(x=ceradsc, y=scale(pstime_Mono3,center=F),fill=ceradsc)) 
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::stat_summary(fun=mean, geom="point", shape=23, size=2)
g <- g + ggplot2::theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
                        legend.text=element_text(size=15)) 
#g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(fill="CERAD\nScore",y="Pseudotime",x="CERAD score")
g
dev.off()


###COGDX plots

tiff(file='~/AMPAD_Lineage/LH_reviewer_response/DLPFC_celltraj_cogdx_mono3.tiff',height=85,width=100,units='mm',res=300)
g<- plot_cells(cds_dlpfc, 
               color_cells_by="cogdx", 
               cell_size=1.5,
               label_cell_groups=FALSE,
               label_leaves=FALSE,
               label_branch_points=FALSE) + theme(legend.text=element_text(size=6)) + theme(legend.position="right")
#g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="CogDx")
g
dev.off()
dlpfc$cogdx <- as.factor(dlpfc$cogdx)
tiff(file='~/AMPAD_Lineage/LH_reviewer_response/DLPFC_boxplot_cogdx_mono3.tiff',height=85,width=100,units='mm',res=300)
g <- ggplot2::ggplot(dlpfc, aes(x=cogdx, y=scale(pstime_Mono3,center=F),fill=cogdx)) 
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::stat_summary(fun=mean, geom="point", shape=23, size=2)
g <- g + ggplot2::theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
                        legend.text=element_text(size=15)) 
#g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(fill="CogDx\nScore",y="Pseudotime",x="CogDx")
g
dev.off()






