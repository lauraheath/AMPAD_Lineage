#MAYO_CBE_TCX_logCPM.tsv (filtered log counts)
tcxCPMObj <- synapser::synGet('syn8466816')
Dat <- read.delim(tcxCPMObj$path,stringsAsFactors = F)
rownames(Dat) <- Dat$ensembl_gene_id

#Covariates (cell_metadata) file for all patients
tcxCovObj <- synapser::synGet('syn8466814')
Dat2 <- read.delim(tcxCovObj$path,stringsAsFactors = F)

#subsetting genes based on differential expression
de_file <- synapser::synGet('syn8468023')
de_file2 <- synapser::synGet('syn18475579')

de1 <- data.table::fread(de_file$path,data.table=F)
de2 <- data.table::fread(de_file2$path,data.table=F)
de3 <- dplyr::filter(de1,Model=='Diagnosis',Comparison=='AD-CONTROL',Tissue.ref=='TCX')
#AMP_mods <-  read.csv('Data/TCX_DE.csv')
#AMP_mods <- data.frame(GeneID=de3$ensembl_gene_id,logPV= - log(de3$adj.P.Val)/log(10), stringsAsFactors=F)
AMP_mods <- de2[,-1]
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

#Subsetting based on brain region
In_BR <- grep('TCX',Dat2$Tissue.Diagnosis)
DatNorm3 <- DatNorm2[,In_BR]
Dat3 <- Dat2[In_BR,]

#subsetting based on gender
Sex <- 'FEMALE'
In_S <- which(Dat3$Sex == Sex)
DatNorm4 <- DatNorm3[,In_S]
Dat4 <- Dat3[In_S,]


temp <- DatNorm4
temp2 <- Dat4

temp2$Diagnosis <- temp2$Tissue.SourceDiagnosis

temp2$Diagnosis[temp2$Tissue.SourceDiagnosis=='TCX.AD'] <- 'AD'
temp2$Diagnosis[temp2$Tissue.SourceDiagnosis=='TCX.CONTROL'] <- 'Control'
temp2$Diagnosis[temp2$Tissue.SourceDiagnosis=='TCX.PATH_AGE'] <- 'PA'
temp2$Diagnosis[temp2$Tissue.SourceDiagnosis=='TCX.PSP'] <- 'PSP'


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
GeneNamesAD <- as.data.frame(GeneNamesAD)
gene_short_name <- as.data.frame(gene_short_name)
gene_metadata <- cbind(GeneNamesAD, gene_short_name)
DatNorm5<-DatNorm4
rownames(DatNorm5) <- gene_metadata$GeneNamesAD
counts <- Matrix(as.matrix(DatNorm5), sparse=TRUE)
rownames(gene_metadata) <- gene_metadata$GeneNamesAD
rownames(temp2) <- temp2$SampleID

#to save the components of the monocle object for later use:
saveRDS(counts, file="TCX_countmatrix_Mono3.rds")
saveRDS(gene_metadata, file="~/AMPAD_Lineage/LH_reviewer_response/TCX_gene_metadata_Mono3.rds")
saveRDS(temp2, file="~/AMPAD_Lineage/LH_reviewer_response/TCX_cell_metadata_Mono3.rds")
#counts <- readRDS(file="~/AMPAD_Lineage/LH_reviewer_response/TCX_countmatrix_Mono3.rds")
#gene_metadata <- readRDS(file="~/AMPAD_Lineage/LH_reviewer_response/TCX_gene_metadata_Mono3.rds")
#temp2 <- readRDS(file="~/AMPAD_Lineage/LH_reviewer_response/TCX_cell_metadata_Mono3.rds")

#for monocle: expression matrix=counts, cell_metadata=temp2, gene_metadata=gene_metadata
# must first unload synapser because causes multiple definitions of S4Vectors
detach("package:synapser", unload=TRUE)
unloadNamespace("PythonEmbedInR") 
#run monocle and get Monocle Object: cds_tcx
library(monocle3)

cds_tcx <- new_cell_data_set(counts, cell_metadata=temp2, gene_metadata=gene_metadata)

cds_tcx <- preprocess_cds(cds_tcx, method="PCA", norm_method="none", num_dim=30)
plot_pc_variance_explained(cds_tcx)
#Monocle3 uses UMAP by default to reduce dimension (which is a manifold learning method)
cds_tcx <- reduce_dimension(cds_tcx, reduction_method="UMAP")
cds_tcx <- cluster_cells(cds_tcx, resolution=0.05)
cds_tcx <- learn_graph(cds_tcx)
p <- plot_cells(cds_tcx, 
                color_cells_by="Diagnosis", 
                cell_size=2) + theme(legend.text=element_text(size=6)) + theme(legend.position="right")
p

#choose the root node:
cds_tcx <- order_cells(cds_tcx)
p <- plot_cells(cds_tcx, 
                color_cells_by="Diagnosis", 
                cell_size=2) + theme(legend.text=element_text(size=6)) + theme(legend.position="right")
p

pstime_Mono3 <- pseudotime(cds_tcx, reduction_method="UMAP")
pstime_Mono3 <- as.data.frame(pstime_Mono3)
pstime_Mono3$SampleID <- rownames(pstime_Mono3)
write.csv(pstime_Mono3, file="~/AMPAD_Lineage/LH_reviewer_response/TCX_pstimes_Mono3.csv")

tiff(file='~/AMPAD_Lineage/LH_reviewer_response/TCX_celltraj_Mono3.tiff',height=85,width=100,units='mm',res=300)
p <- plot_cells(cds_tcx, 
                color_cells_by="Diagnosis", 
                cell_size=1.5,
                label_cell_groups=FALSE,
                label_leaves=FALSE,
                label_branch_points=FALSE) + theme(legend.text=element_text(size=6)) + theme(legend.position="right")
p
dev.off()


tcx <- dplyr::left_join(pstime_Mono3,temp2,by='SampleID')
tcxAD <- rep(NA,nrow(tcx))
tcxAD[tcx$Tissue.SourceDiagnosis == 'TCX.AD'] <- 1
tcxAD[tcx$Tissue.SourceDiagnosis == 'TCX.CONTROL'] <- 0
tcxDf <- data.frame(diagnosis=tcxAD,
                    pseudotime=tcx$pstime_Mono3,
                    stringsAsFactors = FALSE)
summary(glm(diagnosis ~ pseudotime,tcxDf,family='binomial'))
