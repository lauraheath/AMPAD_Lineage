#compiling source data for synapse


#figure S6a-6d: put adjusted pseudotimes (pmi, rin, pcs, all), plus sampleid, diagnosis
tcxLineageTimesLH <- read.csv("~/AMPAD_Lineage/LH_reviewer_response/TCX_F_PStimes_adjusted.csv", stringsAsFactors = FALSE)
dlpfcLineageTimesLH <- read.csv("~/AMPAD_Lineage/LH_reviewer_response/DLPFC_F_PStimes_adjusted.csv", stringsAsFactors = FALSE)
tcxCovObj <- synapser::synGet("syn8466814")
tcxCov <- data.table::fread(tcxCovObj$path,data.table=F)
dlpfcCovObj <- synapser::synGet("syn11024258")
dlpfcCov <- data.table::fread(dlpfcCovObj$path,data.table=F)
tcxLineageTimesLH <- tcxLineageTimesLH[,-c(1)]
dlpfcLineageTimesLH <- dlpfcLineageTimesLH[,-c(1)]
dlpfc <- dplyr::left_join(dlpfcLineageTimesLH,dlpfcCov,by='SampleID')
tcx <- dplyr::left_join(tcxLineageTimesLH,tcxCov,by='SampleID')

tcx$diagnosis<-tcx$Tissue.SourceDiagnosis
tcx$diagnosis[tcx$diagnosis=='TCX.AD'] <- 'AD'
tcx$diagnosis[tcx$diagnosis=='TCX.CONTROL'] <- 'Control'
tcx$diagnosis[tcx$diagnosis=='TCX.PSP'] <- 'PSP'
tcx$diagnosis[tcx$diagnosis=='TCX.PATH_AGE'] <- 'PA'
tcx$Study<-'MayoRNAseq'
tcxDf <- data.frame(SampleID=tcx$SampleID,
                    Pseudotime_rin_adjusted=tcx$Pseudotime_RIN,
                    Pseudotime_pmi_adjusted=tcx$Pseudotime_PMI,
                    Pseudotime_pc_adjusted=tcx$Pseudotime_PCs,
                    Pseudotime_all_adjusted=tcx$Pseudotime_ALL,
                    Diagnosis=tcx$diagnosis,
                    Study=tcx$Study,
                    stringsAsFactors = FALSE)
tcxDf <- subset(tcxDf, tcxDf$Diagnosis=='AD'|tcxDf$Diagnosis=='Control')


dlpfc$diagnosis<-dlpfc$Diagnosis
dlpfc$diagnosis[dlpfc$diagnosis=='AD'] <- 'AD'
dlpfc$diagnosis[dlpfc$diagnosis=='CONTROL'] <- 'Control'
dlpfc$diagnosis[dlpfc$diagnosis=='OTHER'] <- 'Other'
dlpfc$Study<-'ROSMAP'
dlpfcDf <- data.frame(SampleID=dlpfc$SampleID,
                    Pseudotime_rin_adjusted=dlpfc$Pseudotime_RIN,
                    Pseudotime_pmi_adjusted=dlpfc$Pseudotime_PMI,
                    Pseudotime_pc_adjusted=dlpfc$Pseudotime_PCs,
                    Pseudotime_all_adjusted=dlpfc$Pseudotime_ALL,
                    Diagnosis=dlpfc$diagnosis,
                    Study=dlpfc$Study,
                    stringsAsFactors = FALSE)
dlpfcDf <- subset(dlpfcDf, dlpfcDf$Diagnosis=='AD'|dlpfcDf$Diagnosis=='Control')

#join the two data frames and output with figure title
figureS6 <- rbind(dlpfcDf, tcxDf)

g <- ggplot2::ggplot(figureS6,ggplot2::aes(x=Study,
                                             y=Pseudotime_rin_adjusted,
                                             color=Diagnosis))
g <- g + ggplot2::geom_boxplot() + theme(text = element_text(size = 22)) 
g <- g + ggplot2::geom_point(size=2.5, position=ggplot2::position_jitterdodge())
g <- g + ggplot2::scale_color_manual(values=viridis::viridis(3)[1:2])
g

g <- ggplot2::ggplot(figureS6,ggplot2::aes(x=Study,
                                           y=Pseudotime_pmi_adjusted,
                                           color=Diagnosis))
g <- g + ggplot2::geom_boxplot() + theme(text = element_text(size = 22)) 
g <- g + ggplot2::geom_point(size=2.5, position=ggplot2::position_jitterdodge())
g <- g + ggplot2::scale_color_manual(values=viridis::viridis(3)[1:2])
g

g <- ggplot2::ggplot(figureS6,ggplot2::aes(x=Study,
                                           y=Pseudotime_pc_adjusted,
                                           color=Diagnosis))
g <- g + ggplot2::geom_boxplot() + theme(text = element_text(size = 22)) 
g <- g + ggplot2::geom_point(size=2.5, position=ggplot2::position_jitterdodge())
g <- g + ggplot2::scale_color_manual(values=viridis::viridis(3)[1:2])
g

g <- ggplot2::ggplot(figureS6,ggplot2::aes(x=Study,
                                           y=Pseudotime_all_adjusted,
                                           color=Diagnosis))
g <- g + ggplot2::geom_boxplot() + theme(text = element_text(size = 22)) 
g <- g + ggplot2::geom_point(size=2.5, position=ggplot2::position_jitterdodge())
g <- g + ggplot2::scale_color_manual(values=viridis::viridis(3)[1:2])
g

write.csv(figureS6, file="~/AMPAD_Lineage/LH_reviewer_response/SourceData_Synapse/figureS6_data.csv")

#figure S6A is for pmi:
figureS6A <- subset(figureS6, select=c(SampleID, Pseudotime_pmi_adjusted, Diagnosis, Study))
write.csv(figureS6A, file="~/AMPAD_Lineage/LH_reviewer_response/SourceData_Synapse/figureS6A_data.csv", row.names=FALSE)
file <- File(path='~/AMPAD_Lineage/LH_reviewer_response/SourceData_Synapse/figureS6A_data.csv', parent='syn23262420')
file <- synapser::synStore(file)
#doublecheck the file
s6check <- read.csv(file='~/AMPAD_Lineage/LH_reviewer_response/SourceData_Synapse/figureS6A_data.csv')

#figure S6B is for principal components:
figureS6B <- subset(figureS6, select=c(SampleID, Pseudotime_pc_adjusted, Diagnosis, Study))
write.csv(figureS6B, file="~/AMPAD_Lineage/LH_reviewer_response/SourceData_Synapse/figureS6B_data.csv", row.names=FALSE)
file <- synapser::File(path='~/AMPAD_Lineage/LH_reviewer_response/SourceData_Synapse/figureS6B_data.csv', parent='syn23262420')
file <- synapser::synStore(file)

#figure S6C is for rin:
figureS6C <- subset(figureS6, select=c(SampleID, Pseudotime_rin_adjusted, Diagnosis, Study))
write.csv(figureS6C, file="~/AMPAD_Lineage/LH_reviewer_response/SourceData_Synapse/figureS6C_data.csv", row.names=FALSE)
file <- synapser::File(path='~/AMPAD_Lineage/LH_reviewer_response/SourceData_Synapse/figureS6C_data.csv', parent='syn23262420')
file <- synapser::synStore(file)

#figure S6D is for all:
figureS6D <- subset(figureS6, select=c(SampleID, Pseudotime_all_adjusted, Diagnosis, Study))
write.csv(figureS6D, file="~/AMPAD_Lineage/LH_reviewer_response/SourceData_Synapse/figureS6D_data.csv", row.names=FALSE)
file <- synapser::File(path='~/AMPAD_Lineage/LH_reviewer_response/SourceData_Synapse/figureS6D_data.csv', parent='syn23262420')
file <- synapser::synStore(file)



####figureS12A-D: GWAS LOAD distributions
# tcxLineageTimesLH <- read.csv("~/AMPAD_Lineage/LH_reviewer_response/TCX_F_PStimes_adjusted.csv", stringsAsFactors = FALSE)
# dlpfcLineageTimesLH <- read.csv("~/AMPAD_Lineage/LH_reviewer_response/DLPFC_F_PStimes_adjusted.csv", stringsAsFactors = FALSE)
# tcxCovObj <- synapser::synGet("syn8466814")
# tcxCov <- data.table::fread(tcxCovObj$path,data.table=F)
# dlpfcCovObj <- synapser::synGet("syn11024258")
# dlpfcCov <- data.table::fread(dlpfcCovObj$path,data.table=F)
# tcxLineageTimesLH <- tcxLineageTimesLH[,-c(1)]
# dlpfcLineageTimesLH <- dlpfcLineageTimesLH[,-c(1)]
# dlpfc <- dplyr::left_join(dlpfcLineageTimesLH,dlpfcCov,by='SampleID')
# tcx <- dplyr::left_join(tcxLineageTimesLH,tcxCov,by='SampleID')



#needed to tweak the convert function to a different host (biomaRt stopped working)
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

ad_gwas <- c("CR1",
             "BIN1",
             "INPP5D",
             "HLA-DRB1",
             "TREM2",
             "MEF2C",
             "NME8",
             "CD2AP",
             "NYAP1",
             "EPHA1",
             "PTK2B",
             "CLU",
             "SPI1",
             "MS4A2",
             "PICALM",
             "SORL1",
             "FERMT2",
             "SLC24A4",
             "ABCA7",
             "APOE",
             "CASS4",
             "ECHDC3",
             "ACE",
             "NDUFAF6",
             "ECHDC3",
             "ADAMTS20",
             "SPPL2A",
             "ADAM10",
             "IQCK",
             "MIR142",
             "ACE",
             "ADAMTS1",
             "SUCLG2P4",
             "FST",
             "OARD1",
             "WWOX",
             "MAF",
             "CD55",
             "YOD1",
             "HLA-DRB1",
             "PSMB8",
             "C4A",
             "GPSM3",
             "HLA-DPA1",
             "HLA-DQA1",
             "HLA-DRA",
             "HLA-DRB5",
             "PSMB9",
             "CD2AP",
             "AGFG2",
             "PILRA",
             "EPHB4",
             "C7orf43",
             "GAL3ST4",
             "ZKSCCAN1",
             "FAM131B",
             "PSMC3",
             "ACP2",
             "C1QTNF4",
             "CELF1",
             "MTCH2",
             "NDUFS3",
             "NUP160",
             "MS4A6A",
             "MS4A7",
             "MS4A4A",
             "EED",
             "PICALM",
             "STYX",
             "RIN3",
             "HMHA1",
             "CNN2",
             "WDR18",
             "CASS4")


###### start with pmi-adjusted pseudotimes #######

tcxCPMObj <- synapser::synGet('syn8466816')
#Dat <- read.delim(tcxCPMObj$path,stringsAsFactors = F)
Dat_tcx <- data.table::fread(tcxCPMObj$path,data.table=F)
sampleIds <- colnames(Dat_tcx)[-1]
geneIds <- Dat_tcx$ensembl_gene_id
Dat_tcx <- Dat_tcx[,-1]
Dat_tcx <- t(Dat_tcx)
colnames(Dat_tcx) <- geneIds
Dat_tcx <- data.frame(Dat_tcx,stringsAsFactors=F)
Dat_tcx$sampleId <- sampleIds
tcx <- dplyr::left_join(tcx,Dat_tcx,by=c('SampleID'='sampleId'))
tcx2 <- dplyr::select(tcx,dplyr::starts_with("ENSG"))

corvec <- cor(tcx2,tcx$Pseudotime_PMI,method='spearman')
corDf <- data.frame(geneid=colnames(tcx2),cor=corvec,stringsAsFactors=F)
corDf2 <- corDf
corDf2$external_gene_name <- Make.Gene.Symb(corDf2$geneid)
corDf2$cor <- NULL
corDf2 <- dplyr::left_join(corDf,corDf2,by=c('geneid'))

dlpfcCPMObj <- synapser::synGet('syn8456638')
Dat <- data.table::fread(dlpfcCPMObj$path,data.table=F)
sampleIds <- colnames(Dat)[-1]
geneIds <- Dat$ensembl_gene_id
Dat <- Dat[,-1]
Dat <- t(Dat)
colnames(Dat) <- geneIds
Dat <- data.frame(Dat,stringsAsFactors=F)
Dat$sampleId <- sampleIds
dlpfc <- dplyr::left_join(dlpfc,Dat,by=c('SampleID'='sampleId'))
dlpfc2 <- dplyr::select(dlpfc,dplyr::starts_with("ENSG"))

corvec <- cor(dlpfc2,dlpfc$Pseudotime_PMI,method='spearman')
corDfdlpfc <- data.frame(geneid=colnames(dlpfc2),cor=corvec,stringsAsFactors=F)
corDfdlpfc2 <- corDfdlpfc
corDfdlpfc2$external_gene_name <- Make.Gene.Symb(corDfdlpfc2$geneid)
corDfdlpfc2$cor <- NULL
corDfdlpfc2 <- dplyr::left_join(corDfdlpfc,corDfdlpfc2,by=c('geneid'))


corDf <- corDf2
mean(abs(corDf$cor))
mean(abs(corDf[corDf$external_gene_name %in% ad_gwas,]$cor))

corDfdlpfc <- corDfdlpfc2
mean(abs(corDfdlpfc$cor))
mean(abs(corDfdlpfc[corDfdlpfc$external_gene_name %in% ad_gwas,]$cor))

corDf$adGwas <- corDf$external_gene_name %in% ad_gwas
corDfdlpfc$adGwas <- corDfdlpfc$external_gene_name %in% ad_gwas
corDf$brainRegion <- 'TCX'
corDfdlpfc$brainRegion <- 'DLPFC'
corDfcombined <- rbind(corDf,corDfdlpfc)
colnames(corDfcombined)[4] <- 'LOADGWASGene'
g <- ggplot2::ggplot(corDfcombined,ggplot2::aes(x=brainRegion,y=cor,fill=LOADGWASGene))
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(x = 'Brain Region',y='Correlation with PMI-adj pseudotime',fill='LOAD\nGWAS\nGene')
g

#figure 12A is for PMI-adjusted pseudotimes:
figureS12A<-corDfcombined
write.csv(figureS12A, file="~/AMPAD_Lineage/LH_reviewer_response/SourceData_Synapse/figureS12A_data.csv", row.names=FALSE)
file <- synapser::File(path='~/AMPAD_Lineage/LH_reviewer_response/SourceData_Synapse/figureS12A_data.csv', parent='syn23262420')
file <- synapser::synStore(file)



##### repeat from the correlation calcs for PC-adjusted pseudotimes ######
corvec <- cor(tcx2,tcx$Pseudotime_PCs,method='spearman')
corDf <- data.frame(geneid=colnames(tcx2),cor=corvec,stringsAsFactors=F)
corDf2 <- corDf
corDf2$external_gene_name <- Make.Gene.Symb(corDf2$geneid)
corDf2$cor <- NULL
corDf2 <- dplyr::left_join(corDf,corDf2,by=c('geneid'))

#to look at PC-adjusted pseudotime in dlpfc, need to remove NAs to run the correlation in the dlpfc data
dlpfc_pc <- subset(dlpfc, !is.na(dlpfc$Pseudotime_PCs))
dlpfc2_pc <- dplyr::select(dlpfc_pc,dplyr::starts_with("ENSG"))
corvec <- cor(dlpfc2_pc,dlpfc_pc$Pseudotime_PCs,method='spearman')
corDfdlpfc <- data.frame(geneid=colnames(dlpfc2),cor=corvec,stringsAsFactors=F)
corDfdlpfc2 <- corDfdlpfc
corDfdlpfc2$external_gene_name <- Make.Gene.Symb(corDfdlpfc2$geneid)
corDfdlpfc2$cor <- NULL
corDfdlpfc2 <- dplyr::left_join(corDfdlpfc,corDfdlpfc2,by=c('geneid'))

corDf <- corDf2
mean(abs(corDf$cor))
mean(abs(corDf[corDf$external_gene_name %in% ad_gwas,]$cor))

corDfdlpfc <- corDfdlpfc2
mean(abs(corDfdlpfc$cor))
mean(abs(corDfdlpfc[corDfdlpfc$external_gene_name %in% ad_gwas,]$cor))

corDf$adGwas <- corDf$external_gene_name %in% ad_gwas
corDfdlpfc$adGwas <- corDfdlpfc$external_gene_name %in% ad_gwas
corDf$brainRegion <- 'TCX'
corDfdlpfc$brainRegion <- 'DLPFC'
corDfcombined <- rbind(corDf,corDfdlpfc)
colnames(corDfcombined)[4] <- 'LOADGWASGene'
g <- ggplot2::ggplot(corDfcombined,ggplot2::aes(x=brainRegion,y=cor,fill=LOADGWASGene))
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(x = 'Brain Region',y='Correlation with PC-adj pseudotime',fill='LOAD\nGWAS\nGene')
g

#figure 12B is for PC-adjusted pseudotimes:
figureS12B<-corDfcombined
write.csv(figureS12B, file="~/AMPAD_Lineage/LH_reviewer_response/SourceData_Synapse/figureS12B_data.csv", row.names=FALSE)
file <- synapser::File(path='~/AMPAD_Lineage/LH_reviewer_response/SourceData_Synapse/figureS12B_data.csv', parent='syn23262420')
file <- synapser::synStore(file)



##### repeat from the correlation calcs for RIN-adjusted pseudotimes ######
corvec <- cor(tcx2,tcx$Pseudotime_RIN,method='spearman')
corDf <- data.frame(geneid=colnames(tcx2),cor=corvec,stringsAsFactors=F)
corDf2 <- corDf
corDf2$external_gene_name <- Make.Gene.Symb(corDf2$geneid)
corDf2$cor <- NULL
corDf2 <- dplyr::left_join(corDf,corDf2,by=c('geneid'))

corvec <- cor(dlpfc2,dlpfc$Pseudotime_RIN,method='spearman')
corDfdlpfc <- data.frame(geneid=colnames(dlpfc2),cor=corvec,stringsAsFactors=F)
corDfdlpfc2 <- corDfdlpfc
corDfdlpfc2$external_gene_name <- Make.Gene.Symb(corDfdlpfc2$geneid)
corDfdlpfc2$cor <- NULL
corDfdlpfc2 <- dplyr::left_join(corDfdlpfc,corDfdlpfc2,by=c('geneid'))

corDf <- corDf2
mean(abs(corDf$cor))
mean(abs(corDf[corDf$external_gene_name %in% ad_gwas,]$cor))

corDfdlpfc <- corDfdlpfc2
mean(abs(corDfdlpfc$cor))
mean(abs(corDfdlpfc[corDfdlpfc$external_gene_name %in% ad_gwas,]$cor))

corDf$adGwas <- corDf$external_gene_name %in% ad_gwas
corDfdlpfc$adGwas <- corDfdlpfc$external_gene_name %in% ad_gwas
corDf$brainRegion <- 'TCX'
corDfdlpfc$brainRegion <- 'DLPFC'
corDfcombined <- rbind(corDf,corDfdlpfc)
colnames(corDfcombined)[4] <- 'LOADGWASGene'
g <- ggplot2::ggplot(corDfcombined,ggplot2::aes(x=brainRegion,y=cor,fill=LOADGWASGene))
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(x = 'Brain Region',y='Correlation with RIN-adj pseudotime',fill='LOAD\nGWAS\nGene')
g

#figure 12C is for RIN-adjusted pseudotimes:
figureS12C<-corDfcombined
write.csv(figureS12C, file="~/AMPAD_Lineage/LH_reviewer_response/SourceData_Synapse/figureS12C_data.csv", row.names=FALSE)
file <- synapser::File(path='~/AMPAD_Lineage/LH_reviewer_response/SourceData_Synapse/figureS12C_data.csv', parent='syn23262420')
file <- synapser::synStore(file)


##### repeat from the correlation calcs for ALL-adjusted pseudotimes ######
corvec <- cor(tcx2,tcx$Pseudotime_ALL,method='spearman')
corDf <- data.frame(geneid=colnames(tcx2),cor=corvec,stringsAsFactors=F)
corDf2 <- corDf
corDf2$external_gene_name <- Make.Gene.Symb(corDf2$geneid)
corDf2$cor <- NULL
corDf2 <- dplyr::left_join(corDf,corDf2,by=c('geneid'))

#to look at ALL-adjusted pseudotime in dlpfc, need to remove NAs to run the correlation in the dlpfc data
dlpfc_pc <- subset(dlpfc, !is.na(dlpfc$Pseudotime_ALL))
dlpfc2_pc <- dplyr::select(dlpfc_pc,dplyr::starts_with("ENSG"))
corvec <- cor(dlpfc2_pc,dlpfc_pc$Pseudotime_ALL,method='spearman')
corDfdlpfc <- data.frame(geneid=colnames(dlpfc2),cor=corvec,stringsAsFactors=F)
corDfdlpfc2 <- corDfdlpfc
corDfdlpfc2$external_gene_name <- Make.Gene.Symb(corDfdlpfc2$geneid)
corDfdlpfc2$cor <- NULL
corDfdlpfc2 <- dplyr::left_join(corDfdlpfc,corDfdlpfc2,by=c('geneid'))

corDf <- corDf2
mean(abs(corDf$cor))
mean(abs(corDf[corDf$external_gene_name %in% ad_gwas,]$cor))

corDfdlpfc <- corDfdlpfc2
mean(abs(corDfdlpfc$cor))
mean(abs(corDfdlpfc[corDfdlpfc$external_gene_name %in% ad_gwas,]$cor))

corDf$adGwas <- corDf$external_gene_name %in% ad_gwas
corDfdlpfc$adGwas <- corDfdlpfc$external_gene_name %in% ad_gwas
corDf$brainRegion <- 'TCX'
corDfdlpfc$brainRegion <- 'DLPFC'
corDfcombined <- rbind(corDf,corDfdlpfc)
colnames(corDfcombined)[4] <- 'LOADGWASGene'
g <- ggplot2::ggplot(corDfcombined,ggplot2::aes(x=brainRegion,y=cor,fill=LOADGWASGene))
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(x = 'Brain Region',y='Correlation with ALL-adj pseudotime',fill='LOAD\nGWAS\nGene')
g

#figure 12D is for ALL-adjusted pseudotimes:
figureS12D<-corDfcombined
write.csv(figureS12D, file="~/AMPAD_Lineage/LH_reviewer_response/SourceData_Synapse/figureS12D_data.csv", row.names=FALSE)
file <- synapser::File(path='~/AMPAD_Lineage/LH_reviewer_response/SourceData_Synapse/figureS12D_data.csv', parent='syn23262420')
file <- synapser::synStore(file)



#upload pseudotimes for rosmap males only and males/females combined into synapse, after running the Monocle Function code
MalePseudotimes <- read.csv(file='~/AMPAD_Lineage/LH_reviewer_response/Monly_MF_analyses/DLPFC_M_PStimes.csv')
write.csv(MalePseudotimes, file="~/AMPAD_Lineage/LH_reviewer_response/SourceData_Synapse/figureS8B_Male_pseudotimes_rosmap.csv", row.names=FALSE)
file <- synapser::File(path='~/AMPAD_Lineage/LH_reviewer_response/SourceData_Synapse/figureS8B_Male_pseudotimes_rosmap.csv', parent='syn23262420')
file <- synapser::synStore(file)

ALLPseudotimes <- read.csv(file='~/AMPAD_Lineage/LH_reviewer_response/Monly_MF_analyses/DLPFC_MF_PStimes.csv')
write.csv(ALLPseudotimes, file="~/AMPAD_Lineage/LH_reviewer_response/SourceData_Synapse/figureS10B_ALL_pseudotimes_rosmap.csv", row.names=FALSE)
file <- synapser::File(path='~/AMPAD_Lineage/LH_reviewer_response/SourceData_Synapse/figureS10B_ALL_pseudotimes_rosmap.csv', parent='syn23262420')
file <- synapser::synStore(file)

#upload pseudotimes for mayoRNAseq males only and males/females combined into synapse, after running the Monocle Function code
MalePseudotimes <- read.csv(file='~/AMPAD_Lineage/LH_reviewer_response/Monly_MF_analyses/TCX_M_PStimes.csv')
write.csv(MalePseudotimes, file="~/AMPAD_Lineage/LH_reviewer_response/SourceData_Synapse/figureS8A_Male_pseudotimes_mayoRNAseq.csv", row.names=FALSE)
file <- synapser::File(path='~/AMPAD_Lineage/LH_reviewer_response/SourceData_Synapse/figureS8A_Male_pseudotimes_mayoRNAseq.csv', parent='syn23262420')
file <- synapser::synStore(file)

ALLPseudotimes <- read.csv(file='~/AMPAD_Lineage/LH_reviewer_response/Monly_MF_analyses/TCX_MF_PStimes.csv')
ALLPseudotimes <- ALLPseudotimes[,-c(1)]
write.csv(ALLPseudotimes, file="~/AMPAD_Lineage/LH_reviewer_response/SourceData_Synapse/figureS10A_ALL_pseudotimes_mayoRNAseq.csv", row.names=FALSE)
file <- synapser::File(path='~/AMPAD_Lineage/LH_reviewer_response/SourceData_Synapse/figureS10A_ALL_pseudotimes_mayoRNAseq.csv', parent='syn23262420')
file <- synapser::synStore(file)




#Upload pseudotimes calculated in males only
tcxLineageTimesLH <- read.csv("~/AMPAD_Lineage/LH_reviewer_response/Monly_MF_analyses/TCX_M_PStimes.csv", stringsAsFactors = FALSE)
dlpfcLineageTimesLH <- read.csv("~/AMPAD_Lineage/LH_reviewer_response/Monly_MF_analyses/DLPFC_M_PStimes.csv", stringsAsFactors = FALSE)
tcxCovObj <- synapser::synGet("syn8466814")
tcxCov <- data.table::fread(tcxCovObj$path,data.table=F)
dlpfcCovObj <- synapser::synGet("syn11024258")
dlpfcCov <- data.table::fread(dlpfcCovObj$path,data.table=F)
tcxLineageTimesLH <- tcxLineageTimesLH[,-c(1)]
dlpfc <- dplyr::left_join(dlpfcLineageTimesLH,dlpfcCov,by='SampleID')
tcx <- dplyr::left_join(tcxLineageTimesLH,tcxCov,by='SampleID')


tcx$diagnosis<-tcx$Tissue.SourceDiagnosis
tcx$diagnosis[tcx$diagnosis=='TCX.AD'] <- 'AD'
tcx$diagnosis[tcx$diagnosis=='TCX.CONTROL'] <- 'Control'
tcx$diagnosis[tcx$diagnosis=='TCX.PSP'] <- 'PSP'
tcx$diagnosis[tcx$diagnosis=='TCX.PATH_AGE'] <- 'PA'
tcx$Study<-'MayoRNAseq'
tcxDf <- data.frame(SampleID=tcx$SampleID,
                    Pseudotime=tcx$Pseudotime,
                    Diagnosis=tcx$diagnosis,
                    Study=tcx$Study,
                    stringsAsFactors = FALSE)
tcxDf <- subset(tcxDf, tcxDf$Diagnosis=='AD'|tcxDf$Diagnosis=='Control')


dlpfc$diagnosis<-dlpfc$Diagnosis
dlpfc$diagnosis[dlpfc$diagnosis=='AD'] <- 'AD'
dlpfc$diagnosis[dlpfc$diagnosis=='CONTROL'] <- 'Control'
dlpfc$diagnosis[dlpfc$diagnosis=='OTHER'] <- 'Other'
dlpfc$Study<-'ROSMAP'
dlpfcDf <- data.frame(SampleID=dlpfc$SampleID,
                      Pseudotime=dlpfc$Pseudotime,
                      Diagnosis=dlpfc$diagnosis,
                      Study=dlpfc$Study,
                      stringsAsFactors = FALSE)
dlpfcDf <- subset(dlpfcDf, dlpfcDf$Diagnosis=='AD'|dlpfcDf$Diagnosis=='Control')

#join the two data frames and output with figure title
figureS8C <- rbind(dlpfcDf, tcxDf)

g <- ggplot2::ggplot(figureS8C,ggplot2::aes(x=Study,
                                           y=Pseudotime,
                                           color=Diagnosis))
g <- g + ggplot2::geom_boxplot() + theme(text = element_text(size = 22)) 
g <- g + ggplot2::geom_point(size=2.5, position=ggplot2::position_jitterdodge())
g <- g + ggplot2::scale_color_manual(values=viridis::viridis(3)[1:2])
g

#figure S8C, males only, pseudotime & diagnosis data 
write.csv(figureS8C, file="~/AMPAD_Lineage/LH_reviewer_response/SourceData_Synapse/figureS8C_data.csv", row.names=FALSE)
file <- File(path='~/AMPAD_Lineage/LH_reviewer_response/SourceData_Synapse/figureS8C_data.csv', parent='syn23262420')
file <- synapser::synStore(file)




#Upload pseudotimes calculated in males and females combined
tcxLineageTimesLH <- read.csv("~/AMPAD_Lineage/LH_reviewer_response/Monly_MF_analyses/TCX_MF_PStimes.csv", stringsAsFactors = FALSE)
dlpfcLineageTimesLH <- read.csv("~/AMPAD_Lineage/LH_reviewer_response/Monly_MF_analyses/DLPFC_MF_PStimes.csv", stringsAsFactors = FALSE)
tcxCovObj <- synapser::synGet("syn8466814")
tcxCov <- data.table::fread(tcxCovObj$path,data.table=F)
dlpfcCovObj <- synapser::synGet("syn11024258")
dlpfcCov <- data.table::fread(dlpfcCovObj$path,data.table=F)
tcxLineageTimesLH <- tcxLineageTimesLH[,-c(1)]
dlpfc <- dplyr::left_join(dlpfcLineageTimesLH,dlpfcCov,by='SampleID')
tcx <- dplyr::left_join(tcxLineageTimesLH,tcxCov,by='SampleID')


tcx$diagnosis<-tcx$Tissue.SourceDiagnosis
tcx$diagnosis[tcx$diagnosis=='TCX.AD'] <- 'AD'
tcx$diagnosis[tcx$diagnosis=='TCX.CONTROL'] <- 'Control'
tcx$diagnosis[tcx$diagnosis=='TCX.PSP'] <- 'PSP'
tcx$diagnosis[tcx$diagnosis=='TCX.PATH_AGE'] <- 'PA'
tcx$Study<-'MayoRNAseq'
tcxDf <- data.frame(SampleID=tcx$SampleID,
                    Pseudotime=tcx$Pseudotime,
                    Diagnosis=tcx$diagnosis,
                    Study=tcx$Study,
                    stringsAsFactors = FALSE)
tcxDf <- subset(tcxDf, tcxDf$Diagnosis=='AD'|tcxDf$Diagnosis=='Control')


dlpfc$diagnosis<-dlpfc$Diagnosis
dlpfc$diagnosis[dlpfc$diagnosis=='AD'] <- 'AD'
dlpfc$diagnosis[dlpfc$diagnosis=='CONTROL'] <- 'Control'
dlpfc$diagnosis[dlpfc$diagnosis=='OTHER'] <- 'Other'
dlpfc$Study<-'ROSMAP'
dlpfcDf <- data.frame(SampleID=dlpfc$SampleID,
                      Pseudotime=dlpfc$Pseudotime,
                      Diagnosis=dlpfc$diagnosis,
                      Study=dlpfc$Study,
                      stringsAsFactors = FALSE)
dlpfcDf <- subset(dlpfcDf, dlpfcDf$Diagnosis=='AD'|dlpfcDf$Diagnosis=='Control')

#join the two data frames and output with figure title
figureS10C <- rbind(dlpfcDf, tcxDf)

g <- ggplot2::ggplot(figureS10C,ggplot2::aes(x=Study,
                                            y=Pseudotime,
                                            color=Diagnosis))
g <- g + ggplot2::geom_boxplot() + theme(text = element_text(size = 22)) 
g <- g + ggplot2::geom_point(size=2.5, position=ggplot2::position_jitterdodge())
g <- g + ggplot2::scale_color_manual(values=viridis::viridis(3)[1:2])
g

#figure S10C, males & females combined, pseudotime & diagnosis data 
write.csv(figureS10C, file="~/AMPAD_Lineage/LH_reviewer_response/SourceData_Synapse/figureS10C_data.csv", row.names=FALSE)
file <- File(path='~/AMPAD_Lineage/LH_reviewer_response/SourceData_Synapse/figureS10C_data.csv', parent='syn23262420')
file <- synapser::synStore(file)





####figureS8D: GWAS LOAD distributions in males, figureS10D: GWAS LOAD distributions in males and females
#tcxLineageTimesLH <- read.csv("~/AMPAD_Lineage/LH_reviewer_response/Monly_MF_analyses/TCX_M_PStimes.csv", stringsAsFactors = FALSE)
tcxLineageTimesLH <- read.csv("~/AMPAD_Lineage/LH_reviewer_response/Monly_MF_analyses/TCX_MF_PStimes.csv", stringsAsFactors = FALSE)
#dlpfcLineageTimesLH <- read.csv("~/AMPAD_Lineage/LH_reviewer_response/Monly_MF_analyses/DLPFC_M_PStimes.csv", stringsAsFactors = FALSE)
dlpfcLineageTimesLH <- read.csv("~/AMPAD_Lineage/LH_reviewer_response/Monly_MF_analyses/DLPFC_MF_PStimes.csv", stringsAsFactors = FALSE)
tcxCovObj <- synapser::synGet("syn8466814")
tcxCov <- data.table::fread(tcxCovObj$path,data.table=F)
dlpfcCovObj <- synapser::synGet("syn11024258")
dlpfcCov <- data.table::fread(dlpfcCovObj$path,data.table=F)
tcxLineageTimesLH <- tcxLineageTimesLH[,-c(1)]
dlpfc <- dplyr::left_join(dlpfcLineageTimesLH,dlpfcCov,by='SampleID')
tcx <- dplyr::left_join(tcxLineageTimesLH,tcxCov,by='SampleID')



#needed to tweak the convert function to a different host (biomaRt stopped working)
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

ad_gwas <- c("CR1",
             "BIN1",
             "INPP5D",
             "HLA-DRB1",
             "TREM2",
             "MEF2C",
             "NME8",
             "CD2AP",
             "NYAP1",
             "EPHA1",
             "PTK2B",
             "CLU",
             "SPI1",
             "MS4A2",
             "PICALM",
             "SORL1",
             "FERMT2",
             "SLC24A4",
             "ABCA7",
             "APOE",
             "CASS4",
             "ECHDC3",
             "ACE",
             "NDUFAF6",
             "ECHDC3",
             "ADAMTS20",
             "SPPL2A",
             "ADAM10",
             "IQCK",
             "MIR142",
             "ACE",
             "ADAMTS1",
             "SUCLG2P4",
             "FST",
             "OARD1",
             "WWOX",
             "MAF",
             "CD55",
             "YOD1",
             "HLA-DRB1",
             "PSMB8",
             "C4A",
             "GPSM3",
             "HLA-DPA1",
             "HLA-DQA1",
             "HLA-DRA",
             "HLA-DRB5",
             "PSMB9",
             "CD2AP",
             "AGFG2",
             "PILRA",
             "EPHB4",
             "C7orf43",
             "GAL3ST4",
             "ZKSCCAN1",
             "FAM131B",
             "PSMC3",
             "ACP2",
             "C1QTNF4",
             "CELF1",
             "MTCH2",
             "NDUFS3",
             "NUP160",
             "MS4A6A",
             "MS4A7",
             "MS4A4A",
             "EED",
             "PICALM",
             "STYX",
             "RIN3",
             "HMHA1",
             "CNN2",
             "WDR18",
             "CASS4")



tcxCPMObj <- synapser::synGet('syn8466816')
#Dat <- read.delim(tcxCPMObj$path,stringsAsFactors = F)
Dat_tcx <- data.table::fread(tcxCPMObj$path,data.table=F)
sampleIds <- colnames(Dat_tcx)[-1]
geneIds <- Dat_tcx$ensembl_gene_id
Dat_tcx <- Dat_tcx[,-1]
Dat_tcx <- t(Dat_tcx)
colnames(Dat_tcx) <- geneIds
Dat_tcx <- data.frame(Dat_tcx,stringsAsFactors=F)
Dat_tcx$sampleId <- sampleIds
tcx <- dplyr::left_join(tcx,Dat_tcx,by=c('SampleID'='sampleId'))
tcx2 <- dplyr::select(tcx,dplyr::starts_with("ENSG"))

corvec <- cor(tcx2,tcx$Pseudotime,method='spearman')
corDf <- data.frame(geneid=colnames(tcx2),cor=corvec,stringsAsFactors=F)
corDf2 <- corDf
corDf2$external_gene_name <- Make.Gene.Symb(corDf2$geneid)
corDf2$cor <- NULL
corDf2 <- dplyr::left_join(corDf,corDf2,by=c('geneid'))

dlpfcCPMObj <- synapser::synGet('syn8456638')
Dat <- data.table::fread(dlpfcCPMObj$path,data.table=F)
sampleIds <- colnames(Dat)[-1]
geneIds <- Dat$ensembl_gene_id
Dat <- Dat[,-1]
Dat <- t(Dat)
colnames(Dat) <- geneIds
Dat <- data.frame(Dat,stringsAsFactors=F)
Dat$sampleId <- sampleIds
dlpfc <- dplyr::left_join(dlpfc,Dat,by=c('SampleID'='sampleId'))
dlpfc2 <- dplyr::select(dlpfc,dplyr::starts_with("ENSG"))

corvec <- cor(dlpfc2,dlpfc$Pseudotime,method='spearman')
corDfdlpfc <- data.frame(geneid=colnames(dlpfc2),cor=corvec,stringsAsFactors=F)
corDfdlpfc2 <- corDfdlpfc
corDfdlpfc2$external_gene_name <- Make.Gene.Symb(corDfdlpfc2$geneid)
corDfdlpfc2$cor <- NULL
corDfdlpfc2 <- dplyr::left_join(corDfdlpfc,corDfdlpfc2,by=c('geneid'))


corDf <- corDf2
mean(abs(corDf$cor))
mean(abs(corDf[corDf$external_gene_name %in% ad_gwas,]$cor))

corDfdlpfc <- corDfdlpfc2
mean(abs(corDfdlpfc$cor))
mean(abs(corDfdlpfc[corDfdlpfc$external_gene_name %in% ad_gwas,]$cor))

corDf$adGwas <- corDf$external_gene_name %in% ad_gwas
corDfdlpfc$adGwas <- corDfdlpfc$external_gene_name %in% ad_gwas
corDf$brainRegion <- 'TCX'
corDfdlpfc$brainRegion <- 'DLPFC'
corDfcombined <- rbind(corDf,corDfdlpfc)
colnames(corDfcombined)[4] <- 'LOADGWASGene'
g <- ggplot2::ggplot(corDfcombined,ggplot2::aes(x=brainRegion,y=cor,fill=LOADGWASGene))
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(x = 'Brain Region',y='Correlation with PMI-adj pseudotime',fill='LOAD\nGWAS\nGene')
g

#figure S8D is for male only pseudotimes:
figureS8D<-corDfcombined
write.csv(figureS8D, file="~/AMPAD_Lineage/LH_reviewer_response/SourceData_Synapse/figureS8D_data.csv", row.names=FALSE)
file <- synapser::File(path='~/AMPAD_Lineage/LH_reviewer_response/SourceData_Synapse/figureS8D_data.csv', parent='syn23262420')
file <- synapser::synStore(file)

#figure S10D is for males and females combined pseudotimes:
figureS10D<-corDfcombined
write.csv(figureS10D, file="~/AMPAD_Lineage/LH_reviewer_response/SourceData_Synapse/figureS10D_data.csv", row.names=FALSE)
file <- synapser::File(path='~/AMPAD_Lineage/LH_reviewer_response/SourceData_Synapse/figureS10D_data.csv', parent='syn23262420')
file <- synapser::synStore(file)




#### male only dlpfc neuropath data:
#read in data files to run monocle:
temp <- readRDS(file="~/AMPAD_Lineage/LH_reviewer_response/Monly_MF_analyses/DLPFC_M_countmatrix_Mono.rds")
temp2 <- readRDS(file="~/AMPAD_Lineage/LH_reviewer_response/Monly_MF_analyses/DLPFC_M_cell_metadata_Mono.rds")
gene_short_name <- readRDS(file="~/AMPAD_Lineage/LH_reviewer_response/Monly_MF_analyses/DLPFC_gene_metadata_Mono.rds")

#running monocle
rownames(temp) <- NULL
colnames(temp) <- NULL

#run this function then check to see if pseudotime is ordered correctly. rerun with order reversed if necessary.
RunMonocleTobit <- function(Dat, Labels, max_components=2, meth = 'DDRTree',C_by = NULL, 
                            gene_short_name = NULL){ 
  
  library(monocle)
  
  HSMM_expr_matrix <- Dat
  names(HSMM_expr_matrix)<-seq(1,dim(Dat)[2])
  
  if(is.null(gene_short_name)){
    gene_short_name <- c(1:dim(Dat)[1])
  }
  
  
  gene_short_name <- data.frame(gene_short_name)
  Labels <- data.frame(Labels)
  rownames(Labels) <- seq(1,dim(Dat)[2])
  
  pd <- new("AnnotatedDataFrame", data = Labels)
  fd <- new("AnnotatedDataFrame", data = gene_short_name)
  
  
  
  HSMM <- newCellDataSet(as.matrix(HSMM_expr_matrix),
                         phenoData = pd,
                         featureData = fd,
                         expressionFamily=tobit())
  
  HSMM <- reduceDimension(HSMM, max_components=max_components, reduction_method = meth)
  #HSMM <- orderCells(HSMM)
  HSMM <- orderCells(HSMM, reverse=TRUE)
  if(is.null(C_by)){
    plot_cell_trajectory(HSMM, color_by="Labels")
  }
  else{
    plot_cell_trajectory(HSMM, color_by=C_by)
  }
  
  
  return(HSMM)
  
}
MonRun <- RunMonocleTobit(temp, temp2, C_by = 'Pseudotime',gene_short_name = gene_short_name)
plot_cell_trajectory(MonRun, color_by = "Diagnosis")

#Make new dataframe with SampleID, Pseudotime, braaksc
x <- list()

x$SampleID <- MonRun$SampleID
x$Pseudotime <- MonRun$Pseudotime
x$Braak <- MonRun$braaksc
x$CERAD <- MonRun$ceradsc
x$cogdx <- MonRun$cogdx

x <- as.data.frame(x)
x$SampleID <- as.character(x$SampleID)

write.csv(x,'~/AMPAD_Lineage/LH_reviewer_response/SourceData_Synapse/figureS9DEF.csv', row.names=FALSE)
#write.csv(x,'~/AMPAD_Lineage/LH_reviewer_response/SourceData_Synapse/figureS9DEF.csv', row.names=FALSE)
file <- synapser::File(path='~/AMPAD_Lineage/LH_reviewer_response/SourceData_Synapse/figureS9DEF.csv', parent='syn23262420')
file <- synapser::synStore(file)




#### males and females combined dlpfc neuropath data:
#read in data files to simply run monocle:
gene_short_name <- readRDS(file="~/AMPAD_Lineage/LH_reviewer_response/Monly_MF_analyses/DLPFC_gene_metadata_Mono.rds")
temp <- readRDS(file="~/AMPAD_Lineage/LH_reviewer_response/Monly_MF_analyses/DLPFC_MF_countmatrix_Mono.rds")
temp2 <- readRDS(file="~/AMPAD_Lineage/LH_reviewer_response/Monly_MF_analyses/DLPFC_MF_cell_metadata_Mono.rds")

#running monocle
rownames(temp) <- NULL
colnames(temp) <- NULL

#run this function then check to see if pseudotime is ordered correctly. rerun with order reversed if necessary.
RunMonocleTobit <- function(Dat, Labels, max_components=2, meth = 'DDRTree',C_by = NULL, 
                            gene_short_name = NULL){ 
  
  library(monocle)
  
  HSMM_expr_matrix <- Dat
  names(HSMM_expr_matrix)<-seq(1,dim(Dat)[2])
  
  if(is.null(gene_short_name)){
    gene_short_name <- c(1:dim(Dat)[1])
  }
  
  
  gene_short_name <- data.frame(gene_short_name)
  Labels <- data.frame(Labels)
  rownames(Labels) <- seq(1,dim(Dat)[2])
  
  pd <- new("AnnotatedDataFrame", data = Labels)
  fd <- new("AnnotatedDataFrame", data = gene_short_name)
  
  
  
  HSMM <- newCellDataSet(as.matrix(HSMM_expr_matrix),
                         phenoData = pd,
                         featureData = fd,
                         expressionFamily=tobit())
  
  HSMM <- reduceDimension(HSMM, max_components=max_components, reduction_method = meth)
  #HSMM <- orderCells(HSMM)
  HSMM <- orderCells(HSMM, reverse=TRUE)
  if(is.null(C_by)){
    plot_cell_trajectory(HSMM, color_by="Labels")
  }
  else{
    plot_cell_trajectory(HSMM, color_by=C_by)
  }
  
  
  return(HSMM)
  
}
MonRun <- RunMonocleTobit(temp, temp2, C_by = 'Pseudotime',gene_short_name = gene_short_name)
plot_cell_trajectory(MonRun, color_by = "Diagnosis")

#Make new dataframe with SampleID, Pseudotime, braaksc, cerad, cogdx
x <- list()

x$SampleID <- MonRun$SampleID
x$Pseudotime <- MonRun$Pseudotime
x$Braak <- MonRun$braaksc
x$CERAD <- MonRun$ceradsc
x$cogdx <- MonRun$cogdx

x <- as.data.frame(x)
x$SampleID <- as.character(x$SampleID)

write.csv(x,'~/AMPAD_Lineage/LH_reviewer_response/SourceData_Synapse/figureS11DEF.csv', row.names=FALSE)
file <- synapser::File(path='~/AMPAD_Lineage/LH_reviewer_response/SourceData_Synapse/figureS11DEF.csv', parent='syn23262420')
file <- synapser::synStore(file)




#rosmap Braak-adjusted pseudotimes (figureS13: reload the data & rerun monocle with braak adjustment in the function:
temp <- readRDS(file="~/AMPAD_Lineage/LH_reviewer_response/DLPFC_countmatrix_Mono.rds")
temp2 <- readRDS(file="~/AMPAD_Lineage/LH_reviewer_response/DLPFC_cell_metadata_Mono.rds")
gene_short_name <- readRDS(file="~/AMPAD_Lineage/LH_reviewer_response/DLPFC_gene_metadata_Mono.rds")

#running monocle
rownames(temp) <- NULL
colnames(temp) <- NULL


#run this function then check to see if pseudotime is ordered correctly. rerun with order reversed if necessary.
RunMonocleTobit <- function(Dat, Labels, max_components=2, meth = 'DDRTree',C_by = NULL, 
                            gene_short_name = NULL){ 
  
  library(monocle)
  
  HSMM_expr_matrix <- Dat
  names(HSMM_expr_matrix)<-seq(1,dim(Dat)[2])
  
  if(is.null(gene_short_name)){
    gene_short_name <- c(1:dim(Dat)[1])
  }
  
  
  gene_short_name <- data.frame(gene_short_name)
  Labels <- data.frame(Labels)
  rownames(Labels) <- seq(1,dim(Dat)[2])
  
  pd <- new("AnnotatedDataFrame", data = Labels)
  fd <- new("AnnotatedDataFrame", data = gene_short_name)
  
  
  
  HSMM <- newCellDataSet(as.matrix(HSMM_expr_matrix),
                         phenoData = pd,
                         featureData = fd,
                         expressionFamily=tobit())
  
  HSMM <- reduceDimension(HSMM, max_components=max_components, reduction_method = meth, residualModelFormulaStr = "~braaksc")
  #HSMM <- orderCells(HSMM)
  HSMM <- orderCells(HSMM, reverse=TRUE)
  if(is.null(C_by)){
    plot_cell_trajectory(HSMM, color_by="Labels")
  }
  else{
    plot_cell_trajectory(HSMM, color_by=C_by)
  }
  
  
  return(HSMM)
  
}
MonRun <- RunMonocleTobit(temp, temp2, C_by = 'Pseudotime',gene_short_name = gene_short_name)
plot_cell_trajectory(MonRun, color_by = "Diagnosis")

g<- plot_cell_trajectory(MonRun,color_by = "Diagnosis",show_branch_points=F,use_color_gradient = F,cell_size = 1)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="Diagnosis")
g

MonRun$cogdx<-as.character(MonRun$cogdx)
g <- ggplot2::ggplot(MonRun@phenoData@data, aes(x=cogdx, y=scale(Pseudotime,center=F),fill=cogdx)) 
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::stat_summary(fun=mean, geom="point", shape=23, size=2)
g <- g + ggplot2::theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
                        legend.text=element_text(size=15)) 
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(fill="Cognitive\nDiagnosis",y="Pseudotime",x="Cognitive Diagnosis")
g

#Make new dataframe with SampleID, Pseudotime, cogdx for figure S13D
x <- list()
x$SampleID <- MonRun$SampleID
x$Pseudotime_braakadj <- MonRun$Pseudotime
x$cogdx <- MonRun$cogdx

x <- as.data.frame(x)
x$SampleID <- as.character(x$SampleID)

write.csv(x,'~/AMPAD_Lineage/LH_reviewer_response/SourceData_Synapse/figureS13D.csv', row.names=FALSE)
file <- synapser::File(path='~/AMPAD_Lineage/LH_reviewer_response/SourceData_Synapse/figureS13D.csv', parent='syn23262420')
file <- synapser::synStore(file)

#figure S13B: braak adjusted pseudotime-AD status association
dlpfcBraakadjPstimes <- read.csv("~/AMPAD_Lineage/LH_reviewer_response/SourceData_Synapse/figureS13D.csv", stringsAsFactors = FALSE)
dlpfcBraakadjPstimes$cogdx<-NULL
dlpfcCovObj <- synapser::synGet("syn11024258")
dlpfcCov <- data.table::fread(dlpfcCovObj$path,data.table=F)
dlpfc <- dplyr::left_join(dlpfcBraakadjPstimes,dlpfcCov,by='SampleID')

dlpfc$diagnosis<-dlpfc$Diagnosis
dlpfc$diagnosis[dlpfc$diagnosis=='AD'] <- 'AD'
dlpfc$diagnosis[dlpfc$diagnosis=='CONTROL'] <- 'Control'
dlpfc$diagnosis[dlpfc$diagnosis=='OTHER'] <- 'Other'
dlpfc$Study<-'ROSMAP'
dlpfcDf <- data.frame(SampleID=dlpfc$SampleID,
                      Pseudotime_braakadj=dlpfc$Pseudotime_braakadj,
                      Diagnosis=dlpfc$diagnosis,
                      Study=dlpfc$Study,
                      stringsAsFactors = FALSE)
dlpfcDf <- subset(dlpfcDf, dlpfcDf$Diagnosis=='AD'|dlpfcDf$Diagnosis=='Control')

g <- ggplot2::ggplot(dlpfcDf,ggplot2::aes(x=Diagnosis,
                                           y=Pseudotime_braakadj,
                                           color=Diagnosis))
g <- g + ggplot2::geom_boxplot() + theme(text = element_text(size = 22)) 
g <- g + ggplot2::geom_point(size=2.5, position=ggplot2::position_jitterdodge())
g <- g + ggplot2::scale_color_manual(values=viridis::viridis(3)[1:2])
g

write.csv(dlpfcDf, file="~/AMPAD_Lineage/LH_reviewer_response/SourceData_Synapse/figureS13B_data.csv", row.names=FALSE)
file <- File(path='~/AMPAD_Lineage/LH_reviewer_response/SourceData_Synapse/figureS13B_data.csv', parent='syn23262420')
file <- synapser::synStore(file)

######### figure S13E, gwas correlations with braak-adjusted pseudotime
# dlpfcLineageTimesLH <- read.csv("~/AMPAD_Lineage/LH_reviewer_response/DLPFC_F_PStimes_adjusted.csv", stringsAsFactors = FALSE)
# dlpfcCovObj <- synapser::synGet("syn11024258")
# dlpfcCov <- data.table::fread(dlpfcCovObj$path,data.table=F)
# dlpfc <- dplyr::left_join(dlpfcLineageTimesLH,dlpfcCov,by='SampleID')

#needed to tweak the convert function to a different host (biomaRt stopped working)
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

ad_gwas <- c("CR1",
             "BIN1",
             "INPP5D",
             "HLA-DRB1",
             "TREM2",
             "MEF2C",
             "NME8",
             "CD2AP",
             "NYAP1",
             "EPHA1",
             "PTK2B",
             "CLU",
             "SPI1",
             "MS4A2",
             "PICALM",
             "SORL1",
             "FERMT2",
             "SLC24A4",
             "ABCA7",
             "APOE",
             "CASS4",
             "ECHDC3",
             "ACE",
             "NDUFAF6",
             "ECHDC3",
             "ADAMTS20",
             "SPPL2A",
             "ADAM10",
             "IQCK",
             "MIR142",
             "ACE",
             "ADAMTS1",
             "SUCLG2P4",
             "FST",
             "OARD1",
             "WWOX",
             "MAF",
             "CD55",
             "YOD1",
             "HLA-DRB1",
             "PSMB8",
             "C4A",
             "GPSM3",
             "HLA-DPA1",
             "HLA-DQA1",
             "HLA-DRA",
             "HLA-DRB5",
             "PSMB9",
             "CD2AP",
             "AGFG2",
             "PILRA",
             "EPHB4",
             "C7orf43",
             "GAL3ST4",
             "ZKSCCAN1",
             "FAM131B",
             "PSMC3",
             "ACP2",
             "C1QTNF4",
             "CELF1",
             "MTCH2",
             "NDUFS3",
             "NUP160",
             "MS4A6A",
             "MS4A7",
             "MS4A4A",
             "EED",
             "PICALM",
             "STYX",
             "RIN3",
             "HMHA1",
             "CNN2",
             "WDR18",
             "CASS4")


dlpfcCPMObj <- synapser::synGet('syn8456638')
Dat <- data.table::fread(dlpfcCPMObj$path,data.table=F)
sampleIds <- colnames(Dat)[-1]
geneIds <- Dat$ensembl_gene_id
Dat <- Dat[,-1]
Dat <- t(Dat)
colnames(Dat) <- geneIds
Dat <- data.frame(Dat,stringsAsFactors=F)
Dat$sampleId <- sampleIds
dlpfc <- dplyr::left_join(dlpfc,Dat,by=c('SampleID'='sampleId'))
dlpfc2 <- dplyr::select(dlpfc,dplyr::starts_with("ENSG"))

corvec <- cor(dlpfc2,dlpfc$Pseudotime_braakadj,method='spearman')
corDfdlpfc <- data.frame(geneid=colnames(dlpfc2),cor=corvec,stringsAsFactors=F)
corDfdlpfc2 <- corDfdlpfc
corDfdlpfc2$external_gene_name <- Make.Gene.Symb(corDfdlpfc2$geneid)
corDfdlpfc2$cor <- NULL
corDfdlpfc2 <- dplyr::left_join(corDfdlpfc,corDfdlpfc2,by=c('geneid'))

corDfdlpfc <- corDfdlpfc2
mean(abs(corDfdlpfc$cor))
mean(abs(corDfdlpfc[corDfdlpfc$external_gene_name %in% ad_gwas,]$cor))

corDfdlpfc$adGwas <- corDfdlpfc$external_gene_name %in% ad_gwas
corDfdlpfc$brainRegion <- 'DLPFC'
colnames(corDfdlpfc)[4] <- 'LOADGWASGene'
g <- ggplot2::ggplot(corDfdlpfc,ggplot2::aes(x=brainRegion,y=cor,fill=LOADGWASGene))
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(x = 'Brain Region',y='Correlation with braak-adj pseudotime',fill='LOAD\nGWAS\nGene')
g

#figure S13E is for braak-adjusted pseudotimes:
write.csv(corDfdlpfc, file="~/AMPAD_Lineage/LH_reviewer_response/SourceData_Synapse/figureS13E_data.csv", row.names=FALSE)
file <- synapser::File(path='~/AMPAD_Lineage/LH_reviewer_response/SourceData_Synapse/figureS13E_data.csv', parent='syn23262420')
file <- synapser::synStore(file)




## figure S14B: tcx neuropath & pseudotime (Braak & Thal amyloid)
temp <- readRDS(file="~/AMPAD_Lineage/LH_reviewer_response/TCX_countmatrix_Mono.rds")
temp2 <- readRDS(file="~/AMPAD_Lineage/LH_reviewer_response/TCX_cell_metadata_Mono.rds")
tcxneuroobj <- synapser::synGet("syn3817650")
tcxneuro <- data.table::fread(tcxneuroobj$path,data.table=F)
tcxneuro<-subset(tcxneuro, select=c(SampleID, Braak, Thal))
tcxneuro$Braak2<-ifelse(tcxneuro$Braak==0,0,
                       ifelse(tcxneuro$Braak==0.5,1,
                              ifelse(tcxneuro$Braak==1,1,
                                     ifelse(tcxneuro$Braak==1.5,2,
                                            ifelse(tcxneuro$Braak==2,2,
                                                   ifelse(tcxneuro$Braak==2.5,3,
                                                          ifelse(tcxneuro$Braak==3,3,
                                                                 ifelse(tcxneuro$Braak==4.5,5,
                                                                        ifelse(tcxneuro$Braak==5,5,
                                                                               ifelse(tcxneuro$Braak==5.5,6,
                                                                                      ifelse(tcx$neuro$Braak==6,6,'NA')))))))))))
temp2<-dplyr::left_join(temp2, tcxneuro, by='SampleID')
gene_short_name <- readRDS(file="~/AMPAD_Lineage/LH_reviewer_response/TCX_gene_metadata_Mono.rds")

RunMonocleTobit <- function(Dat, Labels, max_components=2, meth = 'DDRTree',C_by = NULL, 
                            gene_short_name = NULL){ 
  
  library(monocle)
  
  HSMM_expr_matrix <- Dat
  names(HSMM_expr_matrix)<-seq(1,dim(Dat)[2])
  
  if(is.null(gene_short_name)){
    gene_short_name <- c(1:dim(Dat)[1])
  }
  
  
  gene_short_name <- data.frame(gene_short_name)
  Labels <- data.frame(Labels)
  rownames(Labels) <- seq(1,dim(Dat)[2])
  
  pd <- new("AnnotatedDataFrame", data = Labels)
  fd <- new("AnnotatedDataFrame", data = gene_short_name)
  
  
  
  HSMM <- newCellDataSet(as.matrix(HSMM_expr_matrix),
                         phenoData = pd,
                         featureData = fd,
                         expressionFamily=tobit())
  
  HSMM <- reduceDimension(HSMM, max_components=max_components, reduction_method = meth)
  HSMM <- orderCells(HSMM)
  #HSMM <- orderCells(HSMM, reverse=TRUE)
  if(is.null(C_by)){
    plot_cell_trajectory(HSMM, color_by="Labels")
  }
  else{
    plot_cell_trajectory(HSMM, color_by=C_by)
  }
  
  
  return(HSMM)
  
}
MonRun <- RunMonocleTobit(temp, temp2, C_by = 'Pseudotime',gene_short_name = gene_short_name)
plot_cell_trajectory(MonRun, color_by = "Tissue.Diagnosis")
MonRun$Braak<-as.character(MonRun$Braak2)
g<- plot_cell_trajectory(MonRun,color_by = "Braak",show_branch_points=F,use_color_gradient = F,cell_size = 1.5)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="Diagnosis")
g
#dev.off()
#make sure pseudotime is oriented in the right direction
g <- ggplot2::ggplot(MonRun@phenoData@data,ggplot2::aes(x=Tissue.Diagnosis,
                                                        y=Pseudotime), fill=Tissue.Diagnosis)
g <- g + ggplot2::geom_boxplot() + theme(text = element_text(size = 22)) 
g


g <- ggplot2::ggplot(MonRun@phenoData@data, aes(x=Braak, y=scale(Pseudotime,center=F),fill=Braak)) 
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::stat_summary(fun=mean, geom="point", shape=23, size=2)
g <- g + ggplot2::theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
                        legend.text=element_text(size=15)) 
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(fill="Braak Score",y="Pseudotime",x="Braak Score")
g

MonRun$Thal<-as.character(MonRun$Thal)
g <- ggplot2::ggplot(MonRun@phenoData@data, aes(x=Thal, y=scale(Pseudotime,center=F),fill=Thal)) 
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::stat_summary(fun=mean, geom="point", shape=23, size=2)
g <- g + ggplot2::theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
                        legend.text=element_text(size=15)) 
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(fill="Thal amyloid",y="Pseudotime",x="Thal amyloid")
g

x <- list()
x$SampleID <- MonRun$SampleID
x$Pseudotime <- MonRun$Pseudotime
x$Braak <- MonRun$Braak
x$Thal <- MonRun$Thal

x <- as.data.frame(x)
x$SampleID <- as.character(x$SampleID)

#figure S14B is for tcx braak & thal pseudotime associations:
write.csv(x, file="~/AMPAD_Lineage/LH_reviewer_response/SourceData_Synapse/figureS14B_data.csv", row.names=FALSE)
file <- synapser::File(path='~/AMPAD_Lineage/LH_reviewer_response/SourceData_Synapse/figureS14B_data.csv', parent='syn23262420')
file <- synapser::synStore(file)




##### FigureS17C: pseudotime-ad associations using UMAP-based pseudotimes from Monocle3
#Upload pseudotimes 
tcxLineageTimesLH <- read.csv("~/AMPAD_Lineage/LH_reviewer_response/TCX_pstimes_Mono3.csv", stringsAsFactors = FALSE)
dlpfcLineageTimesLH <- read.csv("~/AMPAD_Lineage/LH_reviewer_response/DLPFC_pstimes_Mono3.csv", stringsAsFactors = FALSE)
tcxCovObj <- synapser::synGet("syn8466814")
tcxCov <- data.table::fread(tcxCovObj$path,data.table=F)
dlpfcCovObj <- synapser::synGet("syn11024258")
dlpfcCov <- data.table::fread(dlpfcCovObj$path,data.table=F)
tcxLineageTimesLH <- tcxLineageTimesLH[,-c(1)]
dlpfc <- dplyr::left_join(dlpfcLineageTimesLH,dlpfcCov,by='SampleID')
tcx <- dplyr::left_join(tcxLineageTimesLH,tcxCov,by='SampleID')


tcx$diagnosis<-tcx$Tissue.SourceDiagnosis
tcx$diagnosis[tcx$diagnosis=='TCX.AD'] <- 'AD'
tcx$diagnosis[tcx$diagnosis=='TCX.CONTROL'] <- 'Control'
tcx$diagnosis[tcx$diagnosis=='TCX.PSP'] <- 'PSP'
tcx$diagnosis[tcx$diagnosis=='TCX.PATH_AGE'] <- 'PA'
tcx$Study<-'MayoRNAseq'
tcxDf <- data.frame(SampleID=tcx$SampleID,
                    Pseudotime=tcx$pstime_Mono3,
                    Diagnosis=tcx$diagnosis,
                    Study=tcx$Study,
                    stringsAsFactors = FALSE)
tcxDf <- subset(tcxDf, tcxDf$Diagnosis=='AD'|tcxDf$Diagnosis=='Control')


dlpfc$diagnosis<-dlpfc$Diagnosis
dlpfc$diagnosis[dlpfc$diagnosis=='AD'] <- 'AD'
dlpfc$diagnosis[dlpfc$diagnosis=='CONTROL'] <- 'Control'
dlpfc$diagnosis[dlpfc$diagnosis=='OTHER'] <- 'Other'
dlpfc$Study<-'ROSMAP'
dlpfcDf <- data.frame(SampleID=dlpfc$SampleID,
                      Pseudotime=dlpfc$pstime_Mono3,
                      Diagnosis=dlpfc$diagnosis,
                      Study=dlpfc$Study,
                      stringsAsFactors = FALSE)
dlpfcDf <- subset(dlpfcDf, dlpfcDf$Diagnosis=='AD'|dlpfcDf$Diagnosis=='Control')

#join the two data frames and output with figure title
figureS17C <- rbind(dlpfcDf, tcxDf)

g <- ggplot2::ggplot(figureS17C,ggplot2::aes(x=Study,
                                            y=Pseudotime,
                                            color=Diagnosis))
g <- g + ggplot2::geom_boxplot() + theme(text = element_text(size = 22)) 
g <- g + ggplot2::geom_point(size=2.5, position=ggplot2::position_jitterdodge())
g <- g + ggplot2::scale_color_manual(values=viridis::viridis(3)[1:2])
g

#figure S17C, Monocle3-UMAP-derived pseudotime & diagnosis data 
write.csv(figureS17C, file="~/AMPAD_Lineage/LH_reviewer_response/SourceData_Synapse/figureS17C_data.csv", row.names=FALSE)
file <- synapser::File(path='~/AMPAD_Lineage/LH_reviewer_response/SourceData_Synapse/figureS17C_data.csv', parent='syn23262420')
file <- synapser::synStore(file)




####figureS17D: GWAS LOAD distributions for Monocle3 UMAP pseudotimes
# tcxLineageTimesLH <- read.csv("~/AMPAD_Lineage/LH_reviewer_response/TCX_pstimes_Mono3.csv", stringsAsFactors = FALSE)
# dlpfcLineageTimesLH <- read.csv("~/AMPAD_Lineage/LH_reviewer_response/DLPFC_pstimes_Mono3.csv", stringsAsFactors = FALSE)
# tcxCovObj <- synapser::synGet("syn8466814")
# tcxCov <- data.table::fread(tcxCovObj$path,data.table=F)
# dlpfcCovObj <- synapser::synGet("syn11024258")
# dlpfcCov <- data.table::fread(dlpfcCovObj$path,data.table=F)
# tcxLineageTimesLH <- tcxLineageTimesLH[,-c(1)]
# dlpfc <- dplyr::left_join(dlpfcLineageTimesLH,dlpfcCov,by='SampleID')
# tcx <- dplyr::left_join(tcxLineageTimesLH,tcxCov,by='SampleID')
# 


#needed to tweak the convert function to a different host (biomaRt stopped working)
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

ad_gwas <- c("CR1",
             "BIN1",
             "INPP5D",
             "HLA-DRB1",
             "TREM2",
             "MEF2C",
             "NME8",
             "CD2AP",
             "NYAP1",
             "EPHA1",
             "PTK2B",
             "CLU",
             "SPI1",
             "MS4A2",
             "PICALM",
             "SORL1",
             "FERMT2",
             "SLC24A4",
             "ABCA7",
             "APOE",
             "CASS4",
             "ECHDC3",
             "ACE",
             "NDUFAF6",
             "ECHDC3",
             "ADAMTS20",
             "SPPL2A",
             "ADAM10",
             "IQCK",
             "MIR142",
             "ACE",
             "ADAMTS1",
             "SUCLG2P4",
             "FST",
             "OARD1",
             "WWOX",
             "MAF",
             "CD55",
             "YOD1",
             "HLA-DRB1",
             "PSMB8",
             "C4A",
             "GPSM3",
             "HLA-DPA1",
             "HLA-DQA1",
             "HLA-DRA",
             "HLA-DRB5",
             "PSMB9",
             "CD2AP",
             "AGFG2",
             "PILRA",
             "EPHB4",
             "C7orf43",
             "GAL3ST4",
             "ZKSCCAN1",
             "FAM131B",
             "PSMC3",
             "ACP2",
             "C1QTNF4",
             "CELF1",
             "MTCH2",
             "NDUFS3",
             "NUP160",
             "MS4A6A",
             "MS4A7",
             "MS4A4A",
             "EED",
             "PICALM",
             "STYX",
             "RIN3",
             "HMHA1",
             "CNN2",
             "WDR18",
             "CASS4")



tcxCPMObj <- synapser::synGet('syn8466816')
#Dat <- read.delim(tcxCPMObj$path,stringsAsFactors = F)
Dat_tcx <- data.table::fread(tcxCPMObj$path,data.table=F)
sampleIds <- colnames(Dat_tcx)[-1]
geneIds <- Dat_tcx$ensembl_gene_id
Dat_tcx <- Dat_tcx[,-1]
Dat_tcx <- t(Dat_tcx)
colnames(Dat_tcx) <- geneIds
Dat_tcx <- data.frame(Dat_tcx,stringsAsFactors=F)
Dat_tcx$sampleId <- sampleIds
tcx <- dplyr::left_join(tcx,Dat_tcx,by=c('SampleID'='sampleId'))
tcx2 <- dplyr::select(tcx,dplyr::starts_with("ENSG"))

corvec <- cor(tcx2,tcx$pstime_Mono3,method='spearman')
corDf <- data.frame(geneid=colnames(tcx2),cor=corvec,stringsAsFactors=F)
corDf2 <- corDf
corDf2$external_gene_name <- Make.Gene.Symb(corDf2$geneid)
corDf2$cor <- NULL
corDf2 <- dplyr::left_join(corDf,corDf2,by=c('geneid'))

dlpfcCPMObj <- synapser::synGet('syn8456638')
Dat <- data.table::fread(dlpfcCPMObj$path,data.table=F)
sampleIds <- colnames(Dat)[-1]
geneIds <- Dat$ensembl_gene_id
Dat <- Dat[,-1]
Dat <- t(Dat)
colnames(Dat) <- geneIds
Dat <- data.frame(Dat,stringsAsFactors=F)
Dat$sampleId <- sampleIds
dlpfc <- dplyr::left_join(dlpfc,Dat,by=c('SampleID'='sampleId'))
dlpfc2 <- dplyr::select(dlpfc,dplyr::starts_with("ENSG"))

corvec <- cor(dlpfc2,dlpfc$pstime_Mono3,method='spearman')
corDfdlpfc <- data.frame(geneid=colnames(dlpfc2),cor=corvec,stringsAsFactors=F)
corDfdlpfc2 <- corDfdlpfc
corDfdlpfc2$external_gene_name <- Make.Gene.Symb(corDfdlpfc2$geneid)
corDfdlpfc2$cor <- NULL
corDfdlpfc2 <- dplyr::left_join(corDfdlpfc,corDfdlpfc2,by=c('geneid'))


corDf <- corDf2
mean(abs(corDf$cor))
mean(abs(corDf[corDf$external_gene_name %in% ad_gwas,]$cor))

corDfdlpfc <- corDfdlpfc2
mean(abs(corDfdlpfc$cor))
mean(abs(corDfdlpfc[corDfdlpfc$external_gene_name %in% ad_gwas,]$cor))

corDf$adGwas <- corDf$external_gene_name %in% ad_gwas
corDfdlpfc$adGwas <- corDfdlpfc$external_gene_name %in% ad_gwas
corDf$brainRegion <- 'TCX'
corDfdlpfc$brainRegion <- 'DLPFC'
corDfcombined <- rbind(corDf,corDfdlpfc)
colnames(corDfcombined)[4] <- 'LOADGWASGene'
g <- ggplot2::ggplot(corDfcombined,ggplot2::aes(x=brainRegion,y=cor,fill=LOADGWASGene))
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(x = 'Brain Region',y='Correlation with UMAP pseudotime',fill='LOAD\nGWAS\nGene')
g

#figure S17D (UMAP-based pseudotime-gwas correlation)
figureS17D<-corDfcombined
write.csv(figureS17D, file="~/AMPAD_Lineage/LH_reviewer_response/SourceData_Synapse/figureS17D_data.csv", row.names=FALSE)
file <- synapser::File(path='~/AMPAD_Lineage/LH_reviewer_response/SourceData_Synapse/figureS17D_data.csv', parent='syn23262420')
file <- synapser::synStore(file)




########## FigureS17DEF: neuropath-pseudotime associations for umap-derived pseudotime in rosmap
dlpfcLineageTimesLH <- read.csv("~/AMPAD_Lineage/LH_reviewer_response/DLPFC_pstimes_Mono3.csv", stringsAsFactors = FALSE)
dlpfcCovObj <- synapser::synGet("syn11024258")
dlpfcCov <- data.table::fread(dlpfcCovObj$path,data.table=F)
dlpfc <- dplyr::left_join(dlpfcLineageTimesLH,dlpfcCov,by='SampleID')


rosmapObj <- synapser::synGet('syn3191087')
rosmap <- data.table::fread(rosmapObj$path,data.table=F)

#add pseudotime & cogdx to braak score & cerad score
rosmapIdObj <- synapser::synGet('syn3382527')
rosmapId <- data.table::fread(rosmapIdObj$path,data.table=F)
rosmapId <- dplyr::select(rosmapId,projid,rnaseq_id)
rosmapRNAid<-dplyr::left_join(rosmapId,rosmap)
#remove duplicate rows
rosmapRNAid <- unique(rosmapRNAid)
rosmapRNAid2 <- subset(rosmapRNAid, select=c(rnaseq_id,braaksc,ceradsc))
names(rosmapRNAid2)[names(rosmapRNAid2) == "rnaseq_id"] <- "SampleID"

dlpfc3<-dplyr::left_join(dlpfc,rosmapRNAid2, by="SampleID")

#Make new dataframe with SampleID, Pseudotime, braaksc
figureS18DEF <- list()

figureS18DEF$SampleID <- dlpfc3$SampleID
figureS18DEF$Pseudotime <- dlpfc3$pstime_Mono3
figureS18DEF$Braak <- dlpfc3$braaksc
figureS18DEF$CERAD <- dlpfc3$ceradsc
figureS18DEF$cogdx <- dlpfc3$cogdx

figureS18DEF <- as.data.frame(figureS18DEF)

write.csv(figureS18DEF,'~/AMPAD_Lineage/LH_reviewer_response/SourceData_Synapse/figureS18DEF.csv', row.names=FALSE)
file <- synapser::File(path='~/AMPAD_Lineage/LH_reviewer_response/SourceData_Synapse/figureS18DEF.csv', parent='syn23262420')
file <- synapser::synStore(file)
