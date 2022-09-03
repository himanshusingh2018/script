#load all datasets
load('All_62_ARACNE.rda')#load all ARACNe Network
load('ALL_62_ARACNE_READY_EXPMAT.rda')#VST-normalized genes expression value
load('All_62_pancancer_tpm.rda')#TPM-normalized gene expression

#Extract Median Expression of gene from pancan_tpm data
#cbind(pancan_tpm_count[['AdiposeSub']][1:5,1:3], median = matrixStats::rowMedians(pancan_tpm_count[['AdiposeSub']][1:5,1:3]))[,'median']

#Extract median MI (target gene from tissue wise network) and median expression from vst_expression
exp_mi_info <- function(expDATA, tissueNET){
  #Extract Median Expression of gene in specific tissue
  expMed <- cbind(expDATA,
                  median = matrixStats::rowMedians(expDATA))[,'median']
  
  #matrix to find median of target gene present in the network
  mat <- matrix(nrow=0, ncol=2)
  for(tg in names(expMed)){   #iterate each genes with median expression
    mi <- NULL
    for(inet in tissueNET[[1]][,2]){     #iterate over target genes of each TF of the network [,2] is index of network
      if(tg %in% row.names( tissueNET[[2]][[inet]]) ){   #if target gene is present in network
        mi <- append(mi,
                     tissueNET[[2]][[inet]][as.character(tg),'MI']) #extract mutual information
      }
    }
    #matrix with target gene in network and median mutual information
    #It also add gene_id not present in the target gene. In case we need to check. MI value is also gene id
    mat <- rbind(mat, c(tg, median(mi)) )
  }
  row.names(mat) <- mat[,1]#gene_id as row names
  mat <- mat[,-1] #remove gene_id column
  
  #Merge expMed with mutual information matrix 'mat'  
  final_mat <- merge(mat, expMed, by='row.names')#merge MI and median exp of tissue by row names(gene)
  row.names(final_mat) <- final_mat[,'Row.names']#row.names/gene_id as row names
  final_mat <- final_mat[,-1] #remove row.names column
  colnames(final_mat) <- c('MI','medianExp')
  final_mat <- final_mat[final_mat < 1, ]#remove all target genes not found in network
  return(final_mat)
}

#call exp_mi_info for all tissues in varNames
for(tissue in varNames){  
  if( tissue %in% varNames[grep('tcga', varNames)] ){      #iterate over tcga tissues only
    print(paste0('vstexp_mi_',tissue))
    
    assign(paste0('vstexp_mi_ ',tissue),     #var name vstexp_mi + tissue name such as AdiposeSub etc
           exp_mi_info(
             get(paste0('expmat_', tissue) ),
             get(tissue)  )   )    #assign median of mutual information and expression to target gene
  }
  else{          #iterate over gtex tissues only
    print(paste0('vstexp_mi_gtex_',tissue))
    assign(paste0('vstexp_mi_gtex_',tissue),
           exp_mi_info(
             get(paste0( 'expmat_gtex_', tissue)  ),
             get(tissue)  )    )
  }
} 


#tissue wise correlation matrix: correlation between medianMI and medainExp of target gene
exp_mi_tissues <- ls()[grep('vst_exp_mi', ls())] #extract all tissue matrix (mi, exp)
miExpCorMatrix<-matrix(ncol = 3, nrow=0)
for( tissue in  exp_mi_tissues){  #iterate over each vstexp_mi object
  miExpCorMatrix<-rbind(miExpCorMatrix,c(tissue,
                            cor.test(as.numeric(get(tissue)$MI),
                                     as.numeric(get(tissue)$medianExp))$estimate,
                            cor.test(as.numeric(get(tissue)$MI),
                                     as.numeric(get(tissue)$medianExp))$p.value )  )
}
colnames(miExpCorMatrix)<-c('tissue','estimate','pval')
write.table(miExpCorMatrix,file='miExpCorMatrix.txt', sep = "\t", quote = FALSE, row.names = FALSE)
save.image(file='myEnvironment.RData')

################################################
#GENERATE MATRIX TO CALCULATE AVERAGE MEAN
################################################

#create empty matrix col: mean, count
mat <- matrix(ncol=2, nrow=max( as.numeric( names(expMed) ) ), 
              dimnames = list( as.character(rep(1: max( as.numeric( names(expMed) ) ) ) ), c('mean','count') ) )
mat[,] <- 0#add value 0 to mean and count
Sys.time()
#row.names(AdiposeSub[[2]][['5046']])
for(inet in AdiposeSub[[1]][,2]){
  print(inet)
  rowIndex <- row.names(AdiposeSub[[2]][[inet]])
  
  mat[rowIndex, 1] <- mat[rowIndex, 1] + AdiposeSub[[2]][[inet]][rowIndex,2]
  mat[rowIndex, 2] <- mat[rowIndex,'count'] + 1
  
  #break
}
Sys.time()


###################### EXTRA CODE
load('All_62_ARACNE.rda')#load all ARACNe Network
load('ALL_62_ARACNE_READY_EXPMAT.rda')#VST-normalized genes expression value
load('All_62_pancancer_tpm.rda')#pancancer tpm

#Return Matrix of meanMI and meanExp
meanExpMIcorrelation <- function(expDATA,aracneNET){
  mat <- cbind(expDATA,
               expMean = matrixStats::rowMeans2(expDATA))
  #mat = expDATA
  max_id = max(as.numeric(rownames(mat)))
  mean_exp = rep(-1, max_id)
  
  for(i in 1:row(mat))
    mean_exp[as.integer(rownames(mat)[i])] = mean(mat[i,])
  
  net <- aracneNET
  mean_mi = rep(0, max_id)
  int_cnt = rep(0, max_id)
  
  for (j in 1:length(net[[2]])){
    reg = net[[2]][[j]]
    
    for (i in 1:nrow(reg)){
      target_id = as.integer(rownames(reg)[i])
      if ( target_id > 0){
        mi_val = reg[i, "MI"]
        mean_mi[target_id] = mean_mi[target_id] + mi_val
        int_cnt[target_id] = int_cnt[target_id] + 1
      }
    }
  }
  
  for (i in 1:length(mean_mi))
    if (int_cnt[i] != 0)
      mean_mi[i] = mean_mi[i] / int_cnt[i]
  
  mean_mi = mean_mi[as.integer(rownames(mat))]
  mean_exp = mean_exp[as.integer(rownames(mat))]
  names(mean_mi) = names(mean_exp) = rownames(mat)
  mean_mi <- subset(mean_mi,mean_mi>0)
  
  final_mat <- merge(mean_mi,mat[,'expMean'],by='row.names')
  rownames(final_mat) <- final_mat[,'Row.names']
  final_mat[,'Row.names'] <- NULL
  colnames(final_mat) <- c('mean_mi','meanExp')
  
  return(final_mat)
}

############  WRITE VST EXP CORRELATION MATRIX ######################
#call meanExpMIcorrelation for all tissues in varNames
for(tissue in varNames){  
  if( tissue %in% varNames[grep('tcga', varNames)] ){      #iterate over tcga tissues only
    print(paste0('vstexp_mi_',tissue))
    
    assign(paste0('vstexp_mi_',tissue),     #var name vstexp_mi + tissue name such as AdiposeSub etc
           meanExpMIcorrelation(
             get(paste0('expmat_', tissue) ),
             get(tissue)  )   )    #assign median of mutual information and expression to target gene
  }
  else{          #iterate over gtex tissues only
    print(paste0('vstexp_mi_gtex_',tissue))
    assign(paste0('vstexp_mi_gtex_',tissue),
           meanExpMIcorrelation(
             get(paste0( 'expmat_gtex_', tissue)  ),
             get(tissue)  )    )
  }
} 

#tissue wise correlation matrix: correlation between medianMI and medainExp of target gene
vst_exp_mi_tissues <- NULL
vst_exp_mi_tissues <- ls()[grep('vstexp_mi', ls())] #extract all tissue matrix (mi, exp)

miVstExpCorMatrix <- NULL
miVstExpCorMatrix<-matrix(ncol = 3, nrow=0)
for( tissue in  vst_exp_mi_tissues){  #iterate over each vstexp_mi object
  miVstExpCorMatrix<-rbind(miVstExpCorMatrix,c(tissue,
                                               cor.test(as.numeric(get(tissue)$mean_mi),
                                                        as.numeric(get(tissue)$meanExp))$estimate,
                                               cor.test(as.numeric(get(tissue)$mean_mi),
                                                        as.numeric(get(tissue)$meanExp))$p.value )  )
}
colnames(miVstExpCorMatrix)<-c('tissue','estimate','pval')#colnames of matrix
#write vstmat_expression and meanMI correlation
write.table(miVstExpCorMatrix,file='miVstExpCorMatrix.txt', sep = "\t", quote = FALSE, row.names = FALSE)

########### WRITE TPM EXP & MI CORRELATION MATRIX #################
for(tissue in names(pancan_tpm_count)){
  print(paste0('tpmexp_mi_',tissue))
  #print(head(pancan_tpm_count[[tissue]]))
  #break
  assign(paste0('tpmexp_mi_',tissue),
         meanExpMIcorrelation(
           pancan_tpm_count[[tissue]],
           get(tissue)
         ))
}
#call exp_mi_info for all tissues in varNames

#tissue wise correlation matrix: correlation between medianMI and medainExp of target gene
tpm_exp_mi_tissues <- NULL
tpm_exp_mi_tissues <- ls()[grep('tpmexp_mi', ls())] #extract all tissue matrix (mi, exp)

miTpmExpCorMatrix <- NULL
miTpmExpCorMatrix<-matrix(ncol = 3, nrow=0)

for( tissue in  tpm_exp_mi_tissues){  #iterate over each tpmexp_mi object
  miTpmExpCorMatrix<-rbind(miTpmExpCorMatrix,c(tissue,
                                               cor.test(get(tissue)[,'mean_mi'],
                                                        get(tissue)[,'meanExp'])$estimate,
                                               cor.test(get(tissue)[,'mean_mi'],
                                                        get(tissue)[,'meanExp'])$p.value) )  
}

colnames(miTpmExpCorMatrix)<-c('tissue','estimate','pval')#colnames of matrix
#write tpmmat_expression and meanMI correlation
write.table(miTpmExpCorMatrix,file='miTpmExpCorMatrix.txt', sep = "\t", quote = FALSE, row.names = FALSE)



save(vstexp_mi_gtex_AdiposeSub,
     vstexp_mi_gtex_AdiposeVis,
     vstexp_mi_gtex_Adrenal,
     vstexp_mi_gtex_ArteryAor,
     vstexp_mi_gtex_ArteryCor,
     vstexp_mi_gtex_ArteryTib,
     vstexp_mi_gtex_BrainCau,
     vstexp_mi_gtex_BrainCer,
     vstexp_mi_gtex_BrainCerHem,
     vstexp_mi_gtex_BrainCor,
     vstexp_mi_gtex_BrainFroCor,
     vstexp_mi_gtex_BrainNuc,
     vstexp_mi_gtex_Breast,
     vstexp_mi_gtex_CellsEBV,
     vstexp_mi_gtex_CellsTra,
     vstexp_mi_gtex_ColonSig,
     vstexp_mi_gtex_ColonTra,
     vstexp_mi_gtex_EsophagusGasJun,
     vstexp_mi_gtex_EsophagusMuc,
     vstexp_mi_gtex_EsophagusMus,
     vstexp_mi_gtex_HeartAtrApp,
     vstexp_mi_gtex_HeartLefVen,
     vstexp_mi_gtex_Liver,
     vstexp_mi_gtex_Lung,
     vstexp_mi_gtex_Muscle,
     vstexp_mi_gtex_Nerve,
     vstexp_mi_gtex_Pancreas,
     vstexp_mi_gtex_Pituitary,
     vstexp_mi_gtex_Prostate,
     vstexp_mi_gtex_SkinNotSun,
     vstexp_mi_gtex_SkinSunExp,
     vstexp_mi_gtex_Spleen,
     vstexp_mi_gtex_Stomach,
     vstexp_mi_gtex_Testis,
     vstexp_mi_gtex_Thyroid,
     vstexp_mi_gtex_WholeBlood,
     vstexp_mi_tcga_blca,
     vstexp_mi_tcga_brca,
     vstexp_mi_tcga_cesc,
     vstexp_mi_tcga_coad,
     vstexp_mi_tcga_esca,
     vstexp_mi_tcga_gbm,
     vstexp_mi_tcga_hnsc,
     vstexp_mi_tcga_kirc,
     vstexp_mi_tcga_kirp,
     vstexp_mi_tcga_laml,
     vstexp_mi_tcga_lgg,
     vstexp_mi_tcga_lihc,
     vstexp_mi_tcga_luad,
     vstexp_mi_tcga_lusc,
     vstexp_mi_tcga_ov,
     vstexp_mi_tcga_paad,
     vstexp_mi_tcga_pcpg,
     vstexp_mi_tcga_prad,
     vstexp_mi_tcga_read,
     vstexp_mi_tcga_sarc,
     vstexp_mi_tcga_skcm,
     vstexp_mi_tcga_stad,
     vstexp_mi_tcga_tgct,
     vstexp_mi_tcga_thca,
     vstexp_mi_tcga_thym,
     vstexp_mi_tcga_ucec,
     tpmexp_mi_AdiposeSub,
     tpmexp_mi_AdiposeVis,
     tpmexp_mi_Adrenal,
     tpmexp_mi_ArteryAor,
     tpmexp_mi_ArteryCor,
     tpmexp_mi_ArteryTib,
     tpmexp_mi_BrainCau,
     tpmexp_mi_BrainCer,
     tpmexp_mi_BrainCerHem,
     tpmexp_mi_BrainCor,
     tpmexp_mi_BrainFroCor,
     tpmexp_mi_BrainNuc,
     tpmexp_mi_Breast,
     tpmexp_mi_CellsEBV,
     tpmexp_mi_CellsTra,
     tpmexp_mi_ColonSig,
     tpmexp_mi_ColonTra,
     tpmexp_mi_EsophagusGasJun,
     tpmexp_mi_EsophagusMuc,
     tpmexp_mi_EsophagusMus,
     tpmexp_mi_HeartAtrApp,
     tpmexp_mi_HeartLefVen,
     tpmexp_mi_Liver,
     tpmexp_mi_Lung,
     tpmexp_mi_Muscle,
     tpmexp_mi_Nerve,
     tpmexp_mi_Pancreas,
     tpmexp_mi_Pituitary,
     tpmexp_mi_Prostate,
     tpmexp_mi_SkinNotSun,
     tpmexp_mi_SkinSunExp,
     tpmexp_mi_Spleen,
     tpmexp_mi_Stomach,
     tpmexp_mi_tcga_blca,
     tpmexp_mi_tcga_brca,
     tpmexp_mi_tcga_cesc,
     tpmexp_mi_tcga_coad,
     tpmexp_mi_tcga_esca,
     tpmexp_mi_tcga_gbm,
     tpmexp_mi_tcga_hnsc,
     tpmexp_mi_tcga_kirc,
     tpmexp_mi_tcga_kirp,
     tpmexp_mi_tcga_laml,
     tpmexp_mi_tcga_lgg,
     tpmexp_mi_tcga_lihc,
     tpmexp_mi_tcga_luad,
     tpmexp_mi_tcga_lusc,
     tpmexp_mi_tcga_ov,
     tpmexp_mi_tcga_paad,
     tpmexp_mi_tcga_pcpg,
     tpmexp_mi_tcga_prad,
     tpmexp_mi_tcga_read,
     tpmexp_mi_tcga_sarc,
     tpmexp_mi_tcga_skcm,
     tpmexp_mi_tcga_stad,
     tpmexp_mi_tcga_tgct,
     tpmexp_mi_tcga_thca,
     tpmexp_mi_tcga_thym,
     tpmexp_mi_tcga_ucec,
     tpmexp_mi_Testis,
     tpmexp_mi_Thyroid,
     tpmexp_mi_WholeBlood,file='meanExp_meanMI_correlation.Rda')
