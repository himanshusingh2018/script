#load vst_norm interactmes
load('/Volumes/af_lab_data/hs3290/mi_exp_cor/All_62_ARACNE.rda')
#print header of table
cat('tissueName','vst_norm_VS_tpm_count','vst_norm_VS_norm_count','tpm_count_VS_norm_count\n',sep = "\t")
for(i in setdiff( varNames, c('tcga_kirc','ColonSig') ) ){
  cat(paste(i,#tissue name
            #vst_norm vs tpm_count
            countSharedInteractions(i,
                                    load(paste0('/Volumes/af_lab_data/hs3290/aracne_regen_xena/tpm_count/',
                                                i,
                                                '/net.interactome.rda'))
                                    ),
            countSharedInteractions(i,
                                    load(paste0('/Volumes/af_lab_data/hs3290/aracne_regen_xena/norm_count/',
                                                i,
                                                '/net.interactome.rda'))
                                    ),
            #tpm_count vs norm_count
            countSharedInteractions(load(paste0('/Volumes/af_lab_data/hs3290/aracne_regen_xena/tpm_count/',
                                                i,
                                                '/net.interactome.rda')),
                                    load(paste0('/Volumes/af_lab_data/hs3290/aracne_regen_xena/norm_count/',
                                                i,
                                                '/net.interactome.rda'))),
            "\n"),sep="\t")
  #break
}

