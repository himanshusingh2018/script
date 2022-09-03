#################################################################################################################################
#   File name : triplex.hg38.R
#   Author    : Himanshu Singh
#   Date      : August 11, 2021
#   Email     : himanshu720@gmail.com
#   Purpose   : Extract triplex motifs from the genome coordinates
#   Instruction:
#               1. Source("triplex.hg38.R")
#               2. Run the function "triplex_motifs(chrNo,start,end,outputFile)"
#               3. The data will be saved in the folder "outputFile" Folder"
#               4. Criteria: min_score=17,min_len=10, max_len = 15, p_value = 0.001
#               5. Run the function "triplex_motif_gff_to_bed(triplexmotifOUT,chrSTART)
#               6. The BED file will be saved in the folder "outputFile" Folder
#               7. Instructions is in README()
#
#   Example:
#               > source("triplex.hg38.R")
#               > triplex_motifs('chr9',119166629,119370435,'dbc1')
#               > ReadME()


#### Load Libraries ----
if(!require(BSgenome.Hsapiens.UCSC.hg38, quietly = TRUE)) {
  if(!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
  require(BSgenome.Hsapiens.UCSC.hg38, quietly = TRUE)
}
if(!require(triplex, quietly = TRUE)) {
  if(!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("triplex")
  require(triplex, quietly = TRUE)
}
if(!require(rtracklayer, quietly = TRUE)) {
  if(!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("rtracklayer")
  require(rtracklayer, quietly = TRUE)
}
if(!require(bedr, quietly = TRUE)) {
  install.packages("bedr")
  require(bedr, quietly = TRUE)
}

#### Extract Triplex motif regions from the hg38 human genome ----
triplex_motifs <- function(chrNo,start,end,outputFile){
  dir.create(outputFile)
  t <- triplex.search(Hsapiens[[chrNo]][start:end],min_score=17,min_len=10, max_len = 15, p_value = 0.001)
  export(as(t, "GRanges"),paste0(outputFile,'/',outputFile,'.gff'), version="3")
  #Replace 'chr1' to input chrNo
  temp <- gsub(pattern = 'chr1', replace = chrNo, x = readLines(paste0(outputFile,'/',outputFile,'.gff'),encoding = 'UTF-8'))
  writeLines(temp, con=paste0(outputFile,'/',outputFile,'.gff'))
  rm(temp)#remove temp object
  print(paste0('Intra-molecular triplex motif coordinates: ',outputFile,'/',outputFile,'.gff'))
  
  writeXStringSet(as(t, "DNAStringSet"), file=paste0(outputFile,'/',outputFile,".fa"), format="fasta")
  print(paste0('Intra-molecular triplex motif sequence: ',outputFile,'/',outputFile,'.fa'))
  
  #generate all the triplex diagram
  for(i in 1:length(t)){
    jpeg(filename = paste0(outputFile,'/triplex.',i,'.jpeg'), res = 50)#, width = 120, height = 120)
    triplex.diagram(t[i])
    dev.off()
  }
}

#### Generate BED file triplex motif from the output of triplex_motifs function generated file '.gff'
triplex_motif_gff_to_bed <- function(triplexmotifOUT,chrSTART){
  triplex_motif <- read.table(triplexmotifOUT, skip = 3)[,c(1,4,5)]
  triplex_motif$V4 <- triplex_motif$V4 + chrSTART
  triplex_motif$V5 <- triplex_motif$V5 + chrSTART
  #writeLines(triplex_motif, 
  #           con=paste0(substr(triplexmotifOUT,1,nchar(triplexmotifOUT)-4),'.bed'))
  write.table(x = triplex_motif,
              file = paste0(substr(triplexmotifOUT,1,nchar(triplexmotifOUT)-4),'triplex.motif.bed'),
              sep="\t",
              col.names = FALSE,
              row.names = FALSE,
              quote = FALSE)
}

#### Identify triplex motifs overlapping or Nearby cis-Regulatory Regions ----
#triplex_motif_cCREs_overlapping(triplexmotifOUT,cCREs){
  #bedtools path: /Users/hs3290/miniconda3/bin/bedtools
#}

#### set README function ----
README <- function(){
  writeLines(paste(rep("#", 100), collapse = ""))
  writeLines('"triplex_motifs": function will extract the triplex motifs within USER PROVIDED hg38 human genome')
  writeLines('  User Input:\n\tchromosome number e.g. "chrX"\n\tStart e.g. 1\n\tEnd e.g. 1000, outputFile e.g. "dbc1"')
  writeLines('  Criteria: min_score=17,min_len=10, max_len = 15, p_value = 0.001')
  writeLines('  Example:\n\ttriplex_motifs("chr9",119166629,119370435,"dbc1")')
  writeLines('\n"triplex_motif_gff_to_bed": function will convert triplex motif output gff file coordinate to BED format')
  writeLines('  User Input: \n\ttriplex motif GFF file e.g. "dbc1/dbc1.gff" \n\tChr Start Pos e.g. chrSTART=119166629)')
  writeLines('  Example:\n\ttriplex_motif_gff_to_bed(triplexmotifOUT = "dbc1/dbc1.gff",chrSTART=119166629)')
  writeLines(paste(rep("#", 100), collapse = ""))
}
