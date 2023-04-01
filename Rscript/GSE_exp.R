if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("GEOquery")

#Data Analysis
#https://sbc.shef.ac.uk/geo_tutorial/tutorial.nb.html
my_id <- "GSE29532"
my_id <- "GSE62646"
my_id <- "GSE66360"
gse <- getGEO(my_id)

gse <- gse[[1]]

exprs(gse) <- log2(exprs(gse))
annot <- fData(gse)#extract gene symbol for probeset
sample <- pData(gse)#extract sample annotation
exp <- exprs(gse)


##Columns are different for different GPL annotation. First check the columns and then decide to extract
colnames(annot)
View(annot)
#extract two columns; 'ID' and 'gene_assignment
#annot <- annot[,c('ID','gene_assignment')]
annot <- annot[,c('ID','Gene Symbol')]
#annot <- annot[,c('ID','Symbol')]
final_exp <- merge(annot,exp,by="row.names")

write.table(final_exp,file = paste0('/Users/hs3290/Documents/extra/AMI/',my_id,'_exp.txt'),sep="\t",row.names = FALSE)

#OPEN ALL THE EXPRESSION FILES IN EXCEL AND ADD SYMBOL COLUMN BY SPLITTING THE GENE SYMBOL OR GENE_ASSIGNMENT COLUMNS
#AFTER THAT MERGE ALL THE THREE FILES ON THE BASIS OF SYMBOL COLUMNS
df1 <- read.csv('/Users/hs3290/Documents/extra/AMI/GSE29532_exp.txt',sep="\t",header=TRUE)
df2 <- read.csv('/Users/hs3290/Documents/extra/AMI/GSE62646_exp.txt',sep="\t",header=TRUE)
df3 <- read.csv('/Users/hs3290/Documents/extra/AMI/GSE66360_exp.txt',sep="\t",header=TRUE)

df <- merge(df1, df2, by ="Symbol")
View(df)
