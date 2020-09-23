#########################################################################
# Author : Desiree Wilson
# Date   : August 24, 2020
# Purpose: The purpose of this code is to consolidate
#          all probes associated with one gene into 
#          one representative value for the gene in the 
#          GSE65904 dataset. The gene expression 
#          was measured using the following chip:
#          Illlumina Human HT-12V4.0 BeadChip arrays
#          For more info, look at this link:
#          https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65904
#
########################################################################

#loading in libraries:
library(GEOquery)
library(beadarray)
library(illuminaHumanv4.db)
library(readr)
library(ggplot2)
library(WGCNA)


#another way to obtain GEO dataset; directly from website:
#reference: Page31 of following document:
#https://www.bioconductor.org/packages/release/bioc/vignettes/beadarray/inst/doc/beadsummary.pdf
url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE65nnn/GSE65904/matrix/"
filenm <- "GSE65904_series_matrix.txt.gz"
if(!file.exists("GSE65904_series_matrix.txt.gz")) download.file(paste(url, filenm, sep=""), destfile=filenm)
gse <- getGEO(filename=filenm)
head(exprs(gse))


#corresponding feature data according to the same reference listed above:
summaryData <- as(gse, "ExpressionSetIllumina")
head(summaryData)


#using WGCNA to collapse the probe expression data to gene level data:
datET <- exprs(gse)
rowGroup <- fData(summaryData)$SYMBOL
rowID <- rownames(datET)

datasetMaxMean <- collapseRows(datET = datET, rowID = rowID
                               ,rowGroup = rowGroup, method = "MaxMean")
datasetAverage <- collapseRows(datET = datET, rowID = rowID
                               ,rowGroup = rowGroup, method = "Average")
datasetMaxVariance <- collapseRows(datET = datET, rowID = rowID
                                   ,rowGroup = rowGroup, method = "maxRowVariance")

#writing the pheno data to file:
pData <- pData(gse)
directory <- c("C:\\Users\\wilso\\Documents\\Bioinformatics\\BioinformaticsProjects\\KatieHinchee\\MelatoninA\\Results")
write.table(pData, file = paste0(directory,"\\GSE65904_phenoData.tsv")
            ,append = FALSE
            ,quote = FALSE
            ,sep = "\t"
            ,row.names = TRUE
            ,col.names = TRUE)


#writing the collapsed row data to a file:
directory <- c("C:\\Users\\wilso\\Documents\\Bioinformatics\\BioinformaticsProjects\\KatieHinchee\\MelatoninA\\Results")

write.table(datasetMaxMean$datETcollapsed, file = paste0(directory,"\\GSE65904_maxMean.tsv")
            ,append = FALSE
            ,quote = FALSE
            ,sep = "\t"
            ,row.names = TRUE
            ,col.names = TRUE)

write.table(datasetAverage$datETcollapsed, file = paste0(directory,"\\GSE65904_average.tsv")
            ,append = FALSE
            ,quote = FALSE
            ,sep = "\t"
            ,row.names = TRUE
            ,col.names = TRUE)

write.table(datasetMaxVariance$datETcollapsed, file = paste0(directory,"\\GSE65904_maxVariance.tsv")
            ,append = FALSE
            ,quote = FALSE
            ,sep = "\t"
            ,row.names = TRUE
            ,col.names = TRUE)






