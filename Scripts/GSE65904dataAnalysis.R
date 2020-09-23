#########################################################################
# Author : Desiree Wilson
# Date   : April 26, 2020
# Purpose: The purpose of this code is relationship
#          to see if there is any relationship
#          between melatonin-A gene expression and
#          sex in metastatic melanoma in the 
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


#setting working directory:
setwd("C:/Users/wilso/Documents/Bioinformatics/BioinformaticsProjects/KatieHinchee/MelatoninA")


#reading in gene expression dataset:
dat_raw <- read_tsv(file = 'Data/GSE65904_series_matrix.txt'
                    #, skip = 64
                    , comment = '!')

dim(dat_raw)
tail(dat_raw)[,1:3]

#reading in the pheno data from gene expression dataset:
pheno_raw <- read_tsv(file = 'Data/GSE65904_series_matrix.txt'
                      , skip = 26
                      , n_max = 39)


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
summaryData
head(fData(summaryData))
head(pData(summaryData))



