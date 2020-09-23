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
library(ggplot2)


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



#I just downloaded the original manifest file from the Illumina website:
#https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/humanht-12/humanht-12_v4_0_r2_15002873_b.txt.zip
#There is exactly 1/ONE probe for melatonin 1A (genesymbol: MTNR1A). So from here
#I am going to do a t-test:

#identifying the location of the illumina probe:
gs_mtnr1ag_idx <- grep("MTNR1A", fData(summaryData)$ILMN_Gene)
probe_mtnr1a_name <- row.names(fData(summaryData))[gs_mtnr1ag_idx]


#pulling out the illumina probe in the gene expression data:
exprs_probe_mtnr1a_idx <- grep(probe_mtnr1a_name, row.names(exprs(summaryData)))
exprs_probe_mtnr1a_dat <- exprs(summaryData)[exprs_probe_mtnr1a_idx,]


#identifying the samples that are missing sex information:
missing_sex_idx <- which(pData(summaryData)$`gender:ch1` %in% c("NA"))


#need to drop this sample in both the pheno and expression data:
exprs_probe_mtnr1a_dat_filtered <- exprs_probe_mtnr1a_dat[-missing_sex_idx]
sex_filtered <- pData(summaryData)$`gender:ch1`[-missing_sex_idx]


#finally performing the t-test
t_test <- t.test(exprs_probe_mtnr1a_dat_filtered ~ sex_filtered)
summary(t_test)

#plotting the expression values:
normalized_exprs_ln <- exprs_probe_mtnr1a_dat_filtered
sex <- sex_filtered
mtnr1a_exprs_df <- data.frame(normalized_exprs_ln, sex)

exprs_boxplot <- ggplot(mtnr1a_exprs_df, aes(x = sex, y = normalized_exprs_ln)) + 
  geom_boxplot() +
  labs(title = "MTNR1 Normalized Expression"
       ,subtitle = paste0("p.value = ", t_test$p.value)
  )

#printing to file:
pdf_file_name <-'Results/mtnr1a.pdf'
pdf(file = pdf_file_name
    ,height = 8.5
    ,width = 11)
exprs_boxplot
dev.off()

