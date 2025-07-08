############################################################
##    reading the POCROC Proteomics data set              ##
############################################################

#install.packages("tidyverse")
library(tidyverse)

prot_data <- read.csv("./Data_for_paper_figures/prot_quant_fragm_filt_msstats_norepl_nopool_71722.csv")

prot_data <- as.data.frame(prot_data[-1,-1])
prot_data <- separate(prot_data, Protein, into = c("random1", "random2", "Protein"), sep = "\\|")
prot_data <- subset(prot_data, select = -c(random1, random2))
rownames(prot_data) <- NULL
prot_data <- column_to_rownames(prot_data, "Protein")

coldata_prot <- read.table("./Data_for_paper_figures/coldata_prot_newHRD_newReclassify_71822.txt", sep = "\t", header = T)

prot_data <- prot_data[,rownames(coldata_prot)]

all(colnames(prot_data) == rownames(coldata_prot))


prot_data_bpca <- read.csv("./Data_for_paper_figures/prot_data_bpca_71722.csv")
prot_data_bpca <- column_to_rownames(prot_data_bpca, "X")

all(colnames(prot_data) == colnames(prot_data_bpca))

################## Putting 48986 in the HRP category
### SP57 and SP76 go to HRP category, ignoring their somatic mutation status

coldata_prot$HRD_HRP[grep("48986", coldata_prot$Patient_ID)] <- c("HRP", "HRP", "HRP")

all(colnames(prot_data) == colnames(prot_data_bpca))
all(colnames(prot_data) == rownames(coldata_prot))