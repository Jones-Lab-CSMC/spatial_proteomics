# MSstats quantification
df <- read.csv("fragment_NOIMP_reformat_corrected_inputfor_msstats_SDF2_71722.csv")

#BiocManager::install("MSstats")

library(MSstats)
library(tidyverse)

df_msstats <- dataProcess(df, normalization=FALSE, MBimpute=FALSE)

#save(df_msstats, file= "msstats_dataprocess_TvsS_mapdia_filt_SDF1_NOIMP_40221.RData")

sampleQuant_df <- quantification(df_msstats)

write.csv(sampleQuant_df, "prot_quant_fragm_filt_msstats_norepl_nopool_71722.csv")


sampleQuant_df <- as.data.frame(sampleQuant_df[-1,])
sampleQuant_df <- separate(sampleQuant_df, Protein, into = c("random1", "random2", "Protein"), sep = "\\|")
sampleQuant_df <- subset(sampleQuant_df, select = -c(random1, random2))

## saving BPCA imputed data set
set.seed(4242)
data <- pcaMethods::pca(as.matrix(prot_data), nPcs = 10, method = "bpca")
data_df <- as.data.frame(pcaMethods::completeObs(data))

write.csv(data_df,"prot_data_bpca_71722.csv")

write.table(coldata_prot, "coldata_prot_newHRD_newReclassify_71822.txt", quote = F, sep = "\t", row.names = T, col.names = T)
