source('Utils.R')
library(tidyverse)
library(survival)
library(survminer)


load('outcomes/ceRNAnetwork/survival/Survival_data.RData')
mkdir('outcomes/ceRNAnetwork/survival/miRNA_survival')


mysurval <- function(x){
  file = basename(x)
  miRNA_expr_tmp <- condi_matrix[rownames(miRNA_expr) == file,]
  miRNA_expr_tmp <- t(miRNA_expr_tmp)
  ve <- miRNA_expr_tmp[,1]
  df_sur <- dplyr::mutate(df, GENE_condi = ve[bcr_patient_uuid])
  df_sur <- df_sur[!is.na(df_sur$GENE_condi),]
  fit <- survfit(Surv(os, vital_status) ~ GENE_condi, data = df_sur)
  title = paste('Survival analysis of', file)
  p <- ggsurvplot(fit,
                  data = df_sur,
                  title = title,
                  palette = c("#00008B", "#8B0000"),
                  surv.median.line = "hv",
                  pval = T)
  return(p)}



miRNAlist <- rownames(miRNA_expr)

for(i in miRNAlist){
  try(dev.off())
  title <- paste('outcomes/ceRNAnetwork/survival/miRNA_survival/', i, sep = '')
  p <- mysurval(title)
  png(paste(title, '_survival_analysis.png'))
  print(p)
  dev.off()
}



