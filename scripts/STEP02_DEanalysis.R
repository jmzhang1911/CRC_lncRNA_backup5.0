rm(list = ls())
source('Utils.R')


mkdir('outcomes/DEanalysis/')
load('outcomes/inputdata/input.RData')
#===> 批次
gene_count <- input_matrix_count$mRNA_lncRNA_count %>%
  column_to_rownames(var = 'gene_id')
boxplot(log2(gene_count+1))



#===> DESeq2 analysis
if(T){
  condition <- factor(c(rep("control", 3), rep("week2", 3),rep("week4", 3),
                        rep("week7", 3), rep("week10", 3)),
                      levels = c("control", "week2","week4","week7","week10"))
  colData <- data.frame(row.names = colnames(gene_count), condition)
  dds <- DESeqDataSetFromMatrix(gene_count, colData, design= ~ condition)
  dds <- DESeq(dds)
}

if(T){
  tmp <- myddnor(gene_count) #标准化后举证
  tmp2 <- as.data.frame(tmp) %>% rownames_to_column(var = 'gene_id') %>%
    gather(key = 'sample', value = 'count', 2:16)
  p <- ggplot(tmp2, aes(x = sample, y = log10(count+1))) +
    geom_boxplot(fill = '#228B22') +
    theme_classic() +
    labs(x = '',
         title = 'Batch effect after normalization') +
    theme(axis.text.x = element_text(hjust = 1,
                                     vjust = 1,
                                     angle = 45)) 
  ggsave(p, filename = 'outcomes/DEanalysis/Batch_effect.png')
}

DE_res_list <- list(
  wek2_ctrl = results(dds, contrast=c("condition", "week2", "control")),
  wek4_ctrl = results(dds, contrast=c("condition", "week4", "control")),
  wek7_ctrl = results(dds, contrast=c("condition", "week7", "control")),
  wek10_ctrl = results(dds, contrast=c("condition", "week10", "control")),
  wek4_wek2 = results(dds, contrast=c("condition", "week4", "week2")),
  wek7_wek2 = results(dds, contrast=c("condition", "week7", "week2")),
  wek10_wek2 = results(dds, contrast=c("condition", "week10", "week2")),
  wek7_wek4 = results(dds, contrast=c("condition", "week7", "week4")),
  wek10_wek4 = results(dds, contrast=c("condition", "week10", "week4")),
  wek10_wek7 = results(dds, contrast=c("condition", "week10", "week7")))


DE_results_inte <- mydeanalysis_inte(DE_res_list, 
                                     input_matrix_count$mRNA_lncRNA_count)
save(DE_results_inte, file = 'outcomes/DEanalysis/DE_results_inte.RData')
















