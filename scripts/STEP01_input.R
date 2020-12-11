source('Utils.R')


#===> rename matrix
gene_count <- read.table("data/lncRNA-mRNA-kaenno.counts.matrix",
                         header = T, row.names = 1)
gene_count <- rownames_to_column(gene_count, var = 'gene') %>%
  dplyr::select(gene_id = gene, control_1 = X823.feature, control_2 = X824.feature,
                control_3 = X825.feature, week2_1 = X826.feature,
                week2_2 = X828.feature, week2_3 = X829.feature,
                week4_1 = X830.feature, week4_2 = X831.feature,
                week4_3 = X832.feature, week7_1 = X835.feature,
                week7_2 = X836.feature, week7_3 = X839.feature,
                week10_1 = X840.feature, week10_2 = X841.feature,
                week10_3 = X847.feature)



#===> lnc count
lnc_genelist <- read.table('data/final_lnc.bed')$V4
length(unique(lnc_genelist))
lncRNA_count <- dplyr::filter(gene_count,
                              gene_id %in% lnc_genelist) %>%
  mutate(sum = rowSums(.[2:16])) %>%
  dplyr::filter(sum > 5) %>%
  dplyr::select(-sum)



#===> mRNA count
mRNA_count <- dplyr::filter(gene_count,
                            gene_id %!in% lnc_genelist) %>%
  mutate(sum = rowSums(.[2:16])) %>%
  dplyr::filter(sum > 5) %>%
  dplyr::select(-sum)



#===> mRNA_lncRNA count
mRNA_lncRNA_count <- gene_count %>%
  mutate(sum = rowSums(.[2:16])) %>%
  dplyr::filter(sum > 5) %>%
  dplyr::select(-sum)



#===> rename TMM matrix
gene_TMM <- read.table("data/lncRNA-mRNA-kaenno.TMM.EXPR.matrix",
                       header = T, row.names = 1)
colnames(gene_TMM) <- colnames(gene_count)[-1]



#===> rename fpkm
fpkm_count <- read.table("data/merged_fpkm_table.txt",
                         header= T, stringsAsFactors = F) %>%
  dplyr::select(gene_id = gene_id,control_1 = X823.fpkm.count,
                control_2 = X824.fpkm.count, control_3 = X825.fpkm.count,
                week2_1 = X826.fpkm.count, week2_2 =X828.fpkm.count,
                week2_3 = X829.fpkm.count, week4_1 = X830.fpkm.count,
                week4_2 = X831.fpkm.count, week4_3 = X832.fpkm.count,
                week7_1 = X835.fpkm.count, week7_2 =X836.fpkm.count,
                week7_3 = X839.fpkm.count, week10_1 = X840.fpkm.count,
                week10_2 = X841.fpkm.count, week10_3 = X847.fpkm.count)

#===> read exon trans
exon_trans <- read.table("data/count_results.txt") %>%
  dplyr::select(gene_id = V1, trans_id = V2, exon_num = V3, translen = V4)


#===> save
input_matrix_count <- list(
  'mRNA_lncRNA_count' = mRNA_lncRNA_count,
  'lnc_genelist' = lnc_genelist,
  'gene_count' = gene_count,
  'mRNA_count' = mRNA_count,
  'lncRNA_count' = lncRNA_count, 
  'gene_TMM' = gene_TMM,
  'fpkm_count' = fpkm_count,
  'exon_trans' = exon_trans)

mkdir('outcomes/inputdata')
save(input_matrix_count, file = 'outcomes/inputdata/input.RData')


