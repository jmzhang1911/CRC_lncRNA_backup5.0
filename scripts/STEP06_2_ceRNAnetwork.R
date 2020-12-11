rm(list = ls())
source('Utils.R')


load('outcomes/inputdata/input.RData')
load('outcomes/candidate/maxtrans_df.RData')
load('outcomes/ceRNAnetwork/lncRNA_miRNA_res/miRNA_inter_reslist.RData')




#===> miRNA target genes by TargetScan db 
# mmu_scan:all the miRNA target genes
if(T){
  mmu_scan <- read.table('data/db/Predicted_Targets_Info.txt', sep = '\t', header = T) %>%
    mutate(miR.Family = str_replace_all(miR.Family, '/', '='),
           miR.Family = str_replace_all(miR.Family, '\\.', '_'),
           PCT = as.numeric(PCT)) %>%
    dplyr::filter(PCT != 'NULL',
                  PCT >= 0.2) %>%
    dplyr::select(miRNA_family = miR.Family,
                  target_gene_id = Gene.ID)
  #save(mmu_scan, file = 'outcomes/ceRNAnetwork/mmu_scan.RData')
}



#===> miRNAs interact with lncRNA
load('outcomes/ceRNAnetwork/mmu_scan.RData')
load('outcomes/candidate/crc_corgenelist.RData')
load('outcomes/candidate/maxtrans_df.RData')
# filter the miRNAs which have 3 interactions with any of lncRNA at least
prob <- maxtrans_df %>% dplyr::filter(gene_id %in% crc_corgenelist$cor_lnc) %>%
  pull(trans_id) #737

miRNA_inter_final <- miRNA_inter_reslist$miRNA_inter_inte %>%
  mutate(miRNA_family = str_replace_all(miRNA_family, '/', '='),
         miRNA_family = str_replace_all(miRNA_family, '\\.', '_')) %>%
  dplyr::filter(loop == 1, lnc_trans %in% prob) %>%
  group_by(lnc_trans, miRNA_family) %>% summarise(inter_num = n()) %>%
  dplyr::filter(inter_num >= 3)


# 89 candidate miRNAs, next will analysis survival and DE of them 
candidate_miRNA <- miRNA_inter_final %>% 
  group_by(miRNA_family) %>% arrange(miRNA_family, desc(inter_num)) %>%
  top_n(1, inter_num) %>% distinct(miRNA_family, inter_num) %>% 
  rename(maxinternum = inter_num) %>%
  mutate(prob = str_extract(miRNA_family, '^[a-zA-Z]{3,}-[0-9]+')) 

save(miRNA_inter_final, candidate_miRNA, 
     file = 'outcomes/ceRNAnetwork/lncRNA_miRNA_res/lncRNA_miRNA_res.RData')

















  

#===> miRNA target genes enrichment
if(T){
  miRNA_ve <- unique(miRNA_inter_final$miRNA_family)
  scan_anno_d1 <- mmu_scan %>% 
    dplyr::filter(miRNA_family %in% miRNA_ve) 
  table(scan_anno_d1$miRNA_family %in% miRNA_inter_final$miRNA_family)
  
  for(i in miRNA_ve){
    i = 'miR-181-5p'
    tmp <- dplyr::filter(scan_anno_d1, miRNA_family == i) %>% 
      distinct(target_gene_id) %>%
      pull(target_gene_id) %>%
      str_replace_all('\\.[0-9]','')
    if(length(tmp)>0){
      df <- mygetenrichedgenefromkk2(tmp)
      t1 <- paste('outcomes/inte_miRNA/CRC_miRNA2/', i, sep = '')
      t2 <- paste('outcomes/inte_miRNA/not_CRC_miRNA2/', i, sep = '')
      title <- if_else(length(df$ENSEMBL)>0, t1, t2)
      myenrich(tmp, title) 
    }
  }
}


input_matrix_count$mRNA_lncRNA_count$gene_id

load('outcomes/DEanalysis/DE_results_inte.RData')
df <- DE_results_inte$DE_each_inte















