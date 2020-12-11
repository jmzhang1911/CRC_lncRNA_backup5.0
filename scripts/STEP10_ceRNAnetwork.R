source('Utils.R')

load('outcomes/inputdata/input.RData')
load('outcomes/inte_miRNA/miRNA_inter_reslist.RData')
load('outcomes/inte_miRNA/mmu_scan.RData')
gene2trans <- input_matrix_count$exon_trans
lnclist <- input_matrix_count$lnc_genelist

# targetsacan miRNA family
cemiRNA_family <- c('miR-181-5p', 'miR-128-3p', 'miR-425-5p=489-3p',
                    'miR-27-3p', 'miR-183-5p', 'miR-140-3p_2=497',
                    'miR-141-3p=200a-3p', 'miR-212-5p', 'miR-7-5p')

# MirTarBase
cemiRNA <- c('mmu-miR-181a-5p', 'mmu-miR-128-3p', 'mmu-miR-425-5p',
             'mmu-miR-27b-3p', 'mmu-miR-183-5p', 'mmu-miR-140-3p',
             'mmu-miR-141-3p', 'mmu-miR-212-5p', 'mmu-miR-7a-5p')
cemRNAlist <- str_remove_all(cemiRNA, 'mmu-')



#===> miRNAs interact with lncRNA
miRNA_inter_final <- miRNA_inter_reslist$miRNA_inter_inte %>%
  mutate(miRNA_family = str_replace_all(miRNA_family, '/', '='),
         miRNA_family = str_replace_all(miRNA_family, '\\.', '_')) %>%
  dplyr::filter(loop == 1) %>%
  group_by(lnc_trans, miRNA_family) %>% summarise(inter_num = n()) %>% 
  ungroup() %>%
  left_join(gene2trans, by = c('lnc_trans' = 'trans_id')) 

miRNA_lncRNA <- dplyr::filter(miRNA_inter_final, 
              miRNA_family %in% cemiRNA_family,
              inter_num >=3 ) %>%
  dplyr::select(miRNA_family, gene_id, inter_num) %>%
  mutate(miRNA_family = str_replace_all(miRNA_family, 'miR-181-5p', 'miR-181a-5p'),
         miRNA_family = str_replace_all(miRNA_family, 'miR-425-5p=489-3p', 'miR-425-5p'),
         miRNA_family = str_replace_all(miRNA_family, 'miR-27-3p', 'miR-27b-3p'),
         miRNA_family = str_replace_all(miRNA_family, 'miR-140-3p_2=497', 'miR-140-3p'),
         miRNA_family = str_replace_all(miRNA_family, 'miR-141-3p=200a-3p', 'miR-141-3p'),
         miRNA_family = str_replace_all(miRNA_family, 'miR-7-5p', 'miR-7a-5p'))

# filter the miRNA target gene proved by experiments 
miRTarBase <- read_excel("data/miRTarBase_MTI.xlsx")
miRNA_coding <- miRTarBase %>% dplyr::filter(miRNA %in% cemiRNA,
                                         `Support Type` == 'Functional MTI') 

miRNA_coding <- dplyr::select(miRNA_coding, miRNA, `Target Gene`) %>%
  mutate(miRNA = str_remove_all(miRNA, 'mmu-'),
         inter_num = 1) %>%
  distinct(miRNA, `Target Gene`, inter_num) %>%
  dplyr::select(miRNA_family = miRNA,
                gene_id = `Target Gene`,
                inter_num)

ceRNA_network <- rbind(miRNA_coding, miRNA_lncRNA) 
ceRNA_network_anno <- data.frame(gene_id = c(ceRNA_network$miRNA_family,
                                             ceRNA_network$gene_id)) %>%
  mutate(class = if_else(gene_id %in% lnclist, 'lncRNA', 
                         if_else(gene_id %in% cemRNAlist,
                                 'miRNA', 'coding_gene'))) %>%
  distinct(gene_id, class)

# save
write.csv(ceRNA_network, file = 'outcomes/ceRNA_network.csv', quote = F, row.names = F)
write.csv(ceRNA_network_anno, file = 'outcomes/miRNA_network_anno.csv',
          quote = F, row.names = F)







