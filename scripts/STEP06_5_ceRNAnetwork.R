source('Utils.R')
load('outcomes/ceRNAnetwork/key_miRNA_family.RData')
load('outcomes/ceRNAnetwork/lncRNA_miRNA_res/lncRNA_miRNA_res.RData')
load('outcomes/inputdata/input.RData')



#===> filter dynamic coding genes
load('outcomes/candidate/dy_results.RData')
dy_coding_genes <- dynamic_de_df[dynamic_de_df$gene_type == 'coding',]




#===> key miRNAs and lncRNAs
ceRNA_df <- dplyr::filter(miRNA_inter_final, miRNA_family %in% key_miRNA_family) %>% 
  left_join(input_matrix_count$exon_trans, by = c('lnc_trans' = 'trans_id')) %>%
  dplyr::select(-exon_num, -translen) %>%
  mutate(miRNA_family = str_remove_all(miRNA_family,'=.*$|_.*$'))

length(unique(ceRNA_df$miRNA_family)) # 36 key miRNA
length(unique(ceRNA_df$lnc_trans)) # 101 lncRNA trans
length(unique(ceRNA_df$gene_id)) # 101 lncRNA genes



#===> key miRNAs target genes
miRTarBase <- read_excel("data/db/miRTarBase_MTI.xlsx")
prob <- unique(ceRNA_df$miRNA_family)
prob[prob == 'miR-200bc-3p'] <- 'miR-200b-3p'
# section1
miRNA_tar1 <- miRTarBase %>% mutate(miRNA = str_extract(miRNA, 'miR.*')) %>% 
  dplyr::filter(miRNA %in% prob, 
                `Support Type` == 'Functional MTI',
                `Species (Target Gene)` == 'Mus musculus') %>%
  dplyr::select(miRNA, tar_gene = `Target Gene`, 
                species = `Species (miRNA)`,
                tar_gene_entrez = `Target Gene (Entrez ID)`)
a <- unique(miRNA_tar1$miRNA);length(a)
prob[prob %!in% a]


# section2
prob2 <- prob[prob %!in% a]
miRNA_tar2 <- miRTarBase %>% mutate(miRNA_tmp = str_replace_all(miRNA, '[a-d]',''),
                            miRNA_tmp = str_extract(miRNA_tmp, '(miR|let).*')) %>%
  dplyr::filter(miRNA_tmp %in% prob2, 
                `Support Type` == 'Functional MTI',
                `Species (Target Gene)` == 'Mus musculus') %>%
  dplyr::select(miRNA = miRNA_tmp, tar_gene = `Target Gene`, 
                species = `Species (miRNA)`,
                tar_gene_entrez = `Target Gene (Entrez ID)`)
a <- unique(miRNA_tar2$miRNA);length(a)
prob2[prob2 %!in% a]
miRNA_tar <- rbind(miRNA_tar1, miRNA_tar2)
table(miRNA_tar$miRNA)




#===> make dynamic ceRNA network
genelist <- unique(miRNA_tar$tar_gene_entrez)
gene.df <- bitr(genelist, fromType = "ENTREZID",
                toType = c("SYMBOL", "ENSEMBL"),
                OrgDb = org.Mm.eg.db)
myenrich(gene.df$ENSEMBL,'outcomes/ceRNAnetwork/all_miRNA_targenes')

ceRNA_net <- miRNA_tar %>% mutate(tar_gene_entrez = as.character(tar_gene_entrez)) %>% 
  left_join(gene.df, by = c('tar_gene_entrez' = 'ENTREZID')) %>% 
  left_join(ceRNA_df, by = c('miRNA' = 'miRNA_family')) %>%
  rename(tar_lnc_trans = lnc_trans,
         tar_lnc_gene = gene_id) %>%
  dplyr::filter(ENSEMBL %in% dy_coding_genes$gene_id) %>%
  mutate(try = str_c(ENSEMBL, '=', tar_lnc_gene))
myenrich(unique(ceRNA_net$ENSEMBL),'outcomes/ceRNAnetwork/all_miRNA_targenes_dynamic')


length(unique(ceRNA_net$tar_gene_entrez))
length(unique(ceRNA_net$miRNA))
length(unique(ceRNA_net$tar_lnc_gene))




#===> give it a shot
load('outcomes/candidate/crc_corgenelist.RData')
corlnk_gene_lnc = crc_corgenelist$cor_lnc_crcgene %>% mutate(lnk = str_c(source,'=',target))

dplyr::filter(ceRNA_net, try %in% corlnk_gene_lnc$lnk) %>% distinct()




#===> give it to cytoscape
network1 <- ceRNA_net %>% dplyr::select(source = miRNA, 
                                        target = tar_lnc_gene, 
                                        inter = inter_num) %>% distinct() %>% drop_na()
table(network1$source)
network2 <- ceRNA_net %>% dplyr::select(source = miRNA, 
                                        target = SYMBOL) %>% 
  mutate(inter = 1) %>% distinct() 
network <- rbind(network1,network2)
network_anno <- data.frame(Id = unique(c(network$source, network$target))) %>%
  mutate(class = if_else(Id %in% input_matrix_count$lnc_genelist, 'lnc-gene',
                         if_else(str_detect(Id, '^(miR|let)-'), 'miRNA', 'coding-gene')),
         Label = str_remove_all(Id, 'NONMMUG|ENSMUSG'))
write.csv(network, file = 'outcomes/ceRNAnetwork/network.csv',quote = F, row.names = F)
write.csv(network_anno, file = 'outcomes/ceRNAnetwork/network_anno.csv',quote = F, row.names = F)



















#=========================

ceRNA_net <- miRNA_tar %>% 
  mutate(`Target Gene (Entrez ID)` = as.character(`Target Gene (Entrez ID)`)) %>%
  left_join(gene.df, by = c(`Target Gene (Entrez ID)`='ENTREZID')) %>%
  dplyr::select(miRNA, tar_genes_symbol = SYMBOL, tar_gene_ensembl = ENSEMBL) %>%
  mutate(miRNA = str_replace_all(miRNA, 'miR-181a-5p','miR-181-5p'),
         miRNA = str_replace_all(miRNA, 'miR-190a-5p','miR-190-5p'),
         miRNA = str_replace_all(miRNA, 'miR-15a-5p','miR-15-5p'),
         miRNA = str_replace_all(miRNA, 'miR-193a-3p','miR-193-3p'),
         miRNA = str_replace_all(miRNA, 'miR-125a-5p','miR-125-5p'),
         miRNA = str_replace_all(miRNA, 'let-7b-5p','let-7-5p')) %>%
  left_join(ceRNA_df, by = c('miRNA'='miRNA_family')) %>% distinct()
  
# give it a shot
ceRNA_net %>% mutate(lnk = str_c(tar_gene_ensembl,'=',gene_id)) %>%
  dplyr::filter(lnk %in% corlnk_gene_lnc$lnk)

ceRNA_net %>% dplyr::select(miRNA, gene_id) %>%
  distinct() %>%
  group_by(miRNA) %>% summarise(count = n())


#===> give it to cytoscape
network1 <- ceRNA_net %>% dplyr::select(source = miRNA, 
                                        target = gene_id, 
                                        inter = inter_num) %>% distinct()
network2 <- ceRNA_net %>% dplyr::select(source = miRNA, 
                                        target = tar_genes_symbol) %>% 
  mutate(inter = 1) %>% distinct() %>% drop_na()
network <- rbind(network1,network2)
network_anno <- data.frame(source = unique(c(network$source, network$target))) %>%
  mutate(class = if_else(source %in% input_matrix_count$lnc_genelist, 'lnc-gene',
                         if_else(str_detect(source, '^miR-'), 'miRNA', 'coding-gene')))
write.csv(network, file = 'outcomes/ceRNAnetwork/network.csv',quote = F)
write.csv(network_anno, file = 'outcomes/ceRNAnetwork/network_anno.csv',quote = F)







