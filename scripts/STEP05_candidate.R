rm(list = ls())
source('Utils.R')




#===> our dynamic-DE-expression genes = dynamic + DE
#==> 4875 = 807lnc + 4068coding

load('outcomes/inputdata/input.RData')
load('outcomes/DEanalysis/DE_results_inte.RData') # DE genes
load('outcomes/masigpro/dynamicdf.RData') # dynamic genes
mkdir('outcomes/candidate')

dedf <- dplyr::filter(DE_results_inte$DE_each_inte, 
                      direction %in% c('up','down'),
                      abs(log2FoldChange) >= 1) 
table(unique(dynamicdf$gene_id) %in% unique(dedf$gene_id)) 

#4875 dynamic and de genes
dynamic_de_df <- dynamicdf %>% dplyr::filter(gene_id %in% dedf$gene_id) 
#dynamic and DE coding:4068;noncoding:807
table(dynamic_de_df$gene_type)
#get the dy genes
dy_lncdf <- dynamic_de_df[dynamic_de_df$gene_type == 'lnc',];dim(dy_lncdf) #807


# dynamic or not dynamic lnclist
if(T){
  lncRNA_expr_class_df <- input_matrix_count$lncRNA_count %>% dplyr::select(gene_id) %>%
    mutate(expr_class = if_else(gene_id %in% dy_lncdf$gene_id, 'dynamic','not_dynamic'))
  dynamic_lnc <- dy_lncdf$gene_id
  not_dynamic_lnc <- lncRNA_expr_class_df %>% 
    dplyr::filter(expr_class == 'not_dynamic') %>%
    pull(gene_id)
  lncRNA_expr_class_list <- list(lncRNA_expr_class_df = lncRNA_expr_class_df,
                                 dynamic_lnc =dynamic_lnc,
                                 not_dynamic_lnc = not_dynamic_lnc)
  
  write.table(dynamic_lnc, file = 'outcomes/candidate/dy_lnclist.txt',
              quote = F, row.names = F, col.names = F)
  save(lncRNA_expr_class_list, file = 'outcomes/candidate/lncRNA_expr_class_list.RData')
}

# dynamic coding genes enrichment analysis
source('scripts/STEP09_endplot.R')
df <- mygetallpathwaydf(genelist = dynamic_de_df$gene_id, 
                        prob = CRC_related_pathway_prob)
mypathwayplot(df, 'Enrchment of dynamic coding genes ')

View(df)

# CRC dynamic-coding genes list
if(T){
  dy_crcgene <- mygetenrichedgenefromkk2(dynamic_de_df$gene_id) %>% 
    distinct(ENSEMBL) %>% pull(ENSEMBL);length(dy_crcgene) #105
  save(dynamic_de_df, dy_crcgene, file = 'outcomes/candidate/dy_results.RData')
}




#===> calculate the correlation between dy-coding and dy-noncoding
load('outcomes/inputdata/input.RData')
load('outcomes/candidate/dy_results.RData')
load('outcomes/candidate/lncRNA_expr_class_list.RData')

library(Hmisc)
dim(dynamic_de_df)

if(T){
  cor_results <- input_matrix_count$mRNA_lncRNA_count %>%
    dplyr::filter(gene_id %in% dynamic_de_df$gene_id) %>% 
    column_to_rownames(var = 'gene_id') %>% myddnor() %>%
    as.matrix() %>% t() %>%
    mycoranalysis2()
  #save(cor_results, file = 'outcomes/candidate/cor_results.RData') 
}




#===> filter dy-crc-coding relatived dy-lnc
if(T){
  cor_lnc_crcgene <- dplyr::filter(cor_results, source %in% dy_crcgene,
                                   target %in% lncRNA_expr_class_list$dynamic_lnc,
                                   abs(r_value) >=0.7, p_value < 0.05) 
  
  cor_crcgene <- unique(cor_lnc_crcgene$source);length(cor_crcgene) #104
  cor_lnc <- unique(cor_lnc_crcgene$target);length(cor_lnc) #737
  crc_corgenelist <- list(cor_lnc_crcgene = cor_lnc_crcgene, 
                          cor_lnc = cor_lnc, 
                          cor_crcgene = cor_crcgene)
  
  save(crc_corgenelist, file = 'outcomes/candidate/crc_corgenelist.RData')
}





#===> Find the max trans of dy-lncRNA
#===> Find the TSS of dy-lncRNA
###########=====> do it in jobs >>>>>>###########
#===> saved as FindTSS_maxlen.sh in jobs
###########=====> do it in jobs >>>>>>###########
load('outcomes/candidate/dy_results.RData')
maxtrans_df <- read.table('outcomes/candidate/maxtrans.txt', sep = '\t', 
                          comment.char = '#',header = T)
save(maxtrans_df, file = 'outcomes/candidate/maxtrans_df.RData')
write.table(maxtrans_df$trans_id, file = 'outcomes/candidate/maxlentransofdylnc.txt',
            quote = F, row.names = F, col.names = F)











#===> find candidate lnc genes
load('outcomes/coranalysis/cor_analysis.RData')
de_lnc <- dplyr::filter(DE_results_inte$DE_each_inte, 
                      direction %in% c('up','down'),
                      abs(log2FoldChange) >= 1,
                      gene_id %in% input_matrix_count$lnc_genelist) 

candidate_lnc_df <- cor_results %>% 
  dplyr::filter(source %in% dy_de_crc_coding, target %in% dy_de_lnc,
                abs(r_value) >= 0.7, p_value < 0.05,
                target %in% de_lnc$gene_id) 
length(unique(candidate_lnc_df$source)) # 54 dy and de genes
length(unique(candidate_lnc_df$target)) # 659 dy and de lnc FC>=2

# co-expressied with crc-dy-gene lncRNA 
# 659 dy_de_FC>=2 lnc co-expression with 54 dy_de_FC>=2 gene
# candidate_lnc = candidate_lnctrans_help.txt + candidate_lnctrans_help2.txt
co_expr_candidate_lnc <- dplyr::filter(input_matrix_count$exon_trans,
                         gene_id %in% candidate_lnc_df$target) %>%
  group_by(gene_id) %>% top_n(1, translen) %>%
  group_by(gene_id) %>% top_n(1, exon_num) %>%
  group_by(gene_id) %>% top_n(1, trans_id) %>%
  pull(trans_id)

write.table(co_expr_candidate_lnc, file = 'outcomes/candidate_lnc.txt', 
            col.names = F, row.names = F, quote = F)


# all dy-lnc analysis by annolnc2
# old file = candidate_lnctrans_help3.txt
all_candidate_lnc <- dplyr::filter(input_matrix_count$exon_trans,
                               gene_id %in% dy_de_lnc,
                               gene_id %!in% candidate_lnc_df$target) %>%
  group_by(gene_id) %>% top_n(1, translen) %>%
  group_by(gene_id) %>% top_n(1, exon_num) %>%
  group_by(gene_id) %>% top_n(1, trans_id) %>%
  pull(trans_id)
write.table(all_candidate_lnc, file = 'outcomes/all_candidate_lnc.txt', 
            col.names = F, row.names = F, quote = F)




#===> lncRNA classes
lncRNA_classes <- read_table2("data/lncRNA_classes_pbs.txt")
lncRNA_classed_df <- lncRNA_classes %>% 
  dplyr::filter(isBest == 1) %>% 
  mutate(class = if_else(direction == 'sense' & subtype == 'overlapping',
                         'overlapping_sence', 'NA'),
         class = if_else())

table(lncRNA_classes$direction)
table(lncRNA_classes$type)
table(lncRNA_classes$subtype)
table(lncRNA_classes$location)


table4 <- lncRNA_classes %>% dplyr::filter(isBest == 1, direction == 'antisense',
                                           type == 'genic') %>% 
  mutate(lnk = str_c(partnerRNA_gene, lncRNA_gene, sep = '='))


cor_results %>% head()
de_lnc <- DE_results_inte$DE_each_inte %>% 
  dplyr::filter(direction == c('up','down'),
                abs(log2FoldChange) >= 4) %>%
  distinct(gene_id) %>%
  pull(gene_id)
  

colnames(input_matrix_count$mRNA_lncRNA_count)
df1 <- input_matrix_count$mRNA_lncRNA_count %>% dplyr::filter(gene_id %in% de_lnc) %>%
  column_to_rownames(var = 'gene_id') %>% myddnor() %>% as.data.frame() %>%
  rownames_to_column(var = 'gene_id') %>%
  pivot_longer(cols = control_1:week10_3, names_to = 'time', values_to = 'levels') %>%
  mutate(time = str_replace_all(time, '_[0-9]', ''),
         time = factor(time, levels = c('control','week2','week4','week7','week10')))

p1 <- ggplot(df1, aes(x = time, y = log10(levels+1))) +
  geom_point() +
  geom_line(aes(group = gene_id), size = 1.25) 

df2 <- input_matrix_count$mRNA_lncRNA_count %>% dplyr::filter(gene_id %in% dy_de_lnc) %>%
  column_to_rownames(var = 'gene_id') %>% myddnor() %>% as.data.frame() %>%
  rownames_to_column(var = 'gene_id') %>%
  pivot_longer(cols = control_1:week10_3, names_to = 'time', values_to = 'levels') %>%
  mutate(time = str_replace_all(time, '_[0-9]', ''),
         time = factor(time, levels = c('control','week2','week4','week7','week10')))

p2 <- ggplot(df2, aes(x = time, y = log10(levels+1))) +
  geom_point() +
  geom_line(aes(group = gene_id), size = 1.25) 

p1 + p2

table(de_lnc %in% dy_de_lnc)
table(dy_de_lnc %in% de_lnc)

table3 <- cor_results %>% 
  dplyr::filter(source %in% dy_de_crc_coding, target %in% de_lnc,
                                        abs(r_value) >=0.7, p_value < 0.05) %>%
  mutate(lnk = str_c(source, target, sep = '=')) %>%
  dplyr::filter(lnk %in% table4$lnk) 

table5 <- cor_results %>% 
  dplyr::filter(source %in% dy_de_crc_coding, target %in% dy_de_lnc,
                                        abs(r_value) >= 0.7, p_value < 0.05,
                                        target %in% de_lnc) 
length(unique(table5$target))  
length(unique(table5$source))
load('outcomes/inputdata/input.RData')
candidata <- unique(table5$target)
count_len <- input_matrix_count$exon_trans
head(count_len)
dim(count_len)

key_lnc <- dplyr::filter(count_len,
                         gene_id %in% candidata) %>%
  group_by(gene_id) %>% top_n(1, translen) %>%
  group_by(gene_id) %>% top_n(1, exon_num) %>%
  group_by(gene_id) %>% top_n(1, trans_id) %>%
  pull(trans_id)

write.table(key_lnc, file = 'outcomes/candidate_lnctrans.txt', 
            col.names = F, row.names = F,
            quote = F)




tmp_lnc <- table5 %>% distinct(target) %>% pull(target)  
tmp_coding <- lncRNA_classes %>%
  dplyr::filter(isBest == 1, lncRNA_gene %in% tmp_lnc) %>% 
  distinct(partnerRNA_gene) %>% pull(partnerRNA_gene)
tmp_enrich <- mygetenrichedgenefromkk2(tmp_coding)

df <- left_join(table3, table4, by = 'lnk') %>%
  left_join(table, by = c('source' = 'ENSEMBL')) 


