source('Utils.R')


load('outcomes/inputdata/input.RData')
load('outcomes/inputdata/anno.RData');mart
trans2gene <- input_matrix_count$exon_trans
lncbed <- read.table('data/final_lnc.bed')
filepath <-list.files('Annolnc2_results/', full.names = T)
mkdir('outcomes/lnc_proteins')




#===> interaction with TF
if(T){
  TFpath <- paste(filepath, '/', basename(filepath), 
                  '_trans1/transcription_regulation/filter_details.txt', sep = '')
  TF_list <- sapply(filepath, function(x){
    TFpath <- paste(x, '/', basename(x), 
                    '_trans1/transcription_regulation/filter_details.txt', sep = '')
    if(file.info(TFpath)$size > 0 & !is.na(file.info(TFpath)$size)){
      df <- read.table(TFpath, header = T, sep = '\t')
      df$trans_id <- basename(x)
      return(df)}}, USE.NAMES = T)
  
  for(i in names(TF_list)){
    if(is.null(TF_list[[i]])){TF_list[i] <- NULL}
  }
  
  # do.call can bind the list With or without empty elements in list
  TF_df <- do.call(rbind, TF_list);rownames(TF_df) <- NULL
}

# filter the colon related TF
TF_df_res <- TF_df %>% dplyr::filter(str_detect(Cell_type, 'colon|intestin')) %>%
  dplyr::filter(up_TSS == 'Yes' | overlap_TSS == 'Yes')
table(TF_df_res$TF) # 21 TFs

# annotation
TF_anno <- getBM(attributes =c('ensembl_gene_id',
                               'external_gene_name',
                               "description",
                               "entrezgene_id",
                               'chromosome_name',
                               'start_position',
                               'end_position'),
                 filters = 'external_gene_name',
                 values = unique(TF_df_res$TF),
                 mart = mart) 

TF_circol_df <- left_join(TF_df_res, trans2gene, by = 'trans_id') %>%
  dplyr::select(TF, gene_id) %>% distinct() %>%
  left_join(TF_anno, c('TF' = 'external_gene_name')) %>% 
  left_join(lncbed, by = c('gene_id' = 'V4')) %>%
  dplyr::rename(lncchr = V1,
                lncstart = V2,
                lncend = V3,
                chr = chromosome_name, 
                start = start_position, 
                end = end_position, 
                name = TF)

TF_bed <- TF_circol_df %>% dplyr::select(chr, start, end, name) %>%
  mutate(class = 'TF', chr = str_c('chr', chr))
TF_lnc_bed <- TF_circol_df %>% dplyr::select(chr = lncchr, start = lncstart,
                                             end = lncend, name = gene_id) %>%
  mutate(class = 'lnc', chr = str_c('chr', chr))

TF_prob <- dplyr::union_all(TF_bed, TF_lnc_bed) %>% dplyr::count(name) %>%
  dplyr::filter(n > 4) %>% pull(name)
TF_anno_bed <- dplyr::union(TF_bed, TF_lnc_bed) %>% dplyr::filter(name %in% TF_prob) %>%
  mutate(name = str_remove_all(name, 'ENSMUSG|NONMMUG'))

# visualization 
pdf('outcomes/lnc_proteins/lnc_TFS_circols.pdf',width = 8, height = 8)
library(circlize)
circos.initializeWithIdeogram(species = 'mm10')
circos.genomicLabels(TF_anno_bed, labels.column = 4, side = "outside", 
                     col = as.numeric(factor(TF_anno_bed[[5]])), 
                     line_col = as.numeric(factor(TF_anno_bed[[5]])), 
                     labels_height = 0.1,
                     cex = 0.35, connection_height = 0.08)
circos.genomicLink(TF_bed, TF_lnc_bed, col = rand_color(nrow(TF_lnc_bed)), 
                   border = NA)
text(0, 0, "Interactions between dy-lncRNAs and TFs", cex = 0.8)
dev.off()
circos.clear()

s#===> interaction with protein
if(T){
  Pro_list <- sapply(filepath, function(x){
    Propath <- paste(x, '/', basename(x), 
                     '_trans1/protein_interaction/intersected.clip.peak', sep = '')
    if(file.info(Propath)$size > 0 & !is.na(file.info(Propath)$size)){
      df <- read.table(Propath, header = F, sep = '\t')
      df$trans_id <- basename(x)
      return(df)}}
    , USE.NAMES = T)
  for(i in names(Pro_list)){
    if(is.null(Pro_list[[i]])){Pro_list[i] <- NULL}}
  
  Pro_df <- do.call(rbind, Pro_list);rownames(Pro_df) <- NULL 
}

Pro_df <- Pro_df %>% mutate(V7 = tolower(V7),
                            V7 = str_replace_all(V7, 'roquin-1','rc3h1'),
                            V7 = str_replace_all(V7, 'aid','aicda'),
                            V7 = str_replace_all(V7, '^celf$','cebpd'),
                            V7 = str_replace_all(V7, 'dazl\\?','dazl'),
                            V7 = str_replace_all(V7, 'rod1/ptbp3','ptbp3'),
                            V7 = str_replace_all(V7, 'qki5','Qk'),
                            V7 = str_replace_all(V7, '^nova$','Hnrnpk'))

Pro_anno <- getBM(attributes =c('ensembl_gene_id',
                                'external_gene_name',
                                'description',
                                'entrezgene_id',
                                'chromosome_name',
                                'start_position',
                                'end_position'),
                  filters = 'external_gene_name',
                  values = unique(Pro_df$V7),
                  mart = mart) %>%
  mutate(external_gene_name = tolower(external_gene_name)) %>%
  dplyr::filter(str_detect(chromosome_name, '^[0-9]+|^[XY]'))

# proteins are not annoed
length(Pro_anno$external_gene_name) #58
length(unique(tolower(Pro_df$V7))) #58

Pro_circol_df <- Pro_df %>% mutate(V7 = tolower(V7)) %>% distinct(V7, trans_id) %>% 
  left_join(Pro_anno, by = c('V7'='external_gene_name')) %>%
  dplyr::select(chr = chromosome_name,
                start = start_position,
                end = end_position,
                name = V7,
                symbol = ensembl_gene_id,
                trans_id = trans_id) %>%
  left_join(trans2gene, by = c('trans_id'= 'trans_id')) %>%
  dplyr::select(-exon_num, -translen) %>%
  left_join(lncbed, by = c('gene_id' = 'V4')) %>%
  dplyr::rename(lncchr = V1,
                lncstart = V2,
                lncend = V3) %>%
  dplyr::select(-V5, -V6) %>% drop_na() %>%
  dplyr::filter(!grepl('CHR', chr))

Pro_bed <- Pro_circol_df %>% dplyr::select(chr, start, end, name) %>%
  mutate(chr = str_c('chr', chr),
         class = 'pro')
Pro_lnc_bed <- Pro_circol_df %>% dplyr::select(chr = lncchr,
                                               start = lncstart,
                                               end = lncend,
                                               name = gene_id) %>%
  mutate(chr = str_c('chr', chr),
         class = 'lnc') 
Pro_prob <- dplyr::union_all(Pro_bed, Pro_lnc_bed) %>% dplyr::count(name) %>%
  dplyr::filter(n > 5) %>% pull(name)
Pro_anno_bed <- dplyr::union(Pro_bed, Pro_lnc_bed) %>% dplyr::filter(name %in% Pro_prob) %>%
  mutate(name = str_remove_all(name, 'ENSMUSG|NONMMUG'))
table(Pro_anno_bed$class)

# visualization
pdf('outcomes/lnc_proteins/lnc_Pro_circols.pdf',width = 8, height = 8)
circos.initializeWithIdeogram(species = 'mm10')
circos.genomicLabels(Pro_anno_bed, labels.column = 4, side = "outside", 
                     col = as.numeric(factor(Pro_anno_bed[[5]])), 
                     line_col = as.numeric(factor(Pro_anno_bed[[5]])), 
                     labels_height = 0.1,
                     cex = 0.35, connection_height = 0.08)

circos.genomicLink(Pro_bed, Pro_lnc_bed, 
                   col = rand_color(nrow(Pro_lnc_bed)),
                   border = NA)
text(0, 0, "Interactions between dy-lncRNAs and Proteins", cex = 0.8)
dev.off()
circos.clear()




