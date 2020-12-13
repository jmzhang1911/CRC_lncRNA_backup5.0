source('Utils.R')


myenrichplot <- function(df, title, prob = FALSE){
  if(length(prob) > 0){
    ggdata <- df %>%
      mutate(enrich_type = factor(enrich_type,
                                  levels = c('GO-CC-enrichment',
                                             'GO-BP-enrichment',
                                             'KEGG-enrichment'))) %>%
      dplyr::filter(p.adjust < 0.05,
                    Description %in% prob) %>%
      group_by(enrich_type) %>% 
      arrange(enrich_type) 
  }
  else{
    ggdata <- df %>%
      mutate(enrich_type = factor(enrich_type, 
                                  levels = c('GO-CC-enrichment',
                                             'GO-BP-enrichment',
                                             'KEGG-enrichment'))) %>%
      dplyr::filter(p.adjust < 0.05) %>%
      group_by(enrich_type) %>% 
      slice_max(Count, n = 3) %>%
      arrange(enrich_type)
  }
  
  ggdata2 <- ggdata %>% 
    mutate(Description = factor(Description, levels = Description),
           pv = 0.05 - p.adjust)
  
  p <- ggplot(data = ggdata2,
              aes(x = log2(Count+1),
                  y = Description)) +
    geom_col(aes(fill = enrich_type,
                 alpha = p.adjust)) +
    scale_fill_nejm() +
    scale_alpha_continuous(range = c(0.4, 1)) +
    labs(x = 'log2(count+1)',
         y = '',
         title = title) +
    scale_x_continuous(expand = c(0, 0),
                       limits = c(0, 7)) +
    theme_bw()
  
  return(p)
}


myenrichgetdf <- function(genelist, title){
  filename <- basename(title)
  genelist = trygenes$gene_id
  filename = 'try'
  gene.df <- bitr(genelist, fromType = "ENSEMBL",
                  toType = c("SYMBOL", "ENTREZID"),
                  OrgDb = org.Mm.eg.db)
  go_cc <- enrichGO(gene       = gene.df$ENSEMBL,
                    OrgDb      = org.Mm.eg.db,
                    keyType    = 'ENSEMBL',
                    ont        = "CC",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.1,
                    qvalueCutoff = 0.05) %>%
    as.data.frame() %>%
    dplyr::tbl_df() %>%
    mutate(group = filename,
           enrich_type = 'GO-CC-enrichment')
  
  go_bp <- enrichGO(gene       = gene.df$ENSEMBL,
                    OrgDb      = org.Mm.eg.db,
                    keyType    = 'ENSEMBL',
                    ont        = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.1,
                    qvalueCutoff = 0.05) %>%
    as.data.frame() %>%
    dplyr::tbl_df() %>%
    mutate(group = filename,
           enrich_type = 'GO-BP-enrichment')
  go_mf <- enrichGO(gene       = gene.df$ENSEMBL,
                    OrgDb      = org.Mm.eg.db,
                    keyType    = 'ENSEMBL',
                    ont        = "MF",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.1,
                    qvalueCutoff = 0.05) %>%
    as.data.frame() %>%
    dplyr::tbl_df() %>%
    mutate(group = filename,
           enrich_type = 'GO-BP-enrichment')
  
  kk <- enrichKEGG(gene = gene.df$ENTREZID,
                   organism = 'mmu',
                   pvalueCutoff = 0.05) %>% 
    as.data.frame() %>%
    dplyr::tbl_df() %>%
    mutate(group = filename,
           enrich_type = 'KEGG-enrichment')
  tmp <- rbind(go_cc, go_bp)
  tmp2 <- rbind(kk, go_mf)
  tmp3 <- rbind(tmp, tmp2)
  return(tmp2)
}


load('outcomes/masigpro/dynamicdf.RData')

trygenes <- dynamicdf[dynamicdf$cluster == '2',]
myenrich(trygenes$gene_id,title = 'tmp')
tmpdf <- myenrichgetdf(trygenes$gene_id, title = 'tmp')
View(tmpdf)




go_all <- enrichGO(gene       = gene.df$ENSEMBL,
                  OrgDb      = org.Mm.eg.db,
                  keyType    = 'ENSEMBL',
                  ont        = "ALL",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.1,
                  qvalueCutoff = 0.05) 
kk <- enrichKEGG(gene = gene.df$ENTREZID,
                 organism = 'mmu',
                 pvalueCutoff = 0.05) %>% 
  as.data.frame() %>%
  dplyr::tbl_df() %>%
  mutate(enrich_method = 'KEGG') %>%
  dplyr::select(enrich_method,ID:Count)
View(kk)

go_res <- as.data.frame(go_all) %>% as.data.frame() %>% 
  tbl_df() %>% rename(enrich_method = ONTOLOGY)
View(go_res)

df_all <- union_all(kk, go_res)
View(df_all)

tmpprob <- sample(unique(df_all$ID),8)

getdf <- function(genelist, prob){
  genelist = trygenes$gene_id
  prob = tmpprob
  
  gene.df <- bitr(genelist, fromType = "ENSEMBL",
                  toType = c("SYMBOL", "ENTREZID"),
                  OrgDb = org.Mm.eg.db)
  
  kk_res <- enrichKEGG(gene = gene.df$ENTREZID,
                       organism = 'mmu',
                       pvalueCutoff = 0.05) %>%
    as.data.frame() %>%
    dplyr::tbl_df() %>%
    mutate(enrich_method = 'KEGG') %>%
    dplyr::select(enrich_method,ID:Count)
  
  go_res <- enrichGO(gene       = gene.df$ENSEMBL,
                     OrgDb      = org.Mm.eg.db,
                     keyType    = 'ENSEMBL',
                     ont        = "ALL",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.1,
                     qvalueCutoff = 0.05) %>%
    as.data.frame() %>% 
    tbl_df() %>% rename(enrich_method = ONTOLOGY)
  
  df_all <- union_all(kk, go_res)
  
  adf <- dplyr::filter(df_all, ID %in% prob) %>% dplyr::count(enrich_method) %>%
    mutate(get = if_else(n >=10, 0, 10-n))
  
  bdf <- dplyr::filter(df_all, ID %in% prob)
  cdf <- left_join(df_all, adf, by = 'enrich_method') %>% dplyr::select(-n) %>%
    dplyr::filter(ID %!in% prob) %>% group_by(enrich_method) %>% sample_n(get) %>%
    union_all(bdf) %>% 
    mutate(Description = factor(Description, levels = unique(Description)))
}

View(cdf)



ggdata2 <- df_all %>% 
  mutate(Description = factor(Description, levels = unique(Description))) %>%
  group_by(enrich_method) %>% slice_head(n = 10)

ggplot(data = cdf,
            aes(x = log2(Count+1),
                y = Description)) +
  geom_col(aes(fill = enrich_method,
               alpha = p.adjust)) +
  scale_fill_nejm() +
  scale_alpha_continuous(range = c(1, 0.4)) +
  labs(x = 'log2(count+1)',
       y = '') +
  scale_x_continuous(expand = c(0, 0),
                     limits = c(0, 10)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12, face = 'bold'),
        axis.text.y = element_text(size = 8, face = 'bold'))

View(ggdata2)






load('outcomes/inputdata/input.RData')
lnclist <- input_matrix_count$lnc_genelist


cemiRNA_family <- c('miR-181-5p', 'miR-128-3p', 'miR-425-5p=489-3p',
                    'miR-27-3p', 'miR-183-5p', 'miR-140-3p_2=497',
                    'miR-141-3p=200a-3p', 'miR-212-5p', 'miR-7-5p')



cemiRNA <- c('mmu-miR-181a-5p', 'mmu-miR-128-3p', 'mmu-miR-425-5p',
             'mmu-miR-27b-3p', 'mmu-miR-183-5p', 'mmu-miR-140-3p',
             'mmu-miR-141-3p', 'mmu-miR-212-5p', 'mmu-miR-7a-5p')

cemRNAlist <- str_remove_all(cemiRNA, 'mmu-')
miRTarBase <- read_excel("data/miRTarBase_MTI.xlsx")

ceRNA_df <- miRTarBase %>% dplyr::filter(miRNA %in% cemiRNA,
                                        `Support Type` == 'Functional MTI') 
miRNA_targetlist <- unique(ceRNA_df$`Target Gene`)
load('anno.RData')
ceRNA_target_df <- getBM(attributes =c('ensembl_gene_id',
                                    'external_gene_name',
                                    "description",
                                    "entrezgene_id"),
                      filters = 'external_gene_name',
                      values = miRNA_targetlist,
                      mart = mart)
df <- mygetenrichedgenefromkk2(ceRNA_target_df$ensembl_gene_id)
df2 <- myenrichgetdf(ceRNA_target_df$ensembl_gene_id,'miRNA_tar')
prob <- c('Colorectal cancer', 
          'Wnt signaling pathway',
          'JAK-STAT signaling pathway',
          'TGF-beta signaling pathway',
          'Notch signaling pathway',
          'p53 signaling pathway',
          'MAPK signaling pathway',
          'mTOR signaling pathway',
          'ErbB signaling pathway',
          'canonical Wnt signaling pathway',
          'positive regulation of canonical Wnt signaling pathway',
          'transcription regulator complex',
          'RNA polymerase II transcription regulator complex',
          'chromatoid body',
          'NF-kappa B signaling pathway',
          'EGFR tyrosine kinase inhibitor resistance',
          'positive regulation of DNA binding',
          'lymphocyte differentiation',
          'cytokine production involved in inflammatory response')

myenrichplot(df2, title = 'ce-network miRNA target genes', prob = prob)


myenrich(ceRNA_target_df$ensembl_gene_id, title = 'ceRNA_target_gene')


for(i in cemiRNA){
  tmplist <- unique(ceRNA_df[ceRNA_df$miRNA== i,]$`Target Gene`)
  tmpdf <- getBM(attributes =c('ensembl_gene_id',
                               'external_gene_name',
                               "description",
                               "entrezgene_id"),
                 filters = 'external_gene_name',
                 values = tmplist,
                 mart = mart)
  genelist <- tmpdf$ensembl_gene_id
  myenrich(genelist, title = i)
}




