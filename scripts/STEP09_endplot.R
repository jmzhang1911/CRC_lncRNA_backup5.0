source('Utils.R')


#====> plot result of pathway 
mygetallpathwaydf <- function(genelist, prob=FALSE){

  prob = CRC_related_pathway_prob
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
  
  df_all <- union_all(kk_res, go_res) %>% dplyr::filter(p.adjust <= 0.05)
  
  if(!is.logical(prob)){
    # how much to get
    adf <- dplyr::filter(df_all, ID %in% prob) %>% dplyr::count(enrich_method) %>%
      mutate(get = if_else(n >=10, 0, 10-n))
    
    bdf <- dplyr::filter(df_all, ID %in% prob)
    tmpdf <- dplyr::filter(df_all, ID %!in% prob) %>% dplyr::count(enrich_method)

    cdf <- left_join(df_all, adf, by = 'enrich_method') %>% dplyr::select(-n) %>%
      dplyr::filter(ID %!in% prob) %>%
      mutate(get = if_else(is.na(get), 10, get)) %>%
      left_join(tmpdf, by = 'enrich_method') %>%
      mutate(get = if_else(get >= n, as.numeric(n), as.numeric(get))) %>%
      group_by(enrich_method) %>%
      dplyr::sample_n(get) %>% 
      union_all(bdf) %>% 
      mutate(Description = factor(Description, levels = unique(Description)))
    return(cdf)
  }
  return(df_all)
}


mypathwayplot <- function(data,title){
  p <- ggplot(data = data,
              aes(x = log10(Count+1),
                  y = Description)) +
    geom_col(aes(fill = enrich_method,
                 alpha = p.adjust)) +
    scale_fill_nejm() +
    scale_alpha_continuous(range = c(1, 0.4)) +
    labs(x = 'log10(count+1)',
         y = '',
         title = title) +
    scale_x_continuous(expand = c(0, 0),
                       limits = c(0, 3)) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 12, face = 'bold'),
          axis.text.y = element_text(size = 8, face = 'bold'))
  return(p)
}


CRC_related_pathway_prob <- c('mmu04520','mmu03015','mmu05210',
                              'mmu04310','mmu04630','mmu04350',
                              'mmu04330','mmu04115','mmu04010',
                              'mmu04150','mmu04012',
                              'GO:0016570','GO:0006338','GO:0032259',
                              'GO:0050821','GO:0042393','GO:0002039',
                              'GO:0035102')
