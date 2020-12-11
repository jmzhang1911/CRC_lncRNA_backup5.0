#===>
options(stringsAsFactors = F)
options(repos="http://mirrors.tuna.tsinghua.edu.cn/CRAN/")
options(BioC_mirror="http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")
options(warn =-1)



#===> my var 
`%!in%` = Negate(`%in%`)

#===> base
mkdir <- function(x){if(!dir.exists(x)){
  dir.create(x,recursive = T)}else(cat('existed dir\n'))}



#===> my functions
superlib <- function(x){
  if(x == 'super'){
    tmp <- c('tidyverse','PCAtools','DESeq2','pheatmap',
            'edgeR','biomaRt','ggsignif','cowplot','ggsci',
            'scales','patchwork','extrafont','Cairo','ggrepel',
            'cowplot','clusterProfiler','org.Mm.eg.db','readxl',
            'MeSH.db','Biostrings','seqinr','RTCGA','RTCGA.miRNASeq',
            'RTCGA.clinical','circlize','ggplotify','Hmisc')
    newtmp <- tmp[tmp %!in% installed.packages()[,'Package']]
    if(length(newtmp) > 0){
      lapply(newtmp, function(x){
        BiocManager::install(x,ask = F,update = F)})}
    lapply(tmp, function(x){suppressMessages(library(x, character.only=TRUE))})
}}
superlib(x='super')



#===> normalization
mycountFPKM <- function(expr_df, genelen_df){
  # both dfs must be matrix
  tmp1 <- t(apply(expr_df, 1, function(x){
    x / (apply(expr_df, 2, sum)/1000000)}))
  gene_lendf <- gene_lendf[match(rownames(tmp1), rownames(genelen_df)),]
  names(gene_lendf) <- rownames(tmp1)
  tmp2 <- apply(t(tmp1), 1, function(x){
    x / (gene_lendf / 1000)})
  return(tmp2)
}

mycountTPM <- function(expr_df, genelen_df){
  # both dfs must be matrix
  lendf <- gene_lendf[match(rownames(expr_df), rownames(genelen_df)),]
  names(lendf) <- rownames(expr_df)
  tmp <- apply(expr_df, 2, function(x){x / (lendf / 1000)})
  tmp2 <- t(apply(tmp, 1, function(x){ x / (apply(tmp, 2, sum) / 1000000)}))
  return(tmp2)
}



#===> DE analysis
mydeanalysis <- function(sample_data){
  results <- rownames_to_column(as.data.frame(sample_data), 
                                var = "gene_id") %>%
    dplyr::select(gene_id, log2FoldChange, pvalue, padj) %>%
    mutate(FC = 2 ** abs(log2FoldChange),
           direction = if_else(padj > 0.05, "ns", if_else(
             abs(log2FoldChange) < 1, "ns", if_else(
               log2FoldChange > 1, "up", "down"))))
  return(results)}



#===> DE_analysis_inte
mydeanalysis_inte <- function(DE_res_list, count_matrix_tb){
  
  # 提供DESeq2差异分析的results结果列表及counts矩阵
  # 1.整合所有vs差异表达结果
  # 2.整合所有vs差异表达结果+counts
  # 3.每个vs的具体结果
  # 4.每个vs的具体结果纵向合并=> 非常推荐
  
  for(i in names(DE_res_list)){
    a = mydeanalysis(DE_res_list[[i]])
    a$vsxx = i
    DE_res_list[[i]] = a
  }
  
  DE_res <- DE_res_list[[1]]
  name_tmp <- colnames(DE_res)[-1]
  
  for(i in c(2:length(DE_res_list))){
    DE_res <- full_join(DE_res, DE_res_list[[i]], by = 'gene_id')
  } 
  
  colnames(DE_res)[-1] <- paste(rep(name_tmp,length(DE_res_list)),
                                rep(names(DE_res_list), each = length(name_tmp)),
                                sep = '_') 
  DE_res_and_count <- full_join(DE_res, count_matrix_tb, by = "gene_id")
  DE_each_inte <- do.call(rbind, DE_res_list)
  rownames(DE_each_inte) <- NULL
  
  DE_result_pro <- list(
    'DE_res' = DE_res,                          #整合的差异表达结果
    'DE_res_and_count' = DE_res_and_count,      #+counts数
    'DE_each' = DE_res_list,
    'DE_each_inte' = DE_each_inte)              #每个差异结果
  
  
  return(DE_result_pro)
}



#===> cor analysis
mycoranalysis2 <- function(data){
  tmp_re <- rcorr(as.matrix(data))
  cor <- tmp_re$r
  p <- tmp_re$P
  ut <- upper.tri(cor) 
  cor_result <- data.frame(source = rownames(cor)[row(cor)[ut]],
                           target = rownames(cor)[col(cor)[ut]],
                           r_value = cor[ut],
                           p_value = p[ut])
  return(cor_result)} 

mycoranalysis <- function(data){
  tmp_res <- rcorr(as.matrix(tmp))
  result_R <- rownames_to_column(as.data.frame(tmp_res$r),
                                 var = "source") %>%
    gather(key = target, value = R_value, -1)
  result_P <- rownames_to_column(as.data.frame(tmp_res$P), 
                                 var = "source") %>%
    gather(key = target, value = P_value, -1)
  tmp_result <- full_join(result_R, result_P, 
                          by = c("source","target")) %>%
    dplyr::filter(P_value < 0.05, R_value > 0.8) %>%
    dplyr::filter(source != target)
  return(tmp_result)
}



#===> enrichment analysis
myenrich <- function(clu_gene,title){
  genelist <- clu_gene
  try({
    try(dev.off())
    gene.df <- bitr(genelist, fromType = "ENSEMBL",
                    toType = c("SYMBOL", "ENTREZID"),
                    OrgDb = org.Mm.eg.db)
    go_cc <- enrichGO(gene       = gene.df$ENSEMBL,
                      OrgDb      = org.Mm.eg.db,
                      keyType    = 'ENSEMBL',
                      ont        = "CC",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.1,
                      qvalueCutoff = 0.05)
    go_bp <- enrichGO(gene       = gene.df$ENSEMBL,
                      OrgDb      = org.Mm.eg.db,
                      keyType    = 'ENSEMBL',
                      ont        = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.1,
                      qvalueCutoff = 0.05)
    kk<-enrichKEGG(gene = gene.df$ENTREZID,
                   organism = 'mmu',
                   pvalueCutoff = 0.05);kk[1:30]
    
    print(paste(title, " ===>", sep = ''))
    filename <- basename(title)
    a<- barplot(go_bp, showCategory = 70,
                title = paste("The GO_BP enrichment analysis of",
                              filename))
    b <- barplot(go_cc, showCategory = 70,
                 title = paste("The GO_CC enrichment analysis of",
                               filename))
    c <- barplot(kk, showCategory = 70, 
                 title = paste("The KEGG enrichment analysis of", 
                               filename))
    
    pdf(paste(title, "_enrichment.pdf", sep = ''),
        width = 10, height = 32)
    print(a / b / c)
    dev.off() 
  })
}

myenrichgetdf <- function(genelist, title){
  filename <- basename(title)
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
  
  kk <- enrichKEGG(gene = gene.df$ENTREZID,
                   organism = 'mmu',
                   pvalueCutoff = 0.05) %>% 
    as.data.frame() %>%
    dplyr::tbl_df() %>%
    mutate(group = filename,
           enrich_type = 'KEGG-enrichment')
  tmp <- rbind(go_cc, go_bp)
  tmp2 <- rbind(tmp, kk)
  return(tmp2)
}


mygetenrichedgenefromkk2 <- function(genelist){
  
  # Colorectal cancer => mmu05210
  # Wnt signaling pathway => mmu04310
  # JAK-STAT signaling pathway => mmu04630
  # TGF-beta signaling pathway => mmu04350
  # Notch signaling pathway => mmu04330
  # p53 signaling pathway => mmu04115
  # MAPK signaling pathway => mmu04010
  # mTOR signaling pathway => mmu04150
  # ErbB signaling pathway => mmu04012
  
  pathwaylist <- c('mmu05210','mmu04310','mmu04630',
                   'mmu04350','mmu04330','mmu04115',
                   'mmu04010','mmu04150','mmu04012')
  gene.df <- bitr(genelist, fromType = "ENSEMBL",
                  toType = c("SYMBOL", "ENTREZID"),
                  OrgDb = org.Mm.eg.db)
  kk <- enrichKEGG(gene = gene.df$ENTREZID,
                   organism = 'mmu',
                   pvalueCutoff = 0.05)
  kk_res <- as.data.frame(kk)
  df <- dplyr::tbl_df(kk_res) %>%
    dplyr::filter(ID %in% pathwaylist) %>%
    separate_rows(geneID, sep = '/') %>%
    dplyr::select(ID, Description, geneID) %>%
    left_join(gene.df, by = c('geneID'='ENTREZID')) 
  return(df)
} #ENSEMBL SYMBOL ENTREZID


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



#===> DESeq2 normalization
myddnor <- function(df){
  condition <- factor(c(rep("control", 3), rep("week2", 3),
                        rep("week4", 3),rep("week7", 3),
                        rep("week10", 3)), 
                      levels = c("control", "week2",
                                 "week4","week7","week10"))
  colData <- data.frame(row.names = colnames(df), condition)
  dds <- DESeqDataSetFromMatrix(df, colData, design= ~ condition)
  dds <- DESeq(dds)
  nor_df <- counts(dds, normalized=TRUE)
  return(nor_df)
}



#===> volcano
myggvol <- function(df, title){
  ggplot(df, aes(x = df[,2], y = -log10(df[,3]))) +
    geom_point(size = 1, 
               aes(color = df[,4]), 
               show.legend = T) +
    scale_color_manual(values = c("#00008B", "#708090", "#8B0000")) +
    labs(x = "Log2FoldChange",
         #y = "-log10(Adjust P-Value)",
         #y = 'mRNA\n-log10(Adjust P-Value)',
         y = 'lncRNA\n-log10(Adjust P-Value)') +
    geom_hline(yintercept = -log10(0.05), 
               linetype = 'dotdash', size = 1) +
    geom_vline(xintercept = c(1, -1), 
               linetype = 'dotdash', color = 'grey30') +
    labs(subtitle = paste("DE analysis of ", title)) +
    scale_y_continuous(expand = c(0, 0)) +#,
    #limits = c(0, 45)) +
    #xlim(c(-15, 15)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 18),
          legend.title = element_blank(),
          legend.position = c(0.01, 0.85))
}



#===> get kegg crc-pathway genes
mygetenrichedgenefromkk <- function(genelist){
  gene.df <- bitr(genelist, fromType = "ENSEMBL",
                  toType = c("SYMBOL", "ENTREZID"),
                  OrgDb = org.Mm.eg.db)
  kk <- enrichKEGG(gene = gene.df$ENTREZID,
                   organism = 'mmu',
                   pvalueCutoff = 0.05)
  kk_res <- as.data.frame(kk)
  CRC_kk_entrezid <- dplyr::tbl_df(kk_res) %>%
    dplyr::filter(ID == 'mmu05210') %>%
    separate_rows(geneID, sep = '/') %>%
    pull(geneID)
  df <- dplyr::filter(
    gene.df, ENTREZID %in% CRC_kk_entrezid) 
  return(df) #ENSEMBL SYMBOL ENTREZID
}

mkdir('jobs')

#===>
cat('###################################################
=============>     Don\'t panic!      =============>
###################################################')




