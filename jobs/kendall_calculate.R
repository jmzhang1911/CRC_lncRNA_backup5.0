source('Utils.R')


#===> calculate correlation by kendall
if(T){
  load('outcomes/inputdata/input.RData')
  mRNA_count <- input_matrix_count$mRNA_count
  lncRNA_count <- input_matrix_count$lncRNA_count
  all_count <- input_matrix_count$mRNA_lncRNA_count
  
  mRNA_expr <- myddnor(mRNA_count %>% column_to_rownames(var = 'gene_id'))
  lncRNA_expr <- myddnor(lncRNA_count %>% column_to_rownames(var = 'gene_id'))
  all_expr <- myddnor(all_count %>% column_to_rownames(var = 'gene_id'))
  
  sample_infor <- data.frame(
    time = rep(c('control','week2','week4','week7','week10'),each = 3))
  rownames(sample_infor) <- colnames(mRNA_expr)
  sample_infor$time <- factor(sample_infor$time, 
                              levels = c('control','week2','week4','week7','week10'))
  
  mRNA_df <- cor(mRNA_expr, method = 'kendall')
  print('===> mRNA_df done')
  
  lncRNA_df <- cor(lncRNA_expr, method = 'kendall')
  print('===> lncRNA_df done')
  
  all_df <- cor(all_expr, method = 'kendall')
  print('===> all_df done')
  
  save(sample_infor, mRNA_df, 
       mRNA_expr, lncRNA_expr,all_expr
       lncRNA_df, all_df,
       file = 'outcomes/pca_heatmap/pheatmap_input.RData')
}

