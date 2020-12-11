rm(list = ls())
source('Utils.R')
mkdir('outcomes/pca_heatmap')


load('outcomes/inputdata/input.RData')
mRNA_count <- input_matrix_count$mRNA_count
lncRNA_count <- input_matrix_count$lncRNA_count
all_count <- input_matrix_count$mRNA_lncRNA_count



###########=====> do it in jobs >>>>>>###########
#===> saved as kendall_calculate.R in jobs
###########=====> do it in jobs >>>>>>###########



#===> plot heatmap
load('outcomes/pca_heatmap//pheatmap_input.RData')
myplotheatmap <- function(df, sample_infor, title){
  p <- pheatmap(df,
                scale = 'row',
                show_rownames = F,
                show_colnames = F,
                cluster_cols = F,
                cluster_rows = F,
                color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                annotation_col = sample_infor,
                annotation_row = sample_infor,
                main = title
                #gaps_row = c(3,9),
                #gaps_col = c(3,9)
  )
  #require(ggplotify)
  p2 <- as.ggplot(p)
  return(p2)
}

plot_mRNA <- myplotheatmap(mRNA_df, sample_infor, 'heatmap of mRNA')
plot_lncRNA <- myplotheatmap(lncRNA_df, sample_infor, 'heatmap of lncRNA')
plot_all <- myplotheatmap(all_df, sample_infor, 'heatmap of all')

p <- ((plot_mRNA | plot_lncRNA | plot_all )) +
  plot_layout(guides = 'collect',nrow = 1) +
  plot_annotation(tag_levels = "A") 
ggsave('outcomes/pca_heatmap/heatmap.png', width = 60, height = 20, units = "cm")


#===> PCA analysis and plot
mygetpacdata <- function(x, pic=F){
  
  # use pic to check the xlabs and ylabs
  
  time <- c(paste('control', 1:3, sep = '_'),
            paste(rep(c('inflammation','cancer'), each = 6), 1:6, sep = '_'))
  colnames(x) <- time
  time2 <- c(rep('control',3), rep('inflammation', 6), rep('cancer', 6)) 
  sample_infor <- data.frame(time = factor(time2))
  rownames(sample_infor) <- colnames(x)
  
  pca <- pca(x, metadata = sample_infor)
  biplot(pca, x = 'PC1', y = 'PC2')
  
  pca_rlt <- rownames_to_column(pca$rotated, var = "sample_name")
  pca_sample <- rownames_to_column(sample_infor, var = "sample_name")
  pca_plot_data <- full_join(pca_rlt, pca_sample, by = 'sample_name')
  
  if(pic == F){
    return(pca_plot_data)
  }else{
    return(biplot(pca, x = 'PC1', y = 'PC2'))
  }
}


myplotpca <- function(df, title, x, y){
  pcadata <- mygetpacdata(df, pic = F)
  pca <- ggplot(data = pcadata, aes(x = PC1, y = PC2)) +
    geom_point(size = 5,
               aes(color = time)) +
    stat_ellipse(aes(color = time)) +
    scale_shape_manual(values = range(c(22, 24))) +
    scale_color_manual(values = c("#00008B", "#708090", "#8B0000")) +
    labs(title = title,
         x = str_c('PCA1 (',x,'% variance explained)'),
         y = str_c('PCA2 (',y,'% variance explained)')) +
    theme_half_open() +
    scale_fill_brewer(palette = 'Set3') +
    theme(#legend.position = c(0.85, 0.2),
      plot.title = element_text(size = 18, hjust = 0.5)) +
    guides(fill = guide_legend(override.aes = list(shape = 21)))
}


pca_mRNA <- myplotpca(mRNA_expr,title ='PCA of mRNA',x='60.59',y='22.23')
pca_lncRNA <- myplotpca(lncRNA_expr,title='PCA of lncRNA',x='91.24',y='3.16')
pca_all <- myplotpca(all_expr,title='PCA of all genes',x='60.64',y='22.25')

p2 <- ((pca_mRNA | pca_lncRNA | pca_all)) +
  plot_layout(guides = 'collect',nrow = 1) +
  plot_annotation(tag_levels = "A") 
ggsave('outcomes/pca_heatmap/pac.png', width = 60, height = 20, units = "cm")

