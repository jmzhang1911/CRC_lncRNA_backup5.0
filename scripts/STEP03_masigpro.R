rm(list = ls())
source('Utils.R')
library(maSigPro)


mkdir('outcomes/masigpro')
load('outcomes/inputdata/input.RData')




#===> input data
mRNA_lncRNA_count <-  input_matrix_count$mRNA_lncRNA_count %>%
  column_to_rownames(var = 'gene_id')
cpm <- edgeR::cpm(mRNA_lncRNA_count)




#===> maSigPro analysis
if(T){
  Time <- c(rep(c(0, 2, 4, 7, 10), each = 3))
  Replicate <- c(rep(1:5, each = 3))
  Control <- c(rep(0, 15))
  Treat <- c(rep(1, 15))
  sample_infor <- data.frame(Time, Replicate, Control, Treat)
  row.names(sample_infor) <- colnames(mRNA_lncRNA_count)
}

if(T){
  design <- make.design.matrix(sample_infor, degree = 4)
  design$edesign
  design$groups.vector
  df = data.matrix(cpm)
  fit <- p.vector(df, design, Q = 0.05, MT.adjust = "BH", min.obs = 6)
  tstep <- T.fit(fit, step.method = "backward", alfa = 0.05)
  save(tstep, file = 'outcomes/masigpro/tstep.RData') 
}

load('outcomes/masigpro/tstep.RData')
get<-get.siggenes(tstep, rsq = 0.3, vars="all")

if(T){
  pdf("outcomes/masigPro/maSigPro_result.pdf")
  cluster_result = see.genes(get$sig.genes,
                             k = 9, newX11 = FALSE)
  save(cluster_result, file = "outcomes/masigpro/cluster_result.RData")
  dev.off()
}

#===> get all the dynamic genes(coding and lnc)
load('outcomes/masigpro/cluster_result.RData')
dynamicdf <- data.frame(
  gene_id = names(cluster_result$cut),
  cluster = cluster_result$cut,
  row.names = NULL
) %>% mutate(gene_type = if_else(gene_id %in% input_matrix_count$lnc_genelist,
                                            'lnc', 'coding'))

table(dynamicdf$cluster)
save(dynamicdf, file = 'outcomes/masigpro/dynamicdf.RData')



#===> 
a <- column_to_rownames(input_matrix_count$fpkm_count, var = 'gene_id')[dynamicdf$gene_id,]
rownames(a) == dynamicdf$gene_id
anno_row <- dynamicdf %>% group_by(cluster, gene_type) %>% summarise(count = n()) %>%
  pivot_wider(names_from = gene_type, values_from = count) %>%
  mutate(total = coding + lnc,
         Cluster = str_c('cluster',cluster,' ',total,' = ',
                      coding,'(coding)', ' + ', lnc, '(lncRNA)')) %>%
  left_join(dynamicdf, by = 'cluster') %>%
  ungroup() %>%
  dplyr::select(gene_id, Cluster) %>% column_to_rownames(var = 'gene_id')
a <- a[match(rownames(anno_row), rownames(a)),]
#rownames(a) <- rownames(anno_row)

rownames(anno_row) == rownames(a)

anno_col= data.frame(
  Group = rep(c('Ctrl','Week2','Week4','Week7','Week10'), each = 3))
    
anno_col$Group <- factor(anno_col$Group, levels = c("Ctrl", "Week2", "Week4", "Week7", "Week10"))
row.names(anno_col) <- colnames(a)

pdf('outcomes/masigpro/heatmap.pdf', width = 8, height = 6)
pheatmap(a,
         cluster_rows = F,
         show_rownames = F,
         cluster_cols = F,
         scale = "row",
         annotation_row = anno_row,
         color =colorRampPalette(c("#87CEFA", "white", "#CC2121"))(100),
         annotation_col = anno_col)
dev.off()

