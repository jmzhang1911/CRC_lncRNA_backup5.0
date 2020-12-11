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

