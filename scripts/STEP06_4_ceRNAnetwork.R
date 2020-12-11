rm(list = ls())
source('Utils.R')


mkdir('outcomes/ceRNAnetwork/miRNA_DEanalysis')




#===> get the miRNA expression data from TCGA
if(F){
  data_category <- 'Transcriptome Profiling'
  data_type <- 'miRNA Expression Quantification'
  cancer_type <- 'TCGA-COAD'
  
  
  library(TCGAbiolinks)
  library(DT)
  library(SummarizedExperiment)
  miRNA_expr <- GDCquery(project = cancer_type,
                         data.category = data_category,
                         data.type = data_type,
                         workflow.type = "BCGSC miRNA Profiling")
  
  miRNA_expr_df <- GDCprepare(query = miRNA_expr)　
  save(miRNA_expr_df, file = 'outcomes/ceRNAnetwork/miRNA_DEanalysis/input.RData')
  
}




#===> miRNA DE analysis by DESeq2
load('outcomes/ceRNAnetwork/miRNA_DEanalysis/input.RData')
if(T){
  miRNA_expr <- miRNA_expr_df %>% column_to_rownames(var = 'miRNA_ID')
  
  miRNA_expr <- miRNA_expr[,seq(1, ncol(miRNA_expr), 3)]
  colnames(miRNA_expr) <- str_replace_all(colnames(miRNA_expr), 'read_count_', '')
  miRNA_expr <- miRNA_expr[apply(miRNA_expr, 1, function(x){
    sum(x) > 100}),]
}



grouplist <- if_else(as.numeric(str_sub(colnames(miRNA_expr), 14, 15)) < 10,'tumor', 'normal')
grouplist <- factor(grouplist, levels = c('tumor','normal'))
table(grouplist)

nor <- colnames(miRNA_expr)[grouplist == 'normal']
tum <- colnames(miRNA_expr)[grouplist =='tumor']
tmpexpr <- miRNA_expr[,c(nor, tum)]
tmpexpr %>% as.data.frame() %>% rownames_to_column(var = 'gene_id') %>%
  gather(key = 'sample', value = 'count', 2:(ncol(tmpexpr)+1)) %>%
  mutate(sample = factor(sample, levels = colnames(tmpexpr))) %>%
  ggplot(aes(x = sample, y = log10(count+1))) +
  geom_boxplot(fill = '#228B22') +
  theme_classic() +
  labs(x = '',
       title = 'Batch effect after normalization') +
  theme(axis.text.x = element_text(hjust = 1,
                                   vjust = 1,
                                   angle = 90)) 

cancer2 <- c('TCGA-DM-A1D9-01A-11H-A154-13',
             'TCGA-CM-6677-01A-11H-1838-13',
             'TCGA-CK-6748-01A-11H-1838-13',
             'TCGA-F4-6807-01A-11H-1838-13',
             'TCGA-DM-A285-01A-11H-A16S-13',
             'TCGA-A6-6650-01B-02R-A27D-13',
             'TCGA-DM-A1D7-01A-11H-A154-13',
             'TCGA-4T-AA8H-01A-11H-A41D-13')

tmpexpr2 <- miRNA_expr[, c(nor, cancer2)]

#=> DESeq2 analysis
colnames(tmpexpr2) <- paste(rep(c('normal','tumor'), each = 8), 1:8, sep = '_')
colData <- data.frame(row.names = colnames(tmpexpr2),
                      condition = rep(c('normal','cancer'), each = 8))

dds <- DESeq2::DESeqDataSetFromMatrix(countData = tmpexpr2,
                                      colData = colData,
                                      design = ~condition)
dds <- DESeq2::DESeq(dds)
res <- DESeq2::results(dds, contrast = c('condition', 'cancer', 'normal'))
nor_df <- counts(dds, normalized=TRUE)
nor_df %>% as.data.frame() %>% rownames_to_column(var = 'gene_id') %>%
  gather(key = 'sample', value = 'count', 2:(ncol(nor_df)+1)) %>%
  mutate(sample = factor(sample, levels = colnames(nor_df))) %>%
  ggplot(aes(x = sample, y = log10(count+1))) +
  geom_boxplot(fill = '#228B22') +
  theme_classic() +
  labs(x = '',
       title = 'Batch effect after normalization') +
  theme(axis.text.x = element_text(hjust = 1,
                                   vjust = 1,
                                   angle = 90)) 
DE_miRNA_res <- mydeanalysis(res)
save(DE_miRNA_res, file = 'outcomes/ceRNAnetwork/miRNA_DEanalysis/DE_miRNA_res.RData')




#===> check it out whether candidata miRNA is DE and survival related
load('outcomes/ceRNAnetwork/miRNA_DEanalysis/DE_miRNA_res.RData')
load('outcomes/ceRNAnetwork/lncRNA_miRNA_res/lncRNA_miRNA_res.RData')

dim(DE_miRNA_res) #839
dim(candidate_miRNA) #89
prob <- candidate_miRNA %>% mutate(prob = str_replace_all(prob, 'miR', 'mir'))
DE_candi_miRNA_df <- DE_miRNA_res %>% dplyr::filter(direction %in% c('up','down')) %>%
  mutate(prob = str_extract(gene_id, '(let|mir)-[0-9]+')) %>%
  left_join(prob, by = 'prob') %>% drop_na()

write.csv(DE_candi_miRNA_df, file = 'outcomes/ceRNAnetwork/DE_candi_miRNA_df.csv')




#===> expression levels of key 
prob <- c('hsa-mir-181a-2', 'hsa-mir-183','hsa-mir-141', 
          'hsa-mir-27a','hsa-mir-140', 'hsa-mir-128-1',
          'hsa-mir-212', 'hsa-mir-425', 'hsa-mir-7-3')

# visualization
nor_df %>% as.data.frame() %>% rownames_to_column(var = 'gene_id') %>%
  dplyr::filter(gene_id %in% prob) %>%
  gather(key = 'class', value = 'count', 2:(ncol(nor_df)+1)) %>%
  mutate(class = str_replace_all(class, '_[0-9]', ''),
         class = factor(class)) %>%
  ggplot(aes(x = class, y = log10(count+1))) +
  geom_boxplot(aes(fill = class)) +
  facet_wrap(gene_id~.) +
  scale_fill_nejm() +
  scale_y_continuous(limits = c(0,5)) +
  labs(x = '',
       y = paste("Expression level log10", "\n", ("(normalized reads count)"))) +
  geom_signif(comparisons = list(c('normal', 'tumor')),
              map_signif_level = T,
              test = t.test,
              y_position = 4.7) +
  theme_bw() +
  theme(axis.text.x = element_blank())

save(nor_df, DE_miRNA_df, file = 'outcomes/miRNA_DEanalysis/results.RData')



# visualization
dim(ggdata)
ggdata <- mydeanalysis(res) %>%
  drop_na()



ggdata2 <- ggdata %>% dplyr::filter(gene_id %in% prob)
dim(ggdata2)
myggvol <- function(df, title, df2){
  ggplot(df, aes(x = df[,2], y = -log10(df[,4]))) +
    geom_point(size = 1, 
               aes(color = df[,6]), 
               show.legend = T) +
    scale_color_manual(values = c("#00008B", "#708090", "#8B0000")) +
    labs(x = "Log2FoldChange",
         y = 'miRNA\n-log10(Adjust P-Value)') +
    geom_hline(yintercept = -log10(0.05), 
               linetype = 'dotdash', size = 1) +
    geom_vline(xintercept = c(1, -1), 
               linetype = 'dotdash', color = 'grey30') +
    labs(subtitle = paste("DE analysis of ", title)) +
    geom_text(data = df2, aes(label = gene_id)) +
    scale_y_continuous(expand = c(0, 0)) +#,
    #limits = c(0, 45)) +
    #xlim(c(-15, 15)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 18),
          legend.title = element_blank(),
          legend.position = c(0.71, 0.85))
}

library(ggrepel)
myggvol(ggdata, title = 'try', df2 =ggdata2)

ggplot(ggdata, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(size = 2, 
             aes(color = direction), 
             show.legend = T) +
  geom_point(data = ggdata2, size = 2, shape = 21) +
  geom_label_repel(data = ggdata2, aes(label = gene_id),
                   fontface="bold", color="black", 
                   segment.colour = "black",
                   box.padding = unit(1, "lines"),
                   point.padding = unit(0.8, "lines"),
                   force = T) +
  scale_color_manual(values = c("#00008B", "#708090", "#8B0000")) +
  labs(x = "Log2FoldChange",
       y = 'miRNA\n-log10(Adjust P-Value)') +
  geom_hline(yintercept = -log10(0.05), 
             linetype = 'dotdash', size = 1) +
  geom_vline(xintercept = c(1, -1), 
             linetype = 'dotdash', color = 'grey30') +
  labs(subtitle = paste("DE analysis of miRNA")) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        legend.title = element_blank(),
        legend.position = c(0.71, 0.85))


