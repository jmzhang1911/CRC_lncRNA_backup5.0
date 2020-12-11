source('Utils.R')


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
save(miRNA_expr_df, file = 'outcomes/miRNA_DEanalysis/input.RData')
miRNA_expr <- miRNA_expr_df %>% column_to_rownames(var = 'miRNA_ID')

miRNA_expr <- miRNA_expr[,seq(1, ncol(miRNA_expr), 3)]
colnames(miRNA_expr) <- str_replace_all(colnames(miRNA_expr), 'read_count_', '')
miRNA_expr <- miRNA_expr[apply(miRNA_expr, 1, function(x){
  sum(x) >100
}),]

dim(miRNA_expr)
group_list <- if_else(as.numeric(str_sub(colnames(miRNA_expr), 14, 15)) < 10,'tumor', 'normal')
group_list <- factor(group_list, levels = c('tumor','normal'))
table(group_list)
a <- colnames(miRNA_expr)[group_list == 'normal']
b <- colnames(miRNA_expr)[group_list =='tumor']
tmpexpr <- miRNA_expr[,b]
tmpexpr2 <- miRNA_expr[,c(a,cancer2)]
tmpexpr3 <- miRNA_expr[,]
boxplot(log10(tmpexpr+1))
boxplot(log10(tmpexpr2+1))

colnames(tmpexpr2) <- c(paste(rep('normal', 8),1:8,sep = '_'),
                        paste(rep('cancer', 8), 1:8, sep = '_'))
colData <- data.frame(row.names = colnames(tmpexpr2),
                      condition = c(rep('normal',8),rep('cancer',8)))



colnames(tmpexpr2) <- paste(rep(c('normal','cancer'), each = 4), 1:4, sep = '_')
colData <- data.frame(row.names = colnames(tmpexpr2),
                      condition = rep(c('normal','cancer'), each = 4))

dds <- DESeq2::DESeqDataSetFromMatrix(countData = tmpexpr2,
                                      colData = colData,
                                      design = ~condition)
dds <- DESeq2::DESeq(dds)
res <- DESeq2::results(dds, contrast = c('condition', 'normal', 'cancer'))
nor_df <- counts(dds, normalized=TRUE)
boxplot(log10(nor_df+1))
DE_miRNA_df <- mydeanalysis(res) %>% 
  dplyr::filter(direction %in% c('up','down')) 


df <- data.frame(
  df = colnames(miRNA_expr)
)
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
ggsave(p, filename = 'Batch_effect.png')


cancer2 <- c('TCGA-DM-A1D9-01A-11H-A154-13',
             'TCGA-CM-6677-01A-11H-1838-13',
             'TCGA-CK-6748-01A-11H-1838-13',
             'TCGA-F4-6807-01A-11H-1838-13',
             'TCGA-DM-A285-01A-11H-A16S-13',
             'TCGA-A6-6650-01B-02R-A27D-13',
             'TCGA-DM-A1D7-01A-11H-A154-13',
             'TCGA-4T-AA8H-01A-11H-A41D-13')

tmpexpr2 <- miRNA_expr[,cancer2]


prob <- c('hsa-mir-181a-2','hsa-mir-183',
          'hsa-mir-193a','hsa-mir-141',
          'hsa-mir-140','hsa-mir-128-1',
          'hsa-mir-212','hsa-mir-129-1',
          'hsa-mir-7-2','hsa-mir-27a',
          'hsa-mir-425')


ggdata <- nor_df %>% as.data.frame() %>% rownames_to_column(var = 'gene_id') %>%
  dplyr::filter(gene_id %in% prob) %>%
  gather(key = 'class', value = 'count', 2:(ncol(nor_df)+1)) %>%
  #mutate(class = factor(class))
  mutate(class = str_replace_all(class, '_[0-9]', ''),
         class = factor(class))


ggplot(ggdata, aes(x = class, y = log10(count+1))) +
  geom_boxplot(aes(fill = class)) +
  facet_wrap(gene_id~.) +
  scale_fill_nejm() +
  scale_y_continuous(limits = c(0,5)) +
  labs(x = '',
       y = paste("Expression level log10", "\n", ("(normalized reads count)"))) +
  geom_signif(comparisons = list(c('normal', 'cancer')),
              map_signif_level = T,
              test = t.test,
              y_position = 4.7) +
  theme_bw() +
  theme(axis.text.x = element_blank())





