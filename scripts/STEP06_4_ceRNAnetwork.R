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
  #save(miRNA_expr_df, file = 'outcomes/ceRNAnetwork/miRNA_DEanalysis/input.RData')
  
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
colnames(tmpexpr2)
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
if(F){
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
}

DE_miRNA_res <- mydeanalysis(res)
#save(DE_miRNA_res, file = 'outcomes/ceRNAnetwork/miRNA_DEanalysis/DE_miRNA_res.RData')




#===> check it out whether candidata miRNA is DE and survival related
load('outcomes/ceRNAnetwork/miRNA_DEanalysis/DE_miRNA_res.RData')
load('outcomes/ceRNAnetwork/lncRNA_miRNA_res/new_lncRNA_miRNA_res.RData')

dim(DE_miRNA_res) #839
dim(new_candidate_miRNA) #89
prob <- new_candidate_miRNA %>% mutate(prob = str_replace_all(prob, 'miR', 'mir'))
DE_candi_miRNA_df <- DE_miRNA_res %>% dplyr::filter(direction %in% c('up','down')) %>%
  mutate(prob = str_extract(gene_id, '(let|mir)-[0-9]+')) %>%
  left_join(prob, by = 'prob') %>% drop_na()

write.csv(DE_candi_miRNA_df, file = 'outcomes/ceRNAnetwork/DE_candi_miRNA_df.csv')



#===> key miRNA family in ceRNA network
key_miRNA_family <- c('let-7-5p=miR-98-5p','miR-101-3p_1',
                      'miR-101a-3p_2=101b-3p_1=101b-3p_2',
                      'miR-103-3p=107-3p','miR-125-5p=351-5p',
                      'miR-128-3p','miR-132-3p=212-3p','miR-139-5p','miR-141-3p=200a-3p',
                      'miR-143-3p','miR-153-3p','miR-155-5p',
                      'miR-15-5p=16-5p=195-5p=322-5p=497-5p',
                      'miR-17-5p=20-5p=93-5p=106-5p','miR-181-5p',
                      'miR-182-5p','miR-183-5p','miR-190-5p',
                      'miR-192-5p=215-5p','miR-193-3p','miR-194-5p','miR-200bc-3p=429-3p',
                      'miR-21-5p','miR-212-5p','miR-218-5p','miR-219-5p','miR-223-3p',
                      'miR-24-3p','miR-25-3p=32-5p=92-3p=363-3p=367-3p','miR-26-5p',
                      'miR-27-3p','miR-338-3p','miR-34-5p=449-5p','miR-375-3p',
                      'miR-425-5p=489-3p','miR-7-5p','miR-10-5p')

#miRNA_family => hsa_mir
hsa_mir <- c('hsa-let-7b','hsa-let-7d','hsa-let-7f-1','hsa-let-7f-2',
            'hsa-let-7g','hsa-mir-101-1','hsa-mir-103a-1','hsa-mir-103a-2',
            'hsa-mir-125b-1','hsa-mir-128-2','hsa-mir-132','hsa-mir-139',
            'hsa-mir-143','hsa-mir-153-1','hsa-mir-155','hsa-mir-15a','hsa-mir-15b',
            'hsa-mir-17','hsa-mir-181a-2','hsa-mir-181c','hsa-mir-181d',
            'hsa-mir-182','hsa-mir-183','hsa-mir-190a','hsa-mir-192',
            'hsa-mir-193a','hsa-mir-193b','hsa-mir-194-1','hsa-mir-194-2',
            'hsa-mir-200b','hsa-mir-21','hsa-mir-212','hsa-mir-218-1','hsa-mir-141',
            'hsa-mir-218-2','hsa-mir-219a-1','hsa-mir-223','hsa-mir-24-1',
            'hsa-mir-24-2','hsa-mir-25','hsa-mir-26a-1','hsa-mir-26a-2',
            'hsa-mir-26b','hsa-mir-27a','hsa-mir-27b','hsa-mir-338',
            'hsa-mir-34a','hsa-mir-375','hsa-mir-425','hsa-mir-7-1','hsa-mir-7-2','hsa-mir-7-3',
            'has-mir-10a','has-mir-10b')
save(key_miRNA_family, file = 'outcomes/ceRNAnetwork/key_miRNA_family.RData')


ggdata <- DE_candi_miRNA_df %>% dplyr::filter(miRNA_family %in% key_miRNA_family,
                                              gene_id %in% hsa_mir) 



prob <- c('hsa-mir-181a-2', 'hsa-mir-17','hsa-mir-143',
          'hsa-mir-7-2','hsa-mir-375','hsa-mir-155',
          'hsa-mir-34a','hsa-mir-125b-1','hsa-mir-26b')

set.seed(123)
prob <- sample(hsa_mir, 16)

# visualization
nor_df %>% as.data.frame() %>% rownames_to_column(var = 'gene_id') %>%
  dplyr::filter(gene_id %in% prob) %>%
  gather(key = 'class', value = 'count', 2:(ncol(nor_df)+1)) %>%
  mutate(class = str_replace_all(class, '_[0-9]', ''),
         class = factor(class)) %>%
  ggplot(aes(x = class, y = log10(count+1))) +
  geom_boxplot(aes(fill = class)) +
  facet_wrap(gene_id~.,scales = 'free') +
  scale_fill_nejm() +
  scale_y_continuous(limits = c(0,7)) +
  labs(x = '',
       y = paste("Expression level log10", "\n", ("(normalized reads count)"))) +
  geom_signif(comparisons = list(c('normal', 'tumor')),
              map_signif_level = T,
              test = t.test,
              
              y_position = 6.5) +
  theme_bw() +
  theme(axis.text.x = element_blank())
ggsave('outcomes/ceRNAnetwork/hsa_mir_DEplot_prob.pdf', height = 8, width = 8, dpi = 300)



ggdata2 <- mydeanalysis(res) %>% drop_na()
ggdata3 <- ggdata2 %>% dplyr::filter(gene_id %in% hsa_mir)

ggplot(ggdata2, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(size = 2, color = '#708090', show.legend = T, alpha = 0.5) +
  geom_point(data = ggdata3, size = 3, #shape = 21, 
             aes(color = direction)) +
  scale_color_manual(values = c("#00008B", "#8B0000")) +
  labs(x = "Log2FoldChange",
       y = 'miRNA\n-log10(Adjust P-Value)') +
  geom_hline(yintercept = -log10(0.05), 
             linetype = 'dotdash', size = 1) +
  geom_vline(xintercept = c(1, -1), 
             linetype = 'dotdash', color = 'grey30') +
  labs(subtitle = paste("DE analysis of miRNA")) +
  scale_y_continuous(expand = c(0, 0),limits = c(0,120)) +
  theme_minimal(base_size = 17) +
  theme(plot.title = element_text(hjust = 1,size = 2),
        legend.title = element_blank(),
        legend.position = c(0.2, 0.85))
ggsave('outcomes/ceRNAnetwork/hsa_mir_volcano.pdf',height = 10, width = 8, dpi = 300)























