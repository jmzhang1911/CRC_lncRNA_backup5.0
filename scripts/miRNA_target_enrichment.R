source('Utils.R')

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




