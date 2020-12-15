source('Utils.R')

load('outcomes/inputdata/input.RData')
load('outcomes/candidate/cor_results.RData')


head(cor_results)
df <- dplyr::filter(cor_results, target %in% input_matrix_count$lnc_genelist, 
              source %!in% input_matrix_count$lnc_genelist,
              abs(r_value) >= 0.7, p_value < 0.05)

unique(df$target)
set.seed(123)
df_ <- group_by(df, target) %>% sample_n(1)
coding_prob <- data.frame(coding = df_$source)
lnc_prob <- data.frame(lnc = df_$target)

prob <- unique(c(df_$source,df_$target))
df2 <- dplyr::filter(input_matrix_count$mRNA_lncRNA_count, gene_id %in% prob) %>%
  column_to_rownames(var = 'gene_id') %>% myddnor() %>% as.data.frame() %>%
  rownames_to_column(var = 'gene_id') %>%
  pivot_longer(control_1:week10_3, names_to = 'time', values_to = 'expression') %>%
  mutate(time = str_remove_all(time, '_[0-9]'),
         time = factor(time, levels = c('control','week2','week4','week7','week10')))

coding_prob_df <- df2 %>% dplyr::filter(gene_id %in% coding_prob$coding) %>%
  left_join(coding_prob, by = c('gene_id' = 'coding')) %>% rename(coding = expression)

lnc_prob_df <- df2 %>% dplyr::filter(gene_id %in% lnc_prob$lnc) %>%
  left_join(lnc_prob, by = c('gene_id' = 'lnc')) %>% rename(lnc = expression)

ggdata <-  as.data.frame(cbind(coding_prob_df, lnc_prob_df)) %>% dplyr::select(coding, lnc)

ggplot(ggdata, aes(x = log10(coding+1), y = log10(lnc+1))) +
  geom_point() +
  geom_smooth(method = 'lm')



head(df)
head(df2)


