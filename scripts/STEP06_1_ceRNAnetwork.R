rm(list = ls())
source('Utils.R')


load('outcomes/inputdata/input.RData')
filepath <-list.files("Annolnc2_results/", full.names = T)
mkdir('outcomes/ceRNAnetwork/lncRNA_miRNA_res')




#===> miRNA interaction
if(T){
  miRNApath <- paste(filepath, "/", basename(filepath),
                     "_trans1/miRNA_interaction/prediction.out", sep = '')
  miRNA_inter_list <- list()
  count = 0
  
  for(i in miRNApath){
    if(file.exists(i) & file.info(i)$size != 0){
      trans_id <- str_split(i, '/')[[1]][2]
      miRNA_inter <- read_delim(i, '\t', escape_double = FALSE, trim_ws = TRUE) 
      colnames(miRNA_inter)[1:4] <- c('lnc_trans','miRNA_family','match_start','match_end')
      miRNA_inter <- miRNA_inter%>%
        mutate(lnc_trans = str_replace_all(lnc_trans, '_trans1', '')) %>%
        dplyr::select(lnc_trans:match_end)
      miRNA_inter_list[[trans_id]] = miRNA_inter 
      count = count + 1
    }else{next}
  }
  
  miRNA_inter_inte <- do.call(rbind, miRNA_inter_list)
  rownames(miRNA_inter_inte) <- NULL
  print(str_c('done===>',count)) # the maxtrans of 681 lncRNAs can interact with miRNA
}

if(T){
  miRNApath3 <- paste(filepath, "/secondary_structure/seq_struct.dbn", sep = '')
  miRNA_seq_list <- list()
  count = 0
  
  for(i in miRNApath3){
    trans_id <- str_split(i, '/')[[1]][2]
    seqdf <-read.table(i, sep = '\t')
    new_seqdf <- data.frame(lnc_trans = seqdf$V1[1],
                            pos = seqdf$V1[3]) %>%
      mutate(lnc_trans = str_replace_all(lnc_trans, '>',''))
    miRNA_seq_list[[trans_id]] <- new_seqdf
    count = count + 1
  }
  miRNA_seq_inte <- do.call(rbind, miRNA_seq_list)
  rownames(miRNA_seq_inte) <- NULL
  print(str_c('done===>',count))
}

if(T){
  miRNA_inter_inte <- left_join(miRNA_inter_inte, miRNA_seq_inte, by = 'lnc_trans') %>%
    mutate(loop_seq = substr(pos, match_start, match_end),
           loop = if_else(str_count(loop_seq, "\\(|\\)") >= 4, 0, 1)) 
    
  miRNA_inter_reslist <- list(
    miRNA_inter_inte = miRNA_inter_inte,
    miRNA_inter_each = miRNA_inter_list,
    miRNA_seq_inte = miRNA_seq_inte,
    miRNA_seq_each = miRNA_seq_list)
  save(miRNA_inter_reslist, 
       file = 'outcomes/ceRNAnetwork/lncRNA_miRNA_res/miRNA_inter_reslist.RData')
}










