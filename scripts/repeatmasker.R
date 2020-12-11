
filepath <-list.files("Annolnc2_results/", full.names = T)

if(T){
  repeatpath <- paste(filepath, "/", basename(filepath),
                     "_trans1/repeat_elements/repeatmask_details.txt", sep = '')
  repeat_contain_list <- list()
  count = 0
  
  for(i in repeatpath){
    if(file.exists(i) & file.info(i)$size != 0){
      trans_id <- str_split(i, '/')[[1]][2]
      repeat_contain <- read_table2(i, col_names = FALSE) %>%
        dplyr::select(chr = X1,
                      start = X2,
                      end = X3,
                      strand = X6,
                      repeat_name = X4,
                      repeat_class = X11,
                      adjusted_SW_scores = X5)
      repeat_contain$trans_id <- trans_id
      repeat_contain_list[[trans_id]] = repeat_contain 
      count = count + 1
    }else{next}
  }
  
  repeat_contain_inte <- do.call(rbind, repeat_contain_list)
  rownames(repeat_contain_inte) <- NULL
  print(str_c('done===>',count)) # the maxtrans of 586 lncRNAs contains repeat sequence
}


# the maxtrans of 297 lncRNAs contains LTR
dplyr::filter(repeat_contain_inte, str_detect(repeat_class, 'LTR')) %>%
  group_by(trans_id) %>% summarise(count = n()) %>% arrange(desc(count))

# smmmary
tmpdf <- dplyr::filter(repeat_contain_inte) %>%
  group_by(repeat_class) %>% summarise(count = n()) %>% arrange(desc(count))
