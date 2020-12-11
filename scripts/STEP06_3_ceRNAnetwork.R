source('Utils.R')


mkdir('outcomes/ceRNAnetwork/survival')

#===> TCGA CRC miRNA expression and clinical data
library(XML)
library(methods)


filepath <- 'data/TCGA_miRNA_clinical/clinicaldata'
filelist <- list.files(path = filepath, pattern = '*.xml$', recursive = T)
tmp <- lapply(filelist, function(x){
  result <- xmlParse(file = file.path(filepath, x))
  rootnode <- xmlRoot(result)
  xmldataframe <- xmlToDataFrame(rootnode[2])
  return(t(xmldataframe))
})

# some file present wrong format, just delete them
if(T){
  rownames(tmp[[148]]);tmp[[148]] <- NULL
  rownames(tmp[[148]]);tmp[[148]] <- NULL
  rownames(tmp[[343]]);tmp[[343]] <- NULL 
} 

clinical_df <- t(do.call(cbind, tmp))
#save(clinical_df, file = 'data/TCGA_miRNA_clinical/clinical_df.RData')
load('data/TCGA_miRNA_clinical/clinical_df.RData')



#===> miRNA expression data
filepath <- 'data/TCGA_miRNA_clinical/expressiondata/'
filelist <- list.files(path = filepath, pattern = '*.quantification.txt$', recursive = T)

if(T){
  library(rjson)
  result <- fromJSON(file = "data/TCGA_miRNA_clinical/CRC_miRNA.json")
  fls <- unlist(lapply(result,function(x){x[[3]]}))
  cid <- unlist(lapply(result,function(x){x[[2]][[1]][[2]]}))
  id2fls <- data.frame(cid=cid,fls=fls)
  head(id2fls) 
}

miRNA_exprlist <- lapply(filelist, function(x){
  tmp <- read.table(file = file.path(filepath, x), sep = '\t',
                    header = T)[,1:2]
  colnames(tmp)[2] <- basename(x)
  return(tmp)
})

# filter impor messages
if(T){
  miRNA_expr <- t(do.call(cbind, miRNA_exprlist))
  colnames(miRNA_expr) <- miRNA_expr[1,]
  miRNA_expr[1:5,1:5]
  miRNA_expr <- miRNA_expr[seq(2, nrow(miRNA_expr), 2),]
  dim(miRNA_expr)
  miRNA_expr <- t(miRNA_expr)
}


# match the id of expression matrix by id2fls
if(T){
  df <- data.frame(id = filelist)
  match(colnames(miRNA_expr), id2fls$fls)
  miRNA_expr2 <- miRNA_expr[, match(id2fls$fls, colnames(miRNA_expr))]
  colnames(miRNA_expr2) <- id2fls$cid
  miRNA_expr3 <- miRNA_expr2 %>% as.data.frame() %>% rownames_to_column(var = 'gene_id') %>%
    dplyr::select(-gene_id) %>% 
    mutate_if(is.character, as.numeric) 
  rownames(miRNA_expr3) <- rownames(miRNA_expr2)
  dim(miRNA_expr3)
  miRNA_expr <- miRNA_expr3[apply(miRNA_expr3, 1, function(x){sum(x>1)>10}),]
  dim(miRNA_expr)}


# survival meta data
if(T){
  meta <- as.data.frame(clinical_df[,c('bcr_patient_barcode',
                                       'bcr_patient_uuid',
                                       'vital_status',
                                       'days_to_death',
                                       'days_to_last_followup',
                                       'race_list',
                                       'age_at_initial_pathologic_diagnosis',
                                       'gender',
                                       'stage_event')])
  df <- dplyr::select(meta, bcr_patient_uuid,
                      days_to_death,
                      days_to_last_followup,
                      vital_status) %>%
    mutate(os = if_else(vital_status == 'alive',
                        as.numeric(days_to_last_followup),
                        as.numeric(days_to_death)),
           vital_status = if_else(vital_status == 'alive', 1, 2))
  
  tmp <- apply(miRNA_expr, 1, function(x){
    tmp = median(x) 
    if_else(x>=tmp,'high','low')})
  tmp <- as.data.frame(t(tmp))
  rownames(tmp) <- rownames(miRNA_expr)
  colnames(tmp) <- colnames(miRNA_expr)
  condi_matrix <- tmp
  save(miRNA_expr, df, condi_matrix, 
       file = 'outcomes/ceRNAnetwork/survival/Survival_data.RData')}




#===> plot survival
###########=====> do it in jobs >>>>>>###########
#===> saved as miRNA_survival_plot.R in jobs
###########=====> do it in jobs >>>>>>###########


