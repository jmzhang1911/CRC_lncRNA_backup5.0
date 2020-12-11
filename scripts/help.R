theme_classic() +
labs(x = '',
title = 'Batch effect after normalization') +
theme(axis.text.x = element_text(hjust = 1,
vjust = 1,
angle = 90))}
View(miRNA_expr)
grouplist <- if_else(as.numeric(str_sub(colnames(miRNA_expr), 14, 15)) < 10,'tumor', 'normal')
colnames(miRNA_expr) <- grouplist
colData <- data.frame(row.names = colnames(miRNA_expr),
condition = grouplist)
miRNA_expr <- miRNA_expr_df %>% column_to_rownames(var = 'miRNA_ID')
miRNA_expr <- miRNA_expr[,seq(1, ncol(miRNA_expr), 3)]
colnames(miRNA_expr) <- str_replace_all(colnames(miRNA_expr), 'read_count_', '')
# filter
miRNA_expr <- miRNA_expr[apply(miRNA_expr, 1, function(x){sum(x) >20}),]
grouplist <- if_else(as.numeric(str_sub(colnames(miRNA_expr), 14, 15)) < 10,'tumor', 'normal')
grouplist <- factor(grouplist, levels = c('normal','tumor'))
grouplist <- factor(grouplist, levels = c('tumor','normal'))
colData <- data.frame(row.names = colnames(miRNA_expr),
condition = grouplist)
#condition = rep(c('normal','tumor'), each = 8))
dds <- DESeq2::DESeqDataSetFromMatrix(countData = miRNA_expr,
colData = colData,
design = ~condition)
dds <- DESeq2::DESeq(dds)
res <- DESeq2::results(dds, contrast = c('condition', 'tumor', 'normal'))
DE_miRNA_df <- mydeanalysis(res) %>%
dplyr::filter(direction %in% c('up','down')) %>%
mutate(class = if_else(gene_id %in% nor, 'normal', 'tumor'))
prob <- 'hsa-mir-181b-2'
miRNA_expr_nordf %>% rownames_to_column(var = gene_id) %>%
dplyr::filter(gene_id == prob)
miRNA_expr_nordf %>% as.data.frame() %>% rownames_to_column(var = gene_id) %>%
dplyr::filter(gene_id == prob)
miRNA_expr_nordf %>% as.data.frame() %>% rownames_to_column(var = 'gene_id') %>%
dplyr::filter(gene_id == prob)
miRNA_expr_nordf %>% as.data.frame() %>% rownames_to_column(var = 'gene_id') %>%
dplyr::filter(gene_id == prob) %>%
gather(key = 'class', value = 'count', 2:(ncol(miRNA_expr_nordf)+1))
miRNA_expr_nordf %>% as.data.frame() %>% rownames_to_column(var = 'gene_id') %>%
dplyr::filter(gene_id == prob) %>%
gather(key = 'class', value = 'count', 2:(ncol(miRNA_expr_nordf)+1)) %>%
mutate(class = if_else(class %in% nor, 'normal', 'tumor'))
miRNA_expr_nordf %>% as.data.frame() %>% rownames_to_column(var = 'gene_id') %>%
dplyr::filter(gene_id == prob) %>%
gather(key = 'class', value = 'count', 2:(ncol(miRNA_expr_nordf)+1)) %>%
mutate(class = if_else(class %in% nor, 'normal', 'tumor'),
class = factor(class, levels = c('normal','tumor')))
ggdata <- miRNA_expr_nordf %>% as.data.frame() %>% rownames_to_column(var = 'gene_id') %>%
dplyr::filter(gene_id == prob) %>%
gather(key = 'class', value = 'count', 2:(ncol(miRNA_expr_nordf)+1)) %>%
mutate(class = if_else(class %in% nor, 'normal', 'tumor'),
class = factor(class, levels = c('normal','tumor')))
ggplot(ggdata, aes(x = class, y = log10(count+1))) +
geom_boxplot(fill = '#228B22') +
theme_classic() +
labs(x = '',
title = 'Batch effect after normalization') +
theme(axis.text.x = element_text(hjust = 1,
vjust = 1,
angle = 90))
ggplot(ggdata, aes(x = class, y = log10(count+1))) +
geom_boxplot(aes(fill = class)) +
theme_minimal() +
theme(axis.text.x = element_text(hjust = 1,
vjust = 1))
ggplot(ggdata, aes(x = class, y = log10(count+1))) +
geom_boxplot(aes(fill = class)) +
theme_minimal() +
labs(x = '') +
theme(axis.text.x = element_text(hjust = 1,
vjust = 1))
ggplot(ggdata, aes(x = class, y = log10(count+1))) +
geom_boxplot(aes(fill = class)) +
theme_minimal() +
labs(x = '') +
geom_signif(comparisons = list(c('normal','tumor')),
map_signif_level = T) +
theme(axis.text.x = element_text(hjust = 1,
vjust = 1))
ggplot(ggdata, aes(x = class, y = log10(count+1))) +
geom_boxplot(aes(fill = class)) +
theme_minimal() +
labs(x = '') +
geom_signif(comparisons = list(c('normal','tumor')),
map_signif_level = T,
t.test =T) +
theme(axis.text.x = element_text(hjust = 1,
vjust = 1))
?geom_signif
ggplot(ggdata, aes(x = class, y = log10(count+1))) +
geom_boxplot(aes(fill = class)) +
theme_minimal() +
labs(x = '') +
geom_signif(comparisons = list(c('normal','tumor')),
map_signif_level = T,
test.args = 't.test') +
theme(axis.text.x = element_text(hjust = 1,
vsource('Utils.R')


#===> crc relative miRNA in hsa
hsa_miRNA_crc <- list.files('outcomes/survival/crc_miRNA/', pattern = '*.png',
                            full.names = F) %>% str_replace_all(' _survival_analysis.png', '')
write.table(hsa_miRNA_crc, file = 'outcomes/hsa_miRNA_crc.txt', 
            col.names = F, row.names = F, quote = F)



#===> get miRNA seq from targetscan
if(!file.exists('outcomes/inte_miRNA/tarscan_seq.txt')){
  miRNA_scan_df <- read.table('data/miR_Family_Info.txt', sep = '\t', header = T) %>% 
    dplyr::filter(Species.ID == c(9606, 10090)) 
  mmu_miRNA_scan_seq <- miRNA_scan_df %>%
    dplyr::filter(Species.ID == 10090) %>%
    dplyr::select(miR.family, Species.ID, MiRBase.ID, Mature.sequence, MiRBase.Accession) 
  for(i in 1:length(mmu_miRNA_scan_seq$miR.family)){
    ID <- mmu_miRNA_scan_seq[i,]$miR.family
    Seq <- mmu_miRNA_scan_seq[i,]$Mature.sequence
    cat(paste('>', ID,'\n', sep = ''),
        file = 'outcomes/inte_miRNA/tarscan_seq.txt', append = T)
    cat(paste(Seq,'\n', sep = ''), 
        file = 'outcomes/inte_miRNA/tarscan_seq.txt', append = T)
  }
}



#===> TCGA CRC miRNA expression and clinical data
library(XML)
library(methods)


filepath <- 'TCGA_miRNA_clinical/clinicaldata'
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
save(clinical_df, file = 'TCGA_miRNA_clinical/clinical_df.RData')
load('TCGA_miRNA_clinical/clinical_df.RData')



#===> miRNA expression data
filepath <- filepath <- 'TCGA_miRNA_clinical/expressiondata/'
filelist <- list.files(path = filepath, pattern = '*.quantification.txt$', recursive = T)

if(T){
  library(rjson)
  result <- fromJSON(file = "TCGA_miRNA_clinical/CRC_miRNA.json")
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
}


# match the id of expression matrix by id2fls
if(T){
  df <- data.frame(id = filelist)
  match(rownames(miRNA_expr), id2fls$fls)
  miRNA_expr2 <- miRNA_expr[match(id2fls$fls, rownames(miRNA_expr)),]
  rownames(miRNA_expr2) <- id2fls$cid
  miRNA_expr3 <- miRNA_expr2 %>% as.data.frame() %>% rownames_to_column(var = 'gene_id') %>%
    dplyr::select(-gene_id) %>% 
    mutate_if(is.character, as.numeric) 
  rownames(miRNA_expr3) <- rownames(miRNA_expr2)
  dim(miRNA_expr3)
  miRNA_expr <- miRNA_expr3[apply(miRNA_expr3, 1, function(x){sum(x>1)>10}),]
  dim(miRNA_expr)
}


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
                        %PDF-1.4
%âãÏÓ\r
1 0 obj
<<
/CreationDate (D:20201124122230)
/ModDate (D:20201124122230)
/Title (R Graphics Output)
/Producer (R 4.0.2)
/Creator (R)
>>
endobj
2 0 obj
<< /Type /Catalog /Pages 3 0 R >>
endobj
7 0 obj
<< /Type /Page /Parent 3 0 R /Contents 8 0 R /Resources 4 0 R >>
endobj
8 0 obj
<<
/Length 334461 /Filter /FlateDecode
>>
stream
xœä½K-Kr9¿¿"‡u<r÷x¸û”„$€€X€$Ñ¨+Aìª¢X—”ôó{{„?ì[±·G$‹@zÀKìÊ<y§?ÌlÙZËü×_~ù¯¿ÿúÇŸşã×?~-éG¯ÿ®?Öøå]şáÒWr?¶ıë¿|ıç¯?üôo~ıOÿşÏ¿şâ¯~r?œs_ö¿õÿ×ëÛ×ÿşé¯ÿöË}ı×Ÿü×_¾şïïòå¾şÃOëòÃù¯%‡›ÿúıO›ûÂñqY¿~÷Ó¬Ëë£ÿ‘Ãëã¾şØ·×ÇåÇš_£ÿ‘öã›×òÕxüÓÊÇP¾šöòÿË7/éõñõÿ7w|uß_½{ı»Âñå}9>Ç¹üU¯enå³ß~øíøìbù^ÿ¤ıüqùøœlñøz*ÿR¿ì?b>>oÇŸ_×²À%?æw?ıUù]şi¿ª¶±ü_÷eİÒôú¾şêòõõGğíãøúêü¸ô¯×æë~/;Ò¿~~4__\ù½ö¯ŸÍ×_ıbú×Ïæëç?¸}üûÿü·õ÷ów¿¾ùıüúw¯“öïÂë`şö¿ÿÓW*ÿù³ã¿ÇçË¡\Òö#,_¿ıı×oÜÏ_¿ıûŸşíoŸüà®¯3àÃñG_ÇãÛzÙ~ìéøÓÁüi\Ÿß¿ÎÅëÏ­íóïÚçuß~¸rë/‰ûw§èüK—×{ı?¿'óÿê|ò›\Î‡küçü\şÚüµ¾®y>~•ı›¿øİ?ÿúO¿üñõ‡şæ7_ûëDı÷Ÿ¿şÌ»¯ßüò‡_~ıú›Ÿ¿~şÛ¯ßş¥şºøñ±Üí×Ey­ûuŸÊoëëñ/:•»U®XŞÏÕ/?-¯ÃÿùŸüå×_ÿÇ?üáëıü¶¯ßü—ßıó/ıô'_¼úFÕ-üjW»nÊëäÊ7Ô¯÷èõÂ”¯×W­}½~4_?Ÿ¹şõóãøz}÷Ú×ëÇñõú¶¯×æëçËØ¿~~_¯Oeûzı8¾^ßÎöõúq|½=¦ı7T?Ûï8Ÿ×ñçgóõÁíßQ?›ï¨OpÿúÙ~Çù(ï8?›ï¨ÏtÿúÙ|vwìöøï¾SÇŸÿ³µŸâ×\_/Àëá{…¼rşîşğOü‡ßıßşãEĞ±ùÅŞüˆğüG¼ŞİøîG,Ä¾çhËåÃqÿ÷Ï_Ç›ğËÿ^«y÷"¼ıAñTİÇÿ ’—ø?hùÆŠ%*¾ÿAëw–V«õãúÆÒ¼Ç}øğ“¾±¶’Ã-ñÃOŠßYœ%Û§¿µº°üHŸ@üÖê¼¾ÛëbGıIŞ}ky¯0XêO?ê;ë[·óÒøQ‹VßJ¹Sˆ_ö¿ßgçÇ–q|ŠfßøúÛhõ¯øõ·Ñê_ñëbÕ¿öw¼UÿÊßñ>İùÿOÑõº¯«à~,¡üwXÿë‘ßÂXùhÖ_>.cı¯»ë?ÿl_:J€¾şö±%Nõ'÷<)5Ñ8%±¼
cı¯Ï¯Çp¬¿üİ»Yùşİ¬¿|=šõ—Ï¹®Ieç“ßËoÁëú×ë_w¬ÿø8Ö.x¬ßcı+×¿&¬ÿøhÖ_?÷õ¿j³ş]Öüİfıç÷õŸ¿ß±ş­­ß/ËkåÛë(¿Štÿ÷ßsı™ëß¸şÌı÷Üÿûeÿ7®_¹şèeı‰ëßƒìÿ.ûïÛú·#2¤²ÿ)'İÿÈıÏ³õï²ÿ<ÿGjÖï¹ÿYÖŸeıºÿ»œÿMÖ¯çù°şôª®_/ßZNA\"ÖŸK•;Ş¿×ÙÌûç°şòf9óş¹s‹Æûwì`{ÿB‰öxÿìııÏÑÙ÷o?Óşşå!ñıËÎ¾{½íı[êûÑŞ¿úû.ëÏ¯ßÄëo=^÷ú-Ùõ§Ş¼ÿ¯Îÿ×Góş¿~Y‹İÿ#EûŸÎ©íÿë«vÿ_ßdßÿ×ç€óŸjì*I³Ùÿ|¢!}ÿÛç¶ÿé|nûş¿ş|ªëßCYùø¯¯ßÈ†ıwû¿x³ÿ¡ü²ìş¯fı¯`‰øw¬sì¿;Cß_şg³ÿKııŒø—ƒİ_ÏKÛÿå\ŸÙg×ÿŠk‹“øÿú5úûŸö?eìÿº˜ıÏˆÿå7Ìşç‚íı¯?Êî¿]ÿëófÏùáû¿9îÿn×_Îf²ûŸÏß×}ü¿‰ïŸÍúÇOïßññóûwÿ6»ÿí9µñÏæ?ıûMükëŸÅÿ×ûÇóïYÿQ%ó¿ŸoX;ÿ¹üËÆù<ÿáíü/õınçÿu—v{şk¾ÖÏÿv>=ıü¯\ÿëóùŞóïíú_¾ßÿyüÏû7ì²ù_}ÇşOÔØÿ¼`ÿ{ì’üï¼ßcÿ³ìÎÜÿ„ø_ÿ-æısˆé|Éîâÿq-Çú³¾ÿÙÿ×†ÛûŸÎ_–¹ÿöüç3”šûïåşÛøŸ]}¿Ûş»ïêş—m´ûïj¼oûïj>Ôöß±ú.ş×ùŸûšå?Ìÿ'Æärÿó×4ÿ•üg“üo[¹ÿgş`î¿=ÿı³¹ÿ-ÿ™ÄÿRÿí¨ÿìş—’Î™ûÿºsÉÜÿ¥†ƒzÿÄÿöqä?Ùî¹¯¶ş)õ]ùã’ÿ,‘÷ùïëóâxÿórÿ}ÂùßmıW_¸¶ş×?İ¾çÿõÕ`öÿõÕÅ¾¾å»}ıçzÛú_ïaæúÏóİÖêş¶õûZ·õ‡ó¯ëë÷g%ÿ#ïÿÆü×¯Xÿb÷ßGºÿ|&KãşmäŞñşoxÿıùÛëë÷õ~ŒüÇy®ç¿¼×ßë¿›úŸñãû¿1ş§iı¼8¦şaı»Ùø×ë•qÿ“İÿÌı—ú/Iı“äş§'ñ?•»Íÿ<×oñŸ×WÏ„wÄ?®ï*7Ë®?cı©Ş×‘ÿkşwæmıéÄ3Lü;Ÿ#ş­X­FîâæşçE÷ß›øß²£qş½­ÿêÇqş³¹Ö3ãü'›ÿ•×cáúî¿kxXõO®õÓˆ®åÿóú—óÏúgüË‚0ÿ‰<ÿ;÷—ü'Jıç°ÿ©Ö3&ÿY¹ÿç{`òÏóŸÛû7‹ÿ©”eöü/6ÿOç£:ògŞÿ–büg3õ9°æıïçsœÿßõÏQ.™ı_lük•™©Îû3ößãü¿~^;ÿ³ú?—;üéü—CiëßSÌûoÖ_^`‹8Äÿò“íû_Î'â_¦ıı5¿ï@ı[Áçşş{ÉÿÂùû¸‹ÿG1òŸõo¥íıo	_Ç¿ƒÉÿbİÖ¿É‚©ÿlßiëoïÁ¹şÂÉY‘ÿÖØÓ×Ÿ*ĞÖ_ñºşÕûøD)ÜûşÕçÀÄ?¾ÿrÿoî?ñ¿ˆøåü×èaï¿à¿QëÜÿV/™úçIıŸË.!ÿ±øG½ãşÛı/_µõO>#ºyÿ‰Á?vûş½Î¯—ü'Èû-şYî"ğGü³İ‡›ø"\vıGÖaò?[ÿ×šlÜÿÍæÿ5÷ßäåßgß¿¾ÿPÿ¾şı)qı)pı»İÿòûşUß‡»ø­öıÏ6ÿIgĞë?Ğ!³ÿûd,ãıOÿI…±„ı'şQáB³~/ûï%ÿ÷û¿9Ùÿgø¿ä¿;û»yÿ{Â×ï?Ö_»æşÛü?òı/÷Yğûş÷ÀÜà¿-Ÿ4÷_ğ¯üÿO?ã¿[°ÿÁ¾©6¼züüsÉØô¿Ríoøú§ç·ãıCıWòÉ¢àßëOç}»Çÿ×ë_‰ÿØú·a³~‡õoëß2Ö¿Éúwô?3ñ¿†&ÛüÇÿïŒõ'Á}Ë&ñÿli ÿ²ûŸ‰¸ğûı¬óÿùO@ş—ÿwÜÿöššûlÿ·¿§ãş'ô¿ı™ºİÇÿèg÷?®_üSúŸóü¿Ş÷~ÿ“Ôaçş/šÿcÿÛ÷ı_¤ÿÓóßIü?s¸‘ÿùeø×aş_0õOy‚MşSÈ&ş—ŸLü'!ş…†wükü;ÿXj}0ğ¯¬øòßõDcîëÿWp´÷›Ô?‰øÇñvÚûïìùO' ?îÿ‚ıO<ÿ½ß7ö?bÿ3×ŸRÍ÷FıÃú¯Fòûúñ?sıµ÷ß1ÿq¶ÿQ?üÇÙóŸë¿wä?^î?â_Ç³Æıß±ÿNø/õm2÷¿ãŸ7õ¿›õ?ÿÖèhâ?ñïdã_<°‰ÿÈê—Íş/Rÿo’ÿg¹ÿ^âŸ—şGèüŸÏñ¿”]¶ÿU£F»ÿõ…¯ü3W£ã¹şÕåzõ—¯î£ş[İz>‡çúW·´~÷±şÕm-œ¬0ÎŞJ«ÿ^Ÿ³Á¿ÖÂÜ3ùïêWà_«÷5?¨õßIË|ÿ#ïÿbñŸtá£şñ¶şËÀ¿;ÀôíıoñİÄ?äj³ø‡ôÿéÿ¹ÿxÿÊÏoïß<şâ+ûŸëÏÌÿ—ˆõÇŒõïÄÿøşåïÇúúıßÌıõÀ¸ÿÀÊ÷?‰ÿ™÷?3şg¼ÿŠÿûó…3ø§ÙÿÖõŸÿÏµ|øcş“ìû_"3ò‚ofıë÷çïëÿßÉŠìÿEòÿ"ûÿ‘ùàÿÀÿ®ø'óßXë[ƒÿ+ÿ-ş‘ÿğÿœÇï>ãÿ¤7õÍ/ùÿ:ÿ^â¿ş“*¿m¬Åıoña¼ÿ+øõcÖ¿<‰ÿ±\i»ÿIúÿŒAğßûş–¦ÿ±Êú×¿Hıã#×6®Ñø'ûôÿ3øO­K7î¿³ûŸqÿ;¬¿«ÅÿjÁcø_ÒÿóRÿı~Óÿü/¢ÿ“ëyïÿ†õ×şÉ}ÿ?“ÿ›î¿Ô?öü×7ùÅÿã	ÈşG’ı—ü×Kşëåü{9ÿZÿz©ÿ}{ÿæñgÿoµıïdóŸc‡7î?ë_Å¿ÿ\’ıß¹ÿKæş‡Àı÷û¿zîrÜÿ­½³øÿ]şWÂş£ÿ}Ëÿşï*ïÿ*ø—æ¿ä¿EyÿÚç7ü¯Yü?>»ÿAğÖü‹ÿ–ofÿw%ş£õú½öÿÚ¤şÛÁkŸMÿ·å¿7ıÖ?kÂşoû/üŸû/üÿü¯MêŸİqÿÑÿJQêßÈú¯ã¡ŸùÏ©ásşÿœÿ¾yî?ß¿ıŸí{ü‡•ûŸ¤ÿã¥ÿ¹Iÿ‘õ“ÿÎşâMüÏÖ_X‡–ÿ´€ÿœ÷óIoùo$ÿw?#tË·óÁjùïv&¦ÿ·ÿÏWş7ğ¿(üß$øOã“ş_èıyıÏüw—úŸúÄ÷oçùä?Úşgjúƒÿiş'ùÏ&ü/ùÏ*ùOû¿µó?Çÿ7âÿø×jğ¿†I÷o!şµ1şíìÿxéÿ¯Rÿ­À?scıYúÿŒÿ_ñ¯ç¿3şŸ/¾Å?ù_ÙêBİğÁ³õ_8£Ãà¿Ùşo©ç€yòÿø—öÿÃ‰gş»ş¯şµœïÃ}üß¾ÿìƒÿ)ü§ÄóŸpÿ`eú?™ûï%ş%©œ¼YŞ?ôÿËaDüOçyºÓÿyô¿[×yèÿÀ\€ÿ5ˆbàßË†ı·ëïı‰‡…û¿Úü?¯ïØ÷çÿõu÷/°ş-lŒgñ'ÿagşÿqÿk	oâŸàŸä?ìÿ„ÿ›pÿk%jû¢JÒÿL’ÿå'ü¿£,µëÏ\fıïÙÿ‰¬2ãfş{òÙÆùwÈS½ß¦şş‘ÈÿiH¥­$şúwÚÿßˆ€ÿ’¨¬OºÁ?ÈÄ? ÿkø„éÜèŸ¢ä¿»àßIê¿ø