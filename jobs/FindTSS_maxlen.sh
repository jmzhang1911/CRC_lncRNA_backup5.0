#!/bin/bash
#pwd=getwd()
python jobs/maxlen_FindTSS.py maxtrans -gffdb ./data/ensemble_noncode_gtf.db \
	-genelist ./outcomes/candidate/dy_lnclist.txt > ./outcomes/candidate/maxtrans.txt 

python jobs/maxlen_FindTSS.py findTSS -gffdb ./data/ensemble_noncode_gtf.db \
	-genelist ./outcomes/candidate/dy_lnclist.txt > ./outcomes/candidate/trans_TSS.txt 