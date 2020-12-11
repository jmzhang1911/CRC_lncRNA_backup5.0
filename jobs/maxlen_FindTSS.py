import gffutils
import argparse




def count_transmaxlen(geneclass,statdic,gffdb):
    """input a geneClass return a dict"""
    genetranslen = {}
    
    for i in gffdb.children(geneclass,featuretype='transcript'):
        trans_id = i['transcript_id'][0]
        tmplist = []
        translen = 0
        exonnum = 0
        
        for x in gffdb.children(i,featuretype='exon'):
            translen = translen +  (x.end - x.start + 1) 
            exonnum = exonnum + 1
            
        tmplist.append(exonnum)
        tmplist.append(translen)
        genetranslen[trans_id] = tmplist

    maxtrans = sorted(list(genetranslen.items()), key= lambda x:x[1][1], reverse=True)[0]
    statdic[geneclass.id] = list(maxtrans)
    return statdic


def count_stat(genelist,gffdb):
    statdic = {}
    for i in genelist:
        geneclass = gffdb[i]
        statdic = count_transmaxlen(geneclass,statdic,gffdb)
    
    print('#maxlength trans of genes\n#gene_id\ttrans_id\tmaxlength\texon_num')
    for k, v in statdic.items():
        geneid = k
        maxtransid = v[0]
        trans_len = v[1][1]
        trans_exon = v[1][0]
        print('{}\t{}\t{}\t{}'.format(geneid,maxtransid,trans_len,trans_exon))
        
    return 


def findtss(genelist,gffdb):
    print('#trans TSS of genes\n#gene_id\ttrans_id\tchr\tstart\tend\tstrand\tTSS')
    for i in genelist:
        geneclass = gffdb[i]
        gene_id = i
        for trans in gffdb.children(geneclass, featuretype='transcript'):
            chr_ = trans.seqid
            trans_id = trans.id
            start = trans.start
            end = trans.end
            
            if trans.strand == '-':    
                print(gene_id,'\t',trans_id ,chr_,'\t', start, end, trans.strand, trans.end)
            else:
                print(gene_id,'\t',trans_id ,chr_,'\t', start, end, trans.strand, trans.start)


def main():
    args = argparse.ArgumentParser('\nuse > /path/file to save results')
    args.add_argument('subcommand',type=str,choices=['maxtrans','findTSS'],help='subcommand')
    args.add_argument('-gffdb',type=str,help='GFFdb')
    args.add_argument('-genelist',type=str,help='a genelist')
    args = args.parse_args()
    
    gffdb = gffutils.FeatureDB(args.gffdb)
    with open(args.genelist,'r') as f:
        genelist = f.read().strip('\n').split('\n')

        if args.subcommand == 'maxtrans':
            count_stat(genelist,gffdb)
        else:
            findtss(genelist,gffdb)




if __name__ == '__main__':
    main()

