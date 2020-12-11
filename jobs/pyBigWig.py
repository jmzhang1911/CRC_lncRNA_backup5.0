import pyBigWig
from pathlib import Path
import argparse




def bw2bed(bw,outname):
  of = open(outname, "w")
  for chrom, len in bw.chroms().items():
  intervals = bw.intervals(chrom)
    for interval in intervals:
      of.write('{}\t\{}\t{}\t{}\n'.format(chrom,interval[0],interval[1],interval[2]))
  

def getbw(bwlist):
  for i in bwlist:
    p = Path(i)
    outname = p.stem + '.bed'
    bw = pyBigWig.open(p)
    getbw(bw,outname)




if __name__ == '__main__':
  args = argparse.ArgumentParser('BigWig to bed')
  args.add_argument('filelist',type=str,help='filelist of file.bw')
  args = args.parse_args()
  
  with open(bwlist, 'r') as f:
    bwlist = f.read().strip('\n').split('\n')
    
  getbw(args.filelist)




