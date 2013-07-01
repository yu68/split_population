import sys,argparse,os
from xplib import TableIO
from xplib import DBI
from xplib.Annotation import Bed
import numpy as np

def ParseArg():
    p=argparse.ArgumentParser(description="get intensity of epi-modification on single nucleosome",epilog="library dependency: xplib")
    p.add_argument("-N","--Nucleosome",dest="nucleosome",type=str,help="xls file containing nucleosome location information from Danpos output")
    p.add_argument("-b","--beds",nargs='+',dest="beds",type=str,help="bed/bam files for epigenetic data")
    p.add_argument('-f','--fmt',type=str,default='bed',help='format: bed/bam,default:bed')
    p.add_argument("-l","--length",dest="len",type=int,default=200,help="average length of ChIP-seq fragment,default:200")
    p.add_argument("-n","--name",nargs='+',dest='name',type=str,help='name of each bed sample (to be wrote on the header)')
    p.add_argument("-w","--weightP",nargs='+',dest='weightP',type=int,default=[75,125],help='parameters to calculate the weight for each read,[half_len of core nucleosome region and half_len of whole regions],default: [75,125]')
    p.add_argument("-o","--output",dest="output",type=str,help="output file name (can be .txt)")
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        exit(0)
    return p.parse_args()

def find_center(read,half_len,fmt):

    # read.id in strand for bed with for columns
    if fmt=='bed':
        strand=read.id
    elif fmt=='bam':
        strand=read.strand
    if strand=="+":
        center=read.start+half_len
    elif strand=="-":
        center=read.stop-half_len
    else:
        print >>sys.stderr,read.id,"No strand"
    return center


'''
      mi
    -----
    |    \
    |     \
    |      \
    |       \
    --------
  Mid  ma

'''

global mi,ma,args
args=ParseArg()
mi=args.weightP[0]
ma=args.weightP[1]


def Main():
    args=ParseArg()

    #store bed files with indexing and count information:
    bed={}

    print >>sys.stderr,"Starting index bed files:"
    for i in range(len(args.beds)):
        temp_name=args.name[i]
        print >>sys.stderr,"  #Indexing for bed file of",temp_name,"\r",
        bed[temp_name]=DBI.init(args.beds[i],args.fmt)

    half_len=int(args.len)
    print >>sys.stderr
    print >>sys.stderr,"Reading nucleosome peak xls file from Danpos."
    nucleosomes=TableIO.parse(args.nucleosome,'metabed',header=True)

    print >>sys.stderr,"Start Counting..."
    count_matrix=[]


    out=open(args.output,"w")
    line_head=open(args.nucleosome,'r').readline().strip()
    line_head=line_head+"\t"+"\t".join(str(f) for f in args.name)
    print >>out,line_head
    for i in nucleosomes:
        chrom=i.chr
      
        if chrom == 'chrY' or chrom == 'chrX' or chrom == 'chrM':
            continue
        center=int(i.start+i.end)/2
        count=np.zeros(len(args.beds),dtype="float")
        line=str(i)
        for k,name in enumerate(bed.keys()):
            if args.fmt=='bam':
                query=bed[name].query(Bed([chrom,center-ma-(half_len-75),center+ma+(half_len-75)]),method='fetch')
            else:
                query=bed[name].query(Bed([chrom,center-ma-(half_len-75),center+ma+(half_len-75)]))
            for j in query:
                j_center=find_center(j,half_len,args.fmt)
                weight = max(min(1,(ma-abs(j_center-center))*1.0/(ma-mi)),0)
                count[k]+=weight
        line = line + "\t" + "\t".join(str(f) for f in count)
        print >>out,line
        count_matrix.append(count)

if __name__=="__main__":
    Main() 


        
                        

'''
TableIO.parse("Pn_E14_merge_mm9.sorted.Fnor.ajClonal.smooth.peaks.xls",'metabed',header=True)
'''



