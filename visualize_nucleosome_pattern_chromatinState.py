import sys,os,argparse

from xplib import TableIO
from xplib import DBI
from xplib.Annotation import Bed
from xplib.Struct import binindex
from rpy2 import robjects
from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage,AnnotationBbox
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.gridspec as gridspec

# necessary since rpy 2.2.x, see http://stackoverflow.com/questions/2447454/converting-python-objects-for-rpy2
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()


def ParseArg():
    p=argparse.ArgumentParser( description='Visualization of epi- patterns and chromatin states in single nucleosome for populations given regions.', epilog="Library dependency: matplotlib, xplib, numpy, rpy2")
    p.add_argument('-i','--input',dest='input',type=str,required=True,help='input epi-pattern file, which is the output of [split_population.py]')
    p.add_argument('-H','--hist_n',dest='hist_n',default=7,type=int,help='number of histone modifications,default:7')
    p.add_argument('-c','--clu_n',dest='clu_n',default=3 ,type=int,help='number of clusters (populations),default:3')
    p.add_argument('-r','--region',dest='region',type=str,help="Genomic regions to be drawn, example: 'chr13:3220000-3350000'")
    p.add_argument('-o','--output',dest='output',type=str,help='output figure file, can be .pdf/eps/png/jpg...')
    p.add_argument("-g","--gene",type=str,dest="genetab",default='/home/yu68/xiaopeng-BAM2x/ensambl_mm9_gene.tab',  help="Known Gene Tab file (Download From UCSC genome browser, default: mouse mm9 )")
    p.add_argument('-b','--bed',action='store_true',help='If set, print out bed file to be uploaded to UCSC browser with chromatin state information')
    p.add_argument('-e','--emission',dest='emission',type=str,default='/home/yu68/split_population/chromHMM/output_model_E14/emissions_15.txt',help='emission matrix file from chromHMM output')
    p.add_argument('-s','--segment',dest='segment',type=str,default='/home/yu68/split_population/chromHMM/output_model_E14/E14_15_segments.bed',help='segment bed file from chromHMM output, for distribution of all chromatin states')
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        sys.exit(1)
    return p.parse_args()


def histone2state(histone_pattern,count):
    ''' Convert histone patterns into chromatin states based on chromHMM output trained from whole genome data '''
    robjects.globalenv["pattern"]=np.array(histone_pattern)
    #find log(p(y|cluster=c)) for each cluster (chromatin state)
    #robjects.r('apply(emission,1,function(x,y) prod(x^y)*prod((1-x)^(1-y)), y = c(1,1,0,0,0,1,1))')
    prob=robjects.r('log(as.matrix(emission)) %*% pattern + log(as.matrix(1-emission)) %*% (1-pattern)')  
    return np.argmax(np.log(count)+np.array(prob).T[0])+1

def addGeneToFig(gene,ax,start,name=0,bottom=0):

    '''
    add gene to figures
        start is the start of query region
    '''
    if name==1:
        ax.text((gene.start+start)/2,bottom+0.005,gene.id,fontsize=9)
    cds=gene.cds()

    utr5=gene.new_utr5()
    utr3=gene.utr3()
    if cds.stop!=cds.start:
        cds_exons=cds.Exons()
        for cds_exon in cds_exons:
            ax.bar(cds_exon.start,0.02,cds_exon.stop-cds_exon.start,color="blue",edgecolor="blue",alpha=1,bottom=bottom)
    if not utr3 is None:
        for utr3_exon in utr3.Exons():
            ax.bar(utr3_exon.start,0.01,utr3_exon.stop-utr3_exon.start,color="blue",edgecolor="blue",alpha=1,bottom=bottom+0.005)
    if not utr5 is None:
        for utr5_exon in utr5.Exons():
            ax.bar(utr5_exon.start,0.01,utr5_exon.stop-utr5_exon.start,color="blue",edgecolor="blue",alpha=1,bottom=bottom+0.005)
    for intron in gene.Introns():
        if intron.strand=="+":
            ax.arrow(intron.start,bottom+0.01,9*(intron.stop-intron.start)/10,0,color="black",edgecolor="black",alpha=0.7,width=0.002,head_width=0.004,head_length=(intron.stop-intron.start)/10,linestyle="dotted")
        else:
            ax.arrow(intron.stop,bottom+0.01,-9*(intron.stop-intron.start)/10,0,color="black",edgecolor="black",alpha=0.7,width=0.002,head_width=0.0044,head_length=(intron.stop-intron.start)/10)



def Main():
    args=ParseArg()

    hist_n=args.hist_n
    clu_n=args.clu_n
    File=args.input


    #read emission matrix and store in Rpy2
    print "#Reading emission matrix from"
    emission=args.emission
    print '\t'+emission
    robjects.r("emission=read.table('"+emission+"',header=T,sep='\t')")
    robjects.r("emission=emission[c(12,11,13,8,7,10,6,9,4,5,2,1,3,15,14),match(c('H3K4me3','H3K4me2','H3K4me1','H3K27me3','H3K36me3','H3K27ac','H2AZ'),colnames(emission))]")
    state_n=robjects.r("dim(emission)[1]")[0] # number of chromatin state
    
    color_state=['red','pink','purple','DarkOrange','Orange','Gold','yellow','DeepSkyBlue','ForestGreen','Green','Lime','GreenYellow','LightCyan','white','white']


    #Find overall distribution of all chromatin states
    print "Counting distribution of chromatin states..."
    chromHMM_segment = TableIO.parse(args.segment,'bed')
    #count represent overall probability distribution of all chromatin states
    count=np.zeros(state_n)
    num=0
    for segment in chromHMM_segment:
        num=num+1
        i=int(segment.id[1:])
        count[i-1]+=(segment.stop-segment.start)/200
        print 'Reading %d segments... [for distribution of chromatin states]'%(num),'\r',
    print


    ## read and index histone pattern data for single nucleosomes in all populations
    print "Indexing histone pattern data for single nucleosomes in all populations..."
    data=TableIO.parse(File,'metabed',header=True)


    ## generate bed file for chromatin states in nucleosomes to be uploaded in UCSC genome browser
    if args.bed:
        name=os.path.basename(File).split('.')[0]
        outbed=open(name+"_State_browser.bed",'w')
        print "## Start generate BED9 file for uploading..."
        print >>outbed,'track name="ChromatinState" description="'+name+'" visibility=2 itemRgb="On"'
        #print >>outbed,'chr\tstart\tend\t'+'\t'.join('P_%d'%(s+1) for s in range(clu_n))

        for n,i in enumerate(data):
            matrix=np.array(str(i).split('\t')[8:(8+hist_n*clu_n)],dtype="int").reshape(hist_n,clu_n,order="F")  # matrix of histone patterns, row: histone, column: population
            if n % 50000 == 0:
                print "\tWriting %dth nucleosomes into BED9 file,\r"%(n),
            line='\t'.join (str(f) for f in [i.chr,i.start,i.stop])
            for k in range(clu_n):
                state=histone2state(matrix.T[k],count)
                color_code=','.join (str(int(f)) for f in np.array(matplotlib.colors.colorConverter.to_rgb(color_state[state-1]))*255)
                print >>outbed,'\t'.join (str(f) for f in [i.chr,i.start,i.stop,'P_%d_%d'%(k+1,state),0,'.',i.start,i.stop,color_code])
                line=line+'\t%d'%(state)
            #print >>outbed,line
        outbed.close()
        sys.exit(1)


    # read region information
    region=args.region
    chro=region.split(":")[0]
    start=int(region.split(":")[1].split("-")[0])
    end=int(region.split(":")[1].split("-")[1])
    print "#Query region:["+chro+": %d-%d]"%(start,end)


    y_nucle=0.47 #location of nucleosome line

    
    ## query data in region
    dbi=binindex(data)
    query=dbi.query(Bed([chro,start,end]))

    ## initialize figure 
    fig=plt.figure(figsize=(10,6))

    ax = plt.subplot(111,frameon=False,yticks=[])
    ax.set_xlim(start-(end-start)/6,end)
    n=0
    print "##Start draw nucleosomes:"

    #################################################
    ## draw genes from y = y_nucle+0.04*(clu_n+1) 
    
    #### index the gene.tab file

    print "  ## drawing gene track ..."
    print "    ## Indexing gene.tab ..."
    gene_dbi=DBI.init(args.genetab,'genebed')


    print "    ## query regions from gene.tab"
    query_gene=gene_dbi.query(Bed([chro,start,end]))
    #### determine height of gene track    
    bottoms=[0 for i in range(100)]
    max_index=0
    for i in query_gene:
        index=0
        while(1):
            if i.start > bottoms[index]:
                bottoms[index]=i.stop
                if max_index<index: max_index=index
                break
            else:
                index+=1
    gene_track_number=max_index+1
    gene_track_height=0.03*gene_track_number+0.02
    ax.set_ylim(0.05,1+gene_track_height+0.01) 
    
    print "    ## start draw gene track"
    # add frame for gene track
    rect=matplotlib.patches.Rectangle((start,y_nucle+0.04),end-start, gene_track_height, edgecolor='black',fill=False)
    ax.add_patch(rect)
    
    bottoms=[0 for i in range(100)]
    for i in gene_dbi.query(Bed([chro,start,end])):
        index=0
        while(1):
            if i.start > bottoms[index]:
                addGeneToFig(i,ax,start,1,0.03*index+y_nucle+0.05)
                bottoms[index]=i.stop
                break
            index+=1

 
    ################################################# 
    
    top_heatmap_y = 0.71+gene_track_height # the y axis value for bottom of top heatmaps 

    print "##  Draw nucleosome tracks..."
    for i in query:
        n=n+1
        print "  Nucleosome %d\t at "%(n)+chro+": %d-%d"%(i.start,i.stop)
        matrix=np.array(str(i).split('\t')[8:(8+hist_n*clu_n)],dtype="int").reshape(hist_n,clu_n,order="F")  # matrix of histone patterns, row: histone, column: population
        prob=np.array(str(i).split('\t')[(8+hist_n*clu_n):],dtype=float)

        ax.plot([i.smt_pos,i.smt_pos],[y_nucle+0.03,y_nucle],color='r') #red nucleosome midpoint
        rect=matplotlib.patches.Rectangle((i.start,y_nucle), i.stop-i.start, 0.03, color='#EB70AA') #pink nucleosome region
        ax.add_patch(rect)

        for j in range(clu_n):
            state=histone2state(matrix.T[j],count)
            state_rect=matplotlib.patches.Rectangle((i.start,y_nucle+0.04*(j+1)+gene_track_height+0.01), i.stop-i.start, 0.03, color=color_state[state-1])
            ax.add_patch(state_rect)

    
        im = OffsetImage(matrix, interpolation='nearest',zoom=10/(1+gene_track_height+0.01),cmap=plt.cm.binary,alpha=0.5)
        

        if n<=9:
            xybox=((n+0.5)/10.0,top_heatmap_y)
            xy = [i.smt_pos,y_nucle+0.04*clu_n+0.03+gene_track_height+0.01]
            xytext=((n+0.7)/10.0,top_heatmap_y)
            c_style="bar,angle=180,fraction=-0.1"
        elif n<=18:
            xybox=((n-9+0.5)/10.0,0.2)
            xy = [i.smt_pos,y_nucle]
            xytext = ((n-9+0.7)/10.0,0.40)
            c_style="bar,angle=180,fraction=-0.1"
        else:
            print "WARN: nucleosome number larger than 18 in this region, only plot the pattern for first 18 nucleosomes"
            break

        ab = AnnotationBbox(im, xy,
                            xybox=xybox,
                            xycoords='data',
                            boxcoords=("axes fraction", "data"),
                            box_alignment=(0.,0.),
                            pad=0.1)
        ax.annotate("",xy,
                    xytext=xytext,
                    xycoords='data',
                    textcoords=("axes fraction", "data"),
                    arrowprops=dict(arrowstyle="->",connectionstyle=c_style))
                        #arrowprops=None)
    
        ax.add_artist(ab)
        
        # add mark for histone mark and regions with low confidence
        for i in range(hist_n):
            if prob[i]<0.6:
                xy_star=tuple(map(sum,zip(xybox,(0.065,0.03*(hist_n-1-i)-0.01))))
                ax.annotate("*",xy=xy_star,xycoords=("axes fraction", "data"),color='red')


    ax.annotate('Nucleosome:', xy=(start-(end-start)/6, y_nucle),  xycoords='data',size=12)
    ax.annotate('Epigenetic Pattern:', xy=(start-(end-start)/6, 0.23+top_heatmap_y),  xycoords='data',size=12)
    ax.annotate(chro, xy=(start-(end-start)/6, 0.1),  xycoords='data',size=12)

    name=open(File).readline().split('\t')[8:(8+hist_n)]
    for n,i in enumerate(name):
        ax.annotate(i.split("_")[0],xy=(start-(end-start)/8, top_heatmap_y+0.03*(hist_n-1-n)),xycoords='data',size=10)
        ax.annotate(i.split("_")[0],xy=(start-(end-start)/8, 0.2+0.03*(hist_n-1-n)),xycoords='data',size=10)

    # flame for nucleosome and chromatin state tracks
    rect=matplotlib.patches.Rectangle((start,y_nucle),end-start, 0.03, edgecolor='black',fill=False)
    ax.add_patch(rect)    
    for k in range(clu_n):
        rect=matplotlib.patches.Rectangle((start,y_nucle+0.04*(k+1)+gene_track_height+0.01),end-start, 0.03, edgecolor='grey',fill=False)
        ax.add_patch(rect)
        ax.annotate('Population%d'%(k+1),xy=(start-(end-start)/6, y_nucle+0.04*(k+1)+gene_track_height+0.01),xycoords='data',size=12)

    # chromatin state legend
    for s in range(state_n):
        dist=(end-start)*1.0/state_n 
        length=dist*0.75
        rect=matplotlib.patches.Rectangle((start+dist*s,0.1), length, 0.03, color=color_state[s])
        ax.add_patch(rect)
        ax.annotate(s+1,xy=(start+dist*s+length/3,0.075),xycoords='data',size=10) 
    ax.annotate("Chromatin states:",xy=(start,0.14),xycoords='data',size=12)      
    ax.add_patch(matplotlib.patches.Rectangle((start-length/6,0.07),end-start, 0.1, edgecolor='grey',fill=False))

    plt.title("Region: ["+chro+": %d-%d]"%(start,end),size=14)
    plt.savefig(args.output)

if __name__=="__main__":
    Main()

