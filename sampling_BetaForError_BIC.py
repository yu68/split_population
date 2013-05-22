from random import randrange,random
import numpy as np
import rpy2.robjects as robjects
import rpy2
from math import log

import sys,argparse
from functions import *


import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from time import time
robjects.r('library(mgcv)')
# need R library mgcv


global cluster_n,hist_n,nucleosome_n,beta,iteration

def ParseArg():
    p=argparse.ArgumentParser( description = "Simulate data and predict it")
    p.add_argument("-c","--cluster",dest="cluster_n",type=int,default=4,help="cluster number,default:4")
    p.add_argument("-m","--histM",dest="hist_n",type=int,default=2,help="histone marker number,default:2")
    p.add_argument("-n","--nucleosome",dest="nucleosome_n",type=int,default=10,help="nucleosome number,default:10")
    p.add_argument("-b","--beta",dest="beta",type=float,default=500,help="beta distribution parameter beta,default:500")
    p.add_argument("-o","--output",dest="output",type=str,help="output file")
    p.add_argument("-i","--iteration",dest="iteration",type=int,default=5000,help="iteration number,default:5000")
    p.add_argument("-d","--deterministic",action='store_true',help="Choose the S with highest P instead of sampling")
    p.add_argument("-r","--replicates",dest="replicates",type=int,default=100,help="replicate number,default:100")
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        exit(0)
    return p.parse_args()

'''

def SimulateRealData():

    #simulating real data:
    S_real=robjects.r("S_real=matrix(floor(runif(nu_n*hist_n*clu_n,0,2)),nrow=nu_n*hist_n)")
    robjects.r("p_real=runif(clu_n,5*(clu_n-1)/95,1)")
    p_real=robjects.r("p_real=p_real/sum(p_real)")
    print "real p: "
    print p_real
    y=robjects.r("y=as.vector(S_real %*% p_real+2*rbeta(hist_n*nu_n,10000,10000)-1)")
    #print "Obervation:"
    #print (robjects.r("y"))
    return [S_real,p_real,y]
'''

args=ParseArg()

cluster_n=args.cluster_n
hist_n=args.hist_n
nucleosome_n=args.nucleosome_n
beta=args.beta
iteration=args.iteration

robjects.globalenv["clu_n"] = cluster_n
robjects.globalenv["hist_n"] = hist_n
robjects.globalenv["nu_n"] = nucleosome_n
robjects.globalenv["beta"] = beta


matrix_accuracies=[]
p_accuracies=[]

parameters="Clut-%d_Hist-%d_Nucleo-%d_beta-%d"%(cluster_n,hist_n,nucleosome_n,beta)
output=open(parameters+args.output,"w")

if args.deterministic:
    pdf=PdfPages("Clut-%d_sampling_simulation_BIC_determ.pdf"%(cluster_n))
else:
    pdf=PdfPages("Clut-%d_sampling_simulation_BIC_sampling.pdf"%(cluster_n))
fig=plt.figure()

Correct_n=0
y_legend=0
for i in range(0,args.replicates):
    print 
    print "Run number:",i+1
    print >>output,"Run number:",i+1
    [S_real,p_real,y]=SimulateRealData()
    BIC_pool=[]
    print >>output,"cluster\tbestIter\tlog_l\tBIC"
    for clusN in range(2,8):
        [p_pred,S_pred,index,log_l]=Run1Prediction(y,clusN,iteration,args.deterministic)
        BIC=-2*log_l+clusN*log(hist_n*nucleosome_n)
        BIC_pool.append(BIC)
        print "BIC is: ",BIC
        print >>output,"%d\t%d\t%.1f\t%.1f"%(clusN,index,log_l,BIC)

    # Judge if BIC give right cluster number
    if min(BIC_pool)==BIC_pool[cluster_n-2]:
        color='blue'
        Correct_n+=1
    else:
        color='red'
    if max(BIC_pool)>y_legend:
        y_legend=max(BIC_pool)  # Get y-axis coordinate of legend
    plt.plot(range(2,8),BIC_pool,color=color,lw=0.5)

plt.xlabel('cluster number')
plt.ylabel('BIC')
plt.text(5,y_legend-10,"%d / %d are CORRECT"%(Correct_n,args.replicates))
plt.title("cluster:%d,histM:%d,nucleosome:%d,beta:%d"%(cluster_n,hist_n,nucleosome_n,beta))
pdf.savefig(fig)
pdf.close()
print "%d / %d are CORRECT"%(Correct_n,args.replicates)    
output.close()
