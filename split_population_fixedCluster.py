from random import randrange,random
import numpy as np
import rpy2.robjects as robjects
import rpy2
from rpy2.robjects.numpy2ri import numpy2ri

import sys,argparse
from functions import *

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
#from matplotlib.backends.backend_pdf import PdfPages

import pp

from time import time
# need R library mgcv


global cluster_n,hist_n,nucleosome_n,beta,iteration

def ParseArg():
    p=argparse.ArgumentParser( description = "Split ES cells pouplation based on Epi-intensities in Single nucleosome level.",epilog="Library dependency: random, rpy2, matplotlib, pp")
    p.add_argument("input",type=str,help="Input file for relative Epi-intensities in Single nucleosome level")
    p.add_argument("-c","--cluster",dest="cluster_n",type=int,default=4,help="cluster number,default:4")
    #p.add_argument("-m","--histM",dest="hist_n",type=int,default=2,help="histone marker number,default:2")
    #p.add_argument("-n","--nucleosome",dest="nucleosome_n",type=int,default=10,help="nucleosome number,default:10")
    p.add_argument("-b","--beta",dest="beta",type=float,default=500,help="beta distribution parameter beta,default:500")
    p.add_argument("-o","--output",dest="output",type=str,help="output file")
    p.add_argument("-i","--iteration",dest="iteration",type=int,default=5000,help="iteration number,default:5000")
    p.add_argument("-d","--deterministic",action='store_true',help="Choose the S with highest P (using EM) instead of sampling")
    p.add_argument("-r","--replicates",dest="replicates",type=int,default=1,help="replicate number,default:1")
    p.add_argument("-v","--verbose",action='store_true',help="print verbose information and difference change plot")
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        exit(0)
    return p.parse_args()

def ReadInputData(input_file):

    #Input real data:
    print "Reading input data for Epi-intensities..."
    robjects.r("input=read.table('"+input_file+"',sep='\t',header=T)")
    robjects.r("input=input[rowSums(input[,9:dim(input)[2]]>0.1)>0,]")
    hist_n=int(robjects.r("hist_n=dim(input)[2]-8")[0])
    print "Number of histone madifications: ",hist_n
    nucleosome_n=int(robjects.r("nu_n=dim(input)[1]")[0])
    print "Number of nucleosomes: ",nucleosome_n
    y=robjects.r("y=as.vector(as.matrix(input[,9:dim(input)[2]]))")
    return [hist_n,nucleosome_n,y]



args=ParseArg()

cluster_n=args.cluster_n
beta=args.beta
iteration=args.iteration

robjects.globalenv["clu_n"] = cluster_n
robjects.globalenv["beta"] = beta



#R functions for update S
'''
    Update_Sk <- function(o,p) {
                  prob_S=dbeta((o-(p %*% S_m)-1)/2, beta, beta)
                  num=sample(1:(2^clu_n),1,replace=TRUE,prob=prob_S)
                  return (S_o=S_m[,num]) }
'''
'''
robjects.r("Update_Sk_EM <- function(o,p) {  prob_S=dbeta((o-(p %*% S_m)+1)/2, beta, beta); num=which.max(prob_S); return (S_o=S_m[,num]); }")
robjects.r("Update_Sk <- function(o,p) {  prob_S=dbeta((o-(p %*% S_m)+1)/2, beta, beta); num=sample(1:(2^clu_n),1,replace=TRUE,prob=prob_S); return (S_o=S_m[,num]); }")
'''

[hist_n,nucleosome_n,y]=ReadInputData(args.input)
parameters="Hist-%d_Nucleo-%d_beta-%d"%(hist_n,nucleosome_n,beta)
output=parameters+args.output
for i in range(0,args.replicates):
    print 
    print "Run number:",i+1
    
    '''
    ################################
    ## parallel computing attempt ##
    ppservers=()
    t1=time()
    # Creates jobserver with nprocessor workers
    job_server = pp.Server(args.nprocessor,ppservers=ppservers)
    
    job = job_server.submit(Run1Prediction,(y,cluster_n,iteration,args.deterministic,args.verbose,hist_n,nucleosome_n,beta),\
                           (),("from functions import *","import numpy as np",'rpy2','math','time',\
                           'from rpy2 import robjects', 'from math import log',))
    [p_min,S_min,index,log_l] = job()
    print "time required: ",time()-t1
    ################################
    '''
    t1=time()
    [p_min,S_min,index,log_l]=Run1Prediction(y,cluster_n,iteration,args.deterministic,args.verbose,hist_n,nucleosome_n,beta)
    print "time required: %.4fmin"%((time()-t1)/60)
    BIC=-2*log_l+cluster_n*log(hist_n*nucleosome_n)
    print "BIC is: ",BIC

    output='Clut-%d_'%(cluster_n)+parameters+args.output
    ### output resulting matrix
    S_min=np.array(S_min, dtype="int")
    S_min=numpy2ri(S_min)        
    robjects.r.assign("predict_S",S_min)
    robjects.r("predict_S <- matrix(predict_S,ncol=clu_n*hist_n)")
    robjects.r('prob.S <- matrix(prob.S,ncol=hist_n)')
    robjects.r("output_matrix <- cbind(input[,1:8],predict_S,prob.S)")

     ## column number
    robjects.r("name=NULL")
    robjects.r("for (i in 1:clu_n) { name=c(name,paste(colnames(input)[9:(8+hist_n)],'_P',i,sep='')); }")
    robjects.r("colnames(output_matrix)[9:(8+hist_n*clu_n)]=name;  colnames(output_matrix)[3]='stop'; ")
    robjects.r("colnames(output_matrix)[(9+hist_n*clu_n):dim(output_matrix)[2]]=colnames(input)[9:(8+hist_n)];")
     ## print output
    robjects.r("write.table(output_matrix,'"+output+"',sep='\t',quote=F,col.names=T,row.names=F)")

    p_line="# population proportion: ["+",".join(str(f) for f in p_min)+"]."
     ##add populiation proportion at the top
    os.system("sed -i '1i"+p_line+"' "+output)
