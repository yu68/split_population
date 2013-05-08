from random import randrange,random
import numpy as np
import rpy2.robjects as robjects
import rpy2
from math import log

import sys,argparse

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from time import time
robjects.r('library(mgcv)')


#R functions for update S
'''
    Update_Sk <- function(o,p) {
                  y_p=p %*% S_m + 0.01*(0.5-p %*% S_m)  ## avoid 0,1 for predicted y into 0.01 and 0.99
                  prob_S=dbeta(o,beta*y_p^2*(1-y_p),beta*y_p*(1-y_p)^2)
                  num=which.max(prob_S)
                  return (S_o=S_m[,num]) }
'''
robjects.r("Update_Sk_EM <- function(o,p) {  y_p=p %*% S_m + 0.01*(0.5-p %*% S_m); prob_S=dbeta(o, beta*y_p^2*(1-y_p), beta*y_p*(1-y_p)^2); num=which.max(prob_S); return (S_o=S_m[,num]); }")
robjects.r("Update_Sk <- function(o,p) { y_p=p %*% S_m + 0.01*(0.5-p %*% S_m); prob_S=dbeta(o, beta*y_p^2*(1-y_p), beta*y_p*(1-y_p)^2); print (prob_S); print(o); print (y_p); num=sample(1:(2^clusN),1,replace=TRUE,prob=prob_S); return (S_o=S_m[,num]); }")


def SimulateRealData():

    #simulating real data:
    S_real=robjects.r("S_real=matrix(floor(runif(nu_n*hist_n*clu_n,0,2)),nrow=clu_n)")
    robjects.r("p_real=runif(clu_n,5*(clu_n-1)/95,1)")
    p_real=robjects.r("p_real=p_real/sum(p_real)")
    robjects.r("y_p=p_real %*% S_real + 0.01*(0.5-p_real %*% S_real)")
    y=robjects.r("y=as.vector(rbeta(nu_n*hist_n,1000*y_p^2*(1-y_p), 1000*y_p*(1-y_p)^2))")
    return [S_real,p_real,y]



def update_S(determ):
    if determ:
    #robjects.globalenv["p"] = p
        robjects.r("for (i in 1:dim(S)[1]) { S[i,]=Update_Sk_EM(y[i],p); }")
    else:
        robjects.r("for (i in 1:dim(S)[1]) { S[i,]=Update_Sk(y[i],p); }")



def restricted_LS(S,determ):
    #update S
    robjects.globalenv["S"] = S
    update_S(determ)
    S2=robjects.r("S")

    p_ini=robjects.r("p_i=rep(1,clusN)/clusN")

    M = robjects.r("M=list(X=S,p=p_i,Ain=diag(clusN),bin=rep(0,clusN),C=matrix(1,1,clusN),y=y,w=y*0+1)")

    p=robjects.r("p=pcls(M)")

    Variance=robjects.r("function(x) {y_p = as.vector(S %*% x); return (((y-y_p) %*% (y-y_p))/(hist_n*nu_n))}")

    #robjects.globalenv["beta"] = 1.0 /(8*Variance(p)[0])
    Log_L=robjects.r("function(p) {y_r = y - as.vector(S %*% p); prob_all=dbeta((y_r + 1)/2, beta, beta); return (sum(log(prob_all))) }")
    return S2,p,np.sqrt(Variance(p))[0][0],Log_L(p)[0]


def count_matrix(S):
    ''' count the occurancy of each S pattern in predicted S matrix '''
    S=np.array(S)
    count=np.zeros(np.power(2,len(S[0])))
    for i in range(len(S)):
        number=0
        for j in range(len(S[i])):

            number+=S[i][j]*np.power(2,j)
        count[number]+=1
    print count
    return



def Run1Prediction(y,clusN,iteration,determ=False,verbose=False,hist_n=1,nucleosome_n=1):
    t1=time()
    robjects.globalenv["clusN"]=clusN   # clusN is not the real cluster number

    #choose S for each loc and each marker
    robjects.r("S_m=NULL;for (i in 1:clusN) { S_i=floor((0:(2^clusN-1))/2^(i-1)) %%2; S_m=rbind(S_m,S_i);};")


    #initial S
    S1=robjects.r("S=matrix(floor(runif(nu_n*hist_n*clusN,0,2)),nrow=nu_n*hist_n)")
    #initial p
    robjects.r("p=runif(clusN,5*(clusN-1)/95,1)")
    p=robjects.r("p=p/sum(p)")
    print "Initial p: ", robjects.r('p')
        

    diff1=1E8
    diff_min=1E8
    diff_list=[]
    log_l1=0
    log_l_max=0
    index_min=0
    p_min=p
    for i in range(0,iteration):
        if i>1:
            if verbose:
                print "iteration: %d\tdiff_min: %.6f\tdiff_now:%.6f\tlog_likelihood:%.1f"%(i,diff_min,diff1,log_l1)
            else:
                print "iteration: %d\tdiff_min: %.6f\tdiff_now:%.6f\tlog_likelihood:%.1f\r"%(i,diff_min,diff1,log_l1),
        if i<int((iteration)*0.1):
            determin=False
        if i==int((iteration)*0.1):  #start from here, using EM for the best sampling
            S1=S_min
            determin=determ
            robjects.globalenv["p"]=p_min
        try:
            S2,p2,diff2,log_l2=restricted_LS(S1,determin)
        except rpy2.rinterface.RRuntimeError:
            # sampled S is rank deficient give a large diff and didn't change S
            print "Rank deficiency of updated S, change one"
            diff_list.append(np.sqrt(len(y))*0.5)

            robjects.r("p=runif(clusN,5*(clusN-1)/95,1)")
            robjects.r("p=p/sum(p)")
            continue
        if np.exp(np.log(i+1)*(diff1-diff2))>np.random.uniform(0,1,1)[0]:
            S1=S2
            p1=p2
            diff1=diff2
            log_l1=log_l2
        if diff2<diff_min:
            S_min=S2
            p_min=p2
            diff_min=diff2
            log_l_max=log_l2
            index_min=i
        diff_list.append(diff1)
        if diff_min<0.005:
            break

    print
    print "Best diff:",diff_min
    print "time required: ",time()-t1
    print "Max log_likelihood: ",log_l_max
    print "iteration for best:",index_min

    order=[i[0] for i in sorted(enumerate(p_min), key=lambda x:x[1])]
    S_min=np.matrix(S_min)
    S_min=S_min[:,order]
    p_min=np.array(p_min)
    p_min=p_min[order]

    if verbose:
        if determ:
            pdf=PdfPages("Diff_changes_determ_clus%d.pdf"%(clusN))
        else:
            pdf=PdfPages("Diff_changes_sampling_clus%d.pdf"%(clusN))
        fig=plt.figure()
        plt.plot(range(1,(iteration+1)),diff_list,color="blue",lw=0.5)
        plt.title("cluster:%d,histM:%d,nucleosome:%d"%(clusN,hist_n,nucleosome_n))
        pdf.savefig(fig)
        pdf.close()



    print "Predicted P:"
    print p_min

    count_matrix(S_min)
    return [p_min,S_min,index_min,log_l_max]


