# author: Pengfei Yu

import sys,os

class nucleosome():
    '''
    object for single nucleosome with all informations:
        chro, start, end, center, smt_value, fuzziness_score, fuzziness_pval
        
    '''
    def __init__(self,**dict):
        attr=dict["attr"]
        for i in range(len(attr)):
            attr[i]=attr[i].strip()
        self.attr=dict["attr"]
        value=dict["value"]
        if len(value)!=len(attr):
            print >>sys.stderr,"length of value and attr are not the same"
            sys.exit(0)
        for i in range(len(attr)):
            self.__setattr__(attr[i],value[i])
        if not self.__dict__.has_key('chr'):
            if self.__dict__.has_key('chromosome'):
                self.__setattr__('chr',self.chromosome)
            if self.__dict__.has_key('chrom'):
                self.__setattr__('chr',self.chrom)
    def epiCounts(counts):
        self.__setattr__("epiCounts",counts)
    def epiMarks(names):
        self.__setattr__("epiMarks",names)
    def offsets(offsets):
        self.__setattr__("offsets",offsets)
        

        
