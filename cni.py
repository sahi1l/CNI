#!/usr/bin/env python3
"""The following code implements the Cooke-Nieboer index, as shown in the paper by Sam A. Hill found in 
"Measure for characterizing heavy-tailed networks" in Phys. Rev. Research 3, 023257.
The class Stats implements Welford's online algorithm: 
"Note on a method for calculating corrected sums of squares and products". Technometrics. 4 (3): 419–420."""
import math
#========================================
class Stats:
    """Implements a class to keep track of running means and standard errors, etc."""
    count=0
    mean=0
    M2=0
    def init(self):
        self.count=0
        self.mean=0
        self.M2=0
    def add(self,newValue):
        self.count+=1
        delta=newValue-self.mean
        self.mean+=delta/self.count
        delta2=newValue-self.mean
        self.M2+=delta*delta2
    def variance(self):
        """Returns the sample variance"""
        return self.M2/(self.count-1)
    def stdev(self):
        """Returns the sample standard deviation"""
        return math.sqrt(self.variance())
    def sterr(self):
        return self.stdev()/math.sqrt(self.count)

#========================================

from numpy import sign
from random import choices

def cni(alldegs,maxerr=0.01):
    """Given a list of the degrees in a network, calculate the Cooke-Nieboer index using random sampling.
    Stop when the standard error is below maxerr.
    Multiple runs of this function on the same sequence should return values with a standard deviation of maxerr.
    """
    stats=Stats()
    stats.init()
    while True:
        degs=choices(alldegs,k=4)
        val=max(degs)+min(degs)-0.5*sum(degs) #Equivalent to X1+X4-(X2+X3)
        stats.add(sign(val))
        if(stats.count>20 and stats.sterr()<maxerr):
            #Stop if the error is below maxerr, but don't stop too soon.
            #The choice of 20 is arbitrary and can be changed as desired.
            return stats.mean

