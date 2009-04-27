#!/usr/bin/env python
#
#       Copyright 2009 Guillaume Lanquepin <guillaume@himolde.no>
#       version 0.1

import numpy as np
import tools as tl

class ClspHeuristic:

    def __init__(self,dataFile="data_test.csv", _eps=0.1, _cycle=500,verbose=True):
        #problem paramters
        self.t, \
        self.J, \
        self.alpha, \
        self.beta, \
        self.costSetup, \
        self.costStor, \
        self.costProd, \
        self.consProd, \
        self.constraint = tl.importData(dataFile)
        #problem variables
        self.coef = np.zeros(self.t) # lagrangian multiplier
        self.production = np.zeros((self.J,self.t))
        self.price = np.zeros((self.J,self.t))
        self.setup = np.zeros((self.J,self.t))
        self.storage = np.zeros((self.J,self.t))
        #algorithm parameters
        self.eps = _eps #min difference between upper and lower bound (stop condition)
        self.cycle = _cycle #max number of iteration (stop condition)
        self.optUpper = np.zeros(self.t) #vector fo optimun return by Upper
        self.optLower = np.zeros(self.t) #vector fo optimun return by Lower
        if verbose:
            self.__verbose (True)

    def __verbose(self,start):
        if start:
            print "### Data"
            print "Period:"+str(self.t)+"\tProduct:"+str(self.J)
            print "Alpha:",
            tl.printtab(self.alpha)
            print "Beta:",
            tl.printtab(self.beta)
            print "Setup Cost:",
            tl.printtab(self.costSetup)
            print "Storage cost:",
            tl.printtab(self.costStor)
            print "Prodcution cost:",
            tl.printtab(self.costProd)
            print "Product consumption:",
            tl.printtab(self.constraint)
            print "Production constraint",
            tl.printtab(self.constraint)

    def __upper(self):
        for t in xrange(self.t):
            for t0 in xrange(t):
                self.price[t][t0] = (self.alpha[t] + \
                (self.costProd[t0] + sum(self.costStor[t0:(t-1)])+self.coef[t0]*self.consProd[t0]) \
                * self.beta[t]) / (2 * self.beta[t] )
                self.production[t][t0] = self.alpha[t] - self.beta[t] * self.price[t][t0];
        return 0
        
    def __lower(self):
        return 0
    def main(self):

        return 0

if __name__ == '__main__':
    test = ClspHeuristic()

