#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       Copyright 2009 Guillaume Lanquepin <guillaume@himolde.no>
#       version 0.1

import numpy as np
from solvers import thomas,ipopt
from tools import import_data
from plots import show_data

class PCLSP:
    """
    """

    def __init__(self, data_file="data_test.csv", theta=0.3, eps=0.1, cycle=500, verbose=True):
        #problem paramters
        self.time_hor, \
        self.nb_obj, \
        self.alpha, \
        self.beta, \
        self.cost_setup, \
        self.cost_stor, \
        self.cost_prod, \
        self.cons_prod, \
        self.constraint = import_data(data_file)
        #problem variables
        self.coef = np.zeros(((self.nb_obj+1)*self.time_hor), float) # lagrangian multiplier
        self.production = np.zeros((self.nb_obj, self.time_hor), float)
        self.price = np.zeros((self.nb_obj, self.time_hor), float)
        self.setup = np.zeros((self.nb_obj, self.time_hor), float)
        self.storage = np.zeros((self.nb_obj, self.time_hor), float)
        #algorithm parameters
        self.eps = eps #min difference between upper and lower bound
        self.cycle = cycle #max number of iteration (stop condition)
        self.smooth_param = theta
        self.verbose = verbose
        self.opt_upper = []#vector fo optimun return by Upper
        self.opt_lower = []#vector fo optimun return by Lower
        if self.verbose:
            show_data(self)

    def get_data(self):
        return self.alpha, self.beta, self.cost_prod, self.cost_stor, self.cons_prod,
        self.cost_setup, self.coef, self.time_hor, self.nb_obj, self.verbose
        
    def _critere(self):
        """
        
        """
        objective = 0
        for j in xrange(self.nb_obj):
            for t in xrange(self.time_hor):
                objective += (self.alpha[j,t]-self.beta[j,t]*self.price[j,t])*self.price[j,t]\
                - self.cost_setup[j,t]*self.setup[j,t] \
                - self.cost_stor[j,t]*self.storage[j,t] \
                - self.cost_prod[j,t]*self.production[j,t]
        return objective
        
    def _critere_ipopt(self,initial):
        """
            
        """
        objective = initial
        for j in xrange(self.nb_obj):
            for t in xrange(self.time_hor):
                objective -= self.cost_setup[j,t]*self.setup[j,t]
        return objective
        
    def solve(self):
        """
        
        """
        diff = self.eps+ 1.
        count = 0
        while ( (diff > self.eps) & (count < self.cycle) ):
            previous_lambda = self.coef
            self.price, self.setup = thomas(self.get_data ())
            upper = self._critere()
            #try:
            lower, self.coef = ipopt(self.get_data ())
            lower = self._critere_ipopt(lower)
            #except:
            #    print "Ipopt solver unaviable: run heuristic..."
            #    self.price, self.coef = heuristic(self)
            #    lower = self._critere()
            self.opt_upper.append(upper)
            self.opt_lower.append(lower)
            # update lagrangian coefficients
            for i in xrange(len(self.coef)):
                self.coef[i] = self.smooth_param*previous_lambda[i]\
                - (1-self.smooth_param)*self.coef[i] 
            diff = upper-lower
            count += 1
        if self.verbose:
            print "Number of iteration: "+str(count)
            print "Difference between upper and lower bound: "+str(diff)
            print "Last critere value: "+str(lower)
            show_data(self)
        return self._critere ()

if __name__ == '__main__':
    TEST = PCLSP()
    TEST.solve()

