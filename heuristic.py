#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       Copyright 2009 Guillaume Lanquepin <guillaume@himolde.no>
#       version 0.1

import numpy as np
import time
from solvers import thomas, ipopt
from tools import import_data
from plots import show_data, graphic
from math import fabs

class HEURISTIC:
    """
    This class containts everything for compute
    K. Haugen and A. Olstad heuristic for PCLSP
    
    See the profit maximizing capacitated lot-size problem (7 october 2007)
    """

    def __init__(self, theta=0.3, eps=0.1, cycle=500, verbose=1):
        #problem paramters
        self.time_hor = []
        self.nb_obj = []
        self.alpha = []
        self.beta = []
        self.cost_setup = []
        self.cost_stor = []
        self.cost_prod = []
        self.cons_prod = []
        self.constraint = []
        #problem variables
        self.coef = []
        self.production = []
        self.price = []
        self.setup = []
        self.storage = []
        #algorithm parameters
        self.eps = eps #min difference between upper and lower bound
        self.cycle = cycle #max number of iteration (stop condition)
        self.smooth_param = theta
        self.verbose = verbose # 0 not information, 1 input and result, 2 all
        self.opt_upper = [] #vector fo optimun return by Upper
        self.opt_lower = [] #vector fo optimun return by Lower

    def import_data(self, data_file="data_test.csv"):
        """
        Function to import data from csv file
        
        import_data(self,data_file="data_test.csv")
        """
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
        if self.verbose > 0:
            show_data(self)
    
    def set_data(self, alpha, beta, cost_prod, cost_stor,
        cons_prod, constraint, cost_setup, time_hor, nb_obj):
        """
        Function to define data
        
        set_data(self,alpha, beta, cost_prod,cost_stor, cons_prod,constraint
        cost_setup, coef, time_hor, nb_obj, verbose)
        """
        self.time_hor = time_hor
        self.nb_obj = nb_obj
        self.alpha = alpha
        self.beta = beta
        self.cost_setup = cost_setup
        self.cost_stor = cost_stor
        self.cost_prod = cost_prod
        self.cons_prod = cons_prod
        self.constraint = constraint
        #problem variables
        self.coef = np.zeros(((self.nb_obj+1)*self.time_hor), float) # lagrangian multiplier
        self.production = np.zeros((self.nb_obj, self.time_hor), float)
        self.price = np.zeros((self.nb_obj, self.time_hor), float)
        self.setup = np.zeros((self.nb_obj, self.time_hor), float)
        self.storage = np.zeros((self.nb_obj, self.time_hor), float)
        if self.verbose > 0:
            show_data(self)
        
    def get_data(self):
        """
        Function to get back all data in an tuple
        
        get_data(self)
        """
        return [self.alpha, self.beta, self.cost_prod, 
            self.cost_stor, self.cons_prod, 
            self.cost_setup, self.coef, self.time_hor, 
            self.nb_obj, self.verbose]
        
    def _critere(self):
        """
        Private function who compute critere
        
        _critere(self)
        """
        objective = 0
        for j in xrange(self.nb_obj):
            for t in xrange(self.time_hor):
                objective += (self.alpha[j, t]-self.beta[j, t]*self.price[j, t])\
                *self.price[j, t]\
                - self.cost_setup[j, t]*self.setup[j, t]\
                - self.cost_stor[j, t]*self.storage[j, t]\
                - self.cost_prod[j, t]*self.production[j, t]
        return objective
        
    def _critere_ipopt(self, initial):
        """
        Private function who compute critere for lower bound
        
        _critere_ipopt(self,initial)
        """
        objective = initial
        for j in xrange(self.nb_obj):
            for t in xrange(self.time_hor):
                objective -= self.cost_setup[j, t]*self.setup[j, t]
        return objective
        
    def show_convergence(self):
        graphic(self.opt_lower, self.opt_upper, self.cycle)
        
    def solve(self):
        """
        Function who solve the PCLSP problem
        based on K. Haugen article
        
        solve(self)
        """
        diff = self.eps+ 1.
        count = 0
        start = time.clock()
        while ( (diff > self.eps) & (count < self.cycle) ):
            if self.verbose > 2:
                print "\nIterartion: "+str(count+1)
            previous_lambda = self.coef
            # compute upper bound
            self.price, self.setup = thomas(self.alpha,
                self.beta, self.cost_prod, self.cost_stor,
                self.cons_prod, self.cost_setup, self.coef,
                self.time_hor, self.nb_obj, self.verbose)
            upper = self._critere()
            # compute lower bound
            lower, self.coef = ipopt(self.alpha, self.beta,
                self.cost_prod, self.cost_stor, self.cons_prod,
                self.cost_setup, self.constraint, self.time_hor,
                self.nb_obj, self.verbose)
            lower = self._critere_ipopt(lower)
            # update lagrangian coefficients
            self.coef = self.smooth_param*previous_lambda\
                - (1-self.smooth_param)*self.coef
            # stor bounds
            self.opt_upper.append(upper)
            self.opt_lower.append(lower)
            # update 
            diff = fabs(upper-lower)
            count += 1
        end = time.clock()
        if self.verbose > 0:
            print "###RESULTS:"
            print "Number of iteration: "+str(count)
            print "Difference between upper and lower bound: "+str(diff)
            print "Last critere value: "+str(lower)
            print "Duration of computation:"+str(end-start)+"(s)"
        return self._critere ()

if __name__ == '__main__':
    TEST = PCLSP()
    TEST.solve()

