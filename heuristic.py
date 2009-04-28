#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       Copyright 2009 Guillaume Lanquepin <guillaume@himolde.no>
#       version 0.1

import numpy as np
from tools import import_data
from plots import show_data

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
        self.constraint = import_data(dataFile)
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
            show_data(self)

    def main(self):

        return 0

if __name__ == '__main__':
    TEST = ClspHeuristic()

