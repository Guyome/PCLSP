#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       Copyright 2009 Guillaume Lanquepin <guillaume@himolde.no>
#       version 0.1

import numpy as np
from lupb import thomas
from llpb import ipopt
from tools import import_data
from plots import show_data

class clspHeuristic:

    def __init__(self, data_file="data_test.csv", _eps=0.1, _cycle=500,verbose=True):
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
        self.coef = np.zeros(self.time_hor,float) # lagrangian multiplier
        self.production = np.zeros((self.nb_obj, self.time_hor),float)
        self.price = np.zeros((self.nb_obj, self.time_hor),float)
        self.setup = np.zeros((self.nb_obj, self.time_hor),float)
        self.storage = np.zeros((self.nb_obj, self.time_hor),float)
        #algorithm parameters
        self.eps = _eps #min difference between upper and lower bound (stop condition)
        self.cycle = _cycle #max number of iteration (stop condition)
        self.opt_upper = np.zeros(self.time_hor) #vector fo optimun return by Upper
        self.opt_lower = np.zeros(self.time_hor) #vector fo optimun return by Lower
        if verbose:
            show_data(self)

    def main(self):
        upper,delta = thomas(self)
        ipopt(self)
        print upper,delta
        return 0

if __name__ == '__main__':
    TEST = clspHeuristic()
    TEST.main()

