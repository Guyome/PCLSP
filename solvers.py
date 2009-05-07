#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       Copyright 2009 Guillaume Lanquepin <guillaume@himolde.no>
#       version 0.1

import scipy.weave as wv
from scipy.weave import converters
import numpy as np
from os.path import join, split

def thomas(alpha, beta, cost_prod, cost_stor,cons_prod,
    cost_setup,coef, time_hor, nb_obj, verbose):
    """
    This function is an adapted version of J. Thomas's algorithm
    (see price-production decision with deterministic demand 1970)
    for lagrangian relaxatio and multi product 
    """
    opti_price = np.zeros((nb_obj,time_hor), float)
    extra_code = open(join(split(__file__)[0],'LUBP/lupb.cpp')).read()
    thomas_code = """
    int status = thomas(alpha, beta, cost_prod,
        cost_stor, cons_prod, cost_setup,
        opti_price, coef,
        time_hor, nb_obj, 1);
    """
    wv.inline( thomas_code ,
        ['time_hor', 'nb_obj','alpha', 'beta', 'cost_setup',
        'cost_prod', 'cost_stor', 'opti_price','cons_prod','coef','verbose'], 
        support_code=extra_code,
        type_converters=converters.blitz)
    return opti_price, np.array(opti_price > 0, int)
    
def ipopt(alpha, beta, cost_prod, cost_stor,cons_prod,
    cost_setup, constraint, time_hor, nb_obj, verbose):
    """
    This function solve quadratic problem as define
    in 'the profit mazimizing capacited lot-size problem'.
        It is based on ipopt solvers (http://www.coin-or.org/Ipopt/).
    """
    extra_code = open(join(split(__file__)[0],'LLBP/main.cpp')).read()
    results = np.zeros((nb_obj+1)*time_hor+1, float)
    code="""
    int status = ipopt(alpha,beta,cost_prod,
        cost_stor,cons_prod,constraint,
        results,time_hor,nb_obj,verbose);
    """
    wv.inline(code,['time_hor', 'nb_obj','alpha', 'beta', 'cost_setup',
        'cost_prod', 'cost_stor', 'results', 'cons_prod', 'constraint', 'verbose'],
        include_dirs=["LLBP/","/usr/include/coin/"],
        support_code=extra_code,
        libraries=['ipopt','lapack','pthread'],
        sources =['LLBP/LbIpopt.cpp'],
        type_converters=converters.blitz)
    return results[0],results[1:]
    
if __name__ == '__main__':
    time_hor = 3
    nb_obj =1
    alpha = np.array([100.,100.,100.]).reshape(nb_obj,time_hor)
    beta = np.array([1.,1.,1.]).reshape(nb_obj,time_hor)
    cost_prod = np.array([20.,20.,20.]).reshape(nb_obj,time_hor)
    cost_stor = np.array([2.,2.,2.]).reshape(nb_obj,time_hor)
    cons_prod = np.array([0.3,0.3,0.3]).reshape(nb_obj,time_hor)
    cost_setup = np.array([5.,5.,5.]).reshape(nb_obj,time_hor)
    constraint = np.array([40.,40.,40.])
    coef = np.zeros((nb_obj+1)*time_hor, float)
    verbose = 1
    print "THOMAS..."
    print thomas(alpha, beta, cost_prod, cost_stor,cons_prod,
    cost_setup,coef, time_hor, nb_obj, verbose)
    print "OK\n"
    print "IPOPT..."
    print ipopt(alpha, beta, cost_prod, cost_stor,cons_prod,
    cost_setup, constraint, time_hor, nb_obj, verbose)
    print "OK"
