#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       Copyright 2009 Guillaume Lanquepin <guillaume@himolde.no>
#       version 0.1

import scipy.weave as wv
from scipy.weave import converters
import numpy as np
from os.path import join, split

def thomas(alpha, beta, cost_prod, cost_stor, cons_prod, 
    cost_setup, coef, time_hor, nb_obj, verbose):
    """
    This function is an adapted version of J. Thomas's algorithm
    (see price-production decision with deterministic demand 1970)
    for lagrangian relaxatio and multi product 
    
    thomas(alpha, beta, cost_prod, cost_stor,cons_prod,
    cost_setup,coef, time_hor, nb_obj, verbose)
    """
    opti_price = np.zeros((nb_obj,time_hor), float)
    extra_code = open(join(split(__file__)[0],'LUBP','lupb.cpp')).read()
    if verbose > 1:
        print "\tCompute upper bound (Thomas)..."
    thomas_code = """
    int status = thomas(alpha, beta, cost_prod,
        cost_stor, cons_prod, cost_setup,
        opti_price, coef,
        time_hor, nb_obj, verbose);
    """
    wv.inline( thomas_code ,
        ['time_hor', 'nb_obj','alpha', 'beta', 'cost_setup',
        'cost_prod', 'cost_stor', 'opti_price','cons_prod','coef','verbose'], 
        support_code=extra_code,
        type_converters=converters.blitz)
    return opti_price, np.array(opti_price > 0,float)
    
def ipopt(alpha, beta, cost_prod, cost_stor, cons_prod,
    cost_setup, setup, constraint, time_hor, nb_obj, verbose):
    """
    This function solve quadratic problem as define
    in 'the profit mazimizing capacited lot-size problem'.
    It is based on ipopt solvers (http://www.coin-or.org/Ipopt/).
    
    ipopt(alpha, beta, cost_prod, cost_stor,cons_prod,
    cost_setup, constraint, time_hor, nb_obj, verbose)
    """
    extra_code = open(join(split(__file__)[0],'LLBP','main.cpp')).read()
    results = np.zeros((4*nb_obj+1)*time_hor+1, float)
    if len(constraint.shape) > 1 :
        cons = constraint[0]
    else :
        cons = constraint
    if verbose > 1:
        print "\tCompute lower bound (IpOpt)..."
    code="""
    int status = ipopt(alpha,beta,cost_prod,
        cost_stor,cons_prod,setup,cons,
        results,time_hor,nb_obj,verbose);
    """
    wv.inline(code,['time_hor', 'nb_obj','alpha', 'beta', 'cost_setup',
        'cost_prod', 'cost_stor', 'setup', 'results', 'cons_prod', 'cons', 'verbose'],
        include_dirs=["LLBP","/usr/include/coin/"],
        support_code=extra_code,
        libraries=['ipopt','lapack','pthread'],
        sources =[join(split(__file__)[0],'LLBP','LbIpopt.cpp')],
        type_converters=converters.blitz)
    return results[0],results[1:(nb_obj+1)*time_hor+1],\
    results[(nb_obj+1)*time_hor+1:(2*nb_obj+1)*time_hor+1],\
    results[(2*nb_obj+1)*time_hor+1:(3*nb_obj+1)*time_hor+1],\
    results[(3*nb_obj+1)*time_hor+1:]

if __name__ == '__main__':
    time_hor = 3
    nb_obj =1
    alpha = np.array([100.,100.,100.]).reshape(nb_obj,time_hor)
    beta = np.array([1.,1.,1.]).reshape(nb_obj,time_hor)
    cost_prod = np.array([20.,20.,20.]).reshape(nb_obj,time_hor)
    cost_stor = np.array([2.,2.,2.]).reshape(nb_obj,time_hor)
    cons_prod = np.array([0.3,0.3,0.3]).reshape(nb_obj,time_hor)
    cost_setup = np.array([5.,5.,5.]).reshape(nb_obj,time_hor)
    setup = np.array([1.,0.,1.]).reshape(nb_obj,time_hor)
    constraint = np.array([40.,40.,40.])
    coef = np.zeros((nb_obj+1)*time_hor, float)
    verbose = 1
    print "THOMAS..."
    print thomas(alpha, beta, cost_prod, cost_stor,cons_prod,
    cost_setup,coef, time_hor, nb_obj, verbose)
    print "OK\n"
    print "IPOPT..."
    print ipopt(alpha, beta, cost_prod, cost_stor,cons_prod,
    cost_setup, setup, constraint, time_hor, nb_obj, verbose)
    print "OK"
