#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       Copyright 2009 Guillaume Lanquepin <guillaume@himolde.no>
#       version 0.1

import scipy.weave as wv
from scipy.weave import converters
import numpy as np
from os.path import join, split

def thomas(alpha, beta, cost_prod,
    cost_stor, cons_prod, cost_setup,
    opti_price, coef, time_hor, nb_obj, verbose):
    """
    This function is an adapted version of J. Thomas's algorithm
    (see price-production decision with deterministic demand 1970)
    for lagrangian relaxatio and multi product 
    """
    opti_price = np.zeros((clsp.nb_obj,clsp.time_hor), float)
    extra_code = open(join(split(__file__)[0],'LUBP/lupb.cpp')).read()
    thomas_code = """
    int status = ipopt(alpha, beta, cost_prod,
        cost_stor, cons_prod, cost_setup,
        opti_price, coef,
        time_hor, nb_obj, 1);
    """
    wv.inline( thomas_code ,\
    ['time_hor', 'nb_obj','alpha', 'beta', 'cost_setup',\
    'cost_prod', 'cost_stor', 'opti_price','cons_prod','coef','verbose'], \
    support_code=extra_code,\
    type_converters=converters.blitz,verbose=2)
    return opti_price, np.array(opti_price > 0, int)
    
if __name__ == '__main__':
    from heuristic import *
    TEST = PCLSP()
    print thomas(TEST)
