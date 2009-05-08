#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       Copyright 2009 Guillaume Lanquepin <guillaume@himolde.no>
#       version 0.1

import scipy.weave as wv
from scipy.weave import converters
import numpy as np

def thomas(clsp):
    
    hor = clsp.time_hor
    alpha  = clsp.alpha
    beta = clsp.beta
    setup = clsp.cost_setup
    prod = clsp.cost_prod
    stor = clsp.cost_stor
    optiprice = np.zeros(hor, float)
    thomas_code = """
        //problem variable
        float c[hor];
        float f[hor];
        float price[hor][hor];
        float demand[hor][hor];
        // algorithm variable
        int k,j,t0,t = 0,ind[hor];
        float sumstor,sum,summin;
        
        //initiate price,demand since there is no storage
        price[t][t] = (alpha(t) + prod(t) * beta(t)) / (2 * beta(t) );
        demand[t][t] = alpha(t) - beta(t) * price[t][t];
        c[t]= demand[t][t]*( price[t][t] - prod(t))-setup(t);
        f[t+1] = -c[t];
        ind[0] = 0;
        for ( t = 1; t < hor;  t++)
        {
            for (t0 = 0; t0 <= t; t0++)
            {
                // sum of first storage
                sumstor = 0;
                for (j = t0; j < t; j++)
                {
                    sumstor += stor(j);
                }
                // compute price and demand
                price[t][t0] = (alpha(t) + (prod(t0) + sumstor)* beta(t)) / (2 * beta(t) );
                demand[t][t0] = alpha(t) - beta(t) * price[t][t0];
            }
            for (t0 = 0; t0 <= t; t0++)
            {
                sum = 0;
                for (j = t0; j <= t; j++)
                {
                    sumstor = 0;
                    for (k = t0; k < j; k++)
                    {
                        sumstor += stor(j);
                    }
                    sum += (prod(t0) + sumstor - price[j][t0])*demand[j][t0];
                }
                c[t0] = sum + setup(t0);
            }
            // compute minimal criterium
            summin = 0;
            for (t0 = 0; t0 <= t; t0++)
            {
                if (c[t0] + f[t0] < summin)
                {
                    summin = c[t0] + f[t0];
                    ind[t] = t0;
                }
            }
            f[t+1] = summin;
        }
        for( t = 0; t < hor;  t++)
        {
            optiprice(t)=price[t][ind[t]];
        }
    """
    wv.inline( thomas_code , ['hor', 'alpha', 'beta', 'setup', 'prod', 'stor', 'optiprice'], \
    type_converters=converters.blitz)
    return optiprice, np.array(optiprice > 0, int)
