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
    obj = clsp.nb_obj
    alpha  = clsp.alpha
    beta = clsp.beta
    setup = clsp.cost_setup
    prod = clsp.cost_prod
    stor = clsp.cost_stor
    lamdba = clsp.coef
    cons = clsp.cons_prod
    optiprice = np.zeros((obj,hor), float)
    thomas_code = """
        //problem variable
        float c[hor];
        float f[hor];
        float price[hor][hor];
        float demand[hor][hor];
        // algorithm variable
        int k,i,t0,t = 0,ind[hor];
        float sumstor,sum,summin;
        
        //initiate price,demand since there is no storage
        for(int j = 0; j < obj;  j++)
        {
            price[t][t] = (alpha(j,t) + (prod(j,t) + cons(j,t)*lamdba(j,t))* beta(j,t)) / (2 * beta(j,t) );
            demand[t][t] = alpha(j,t) - beta(j,t) * price[t][t];
            c[t]= demand[t][t]*( price[t][t] - prod(j,t))-setup(j,t);
            f[t+1] = -c[t];
            ind[0] = 0;
            for ( t = 1; t < hor;  t++)
            {
                for (t0 = 0; t0 <= t; t0++)
                {
                    // sum of first storage
                    sumstor = 0;
                    for (i = t0; i < t; i++)
                    {
                        sumstor += stor(j,i);
                    }
                    // compute price and demand
                    price[t][t0] = (alpha(j,t) + 
                        (prod(j,t0) + sumstor + cons(j,t0)*lamdba(j,t0))* beta(j,t)) / (2 * beta(j,t));
                    demand[t][t0] = alpha(j,t) - beta(j,t) * price[t][t0];
                }
                for (t0 = 0; t0 <= t; t0++)
                {
                    sum = 0;
                    for (i = t0; i <= t; i++)
                    {
                        sumstor = 0;
                        for (k = t0; k < i; k++)
                        {
                            sumstor += stor(j,i);
                        }
                        sum += (prod(j,t0) + sumstor - price[j][t0])*demand[j][t0];
                    }
                    c[t0] = sum + setup(j,t0);
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
                optiprice(j,t)=price[t][ind[t]];
            }
        }
    """
    wv.inline( thomas_code , ['hor', 'obj','alpha', 'beta', 'setup', 'prod', 'stor', 'optiprice','cons','lamdba'], \
    type_converters=converters.blitz)
    return optiprice, np.array(optiprice > 0, int)
