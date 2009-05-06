#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       Copyright 2009 Guillaume Lanquepin <guillaume@himolde.no>
#       version 0.1

import scipy.weave as wv
from scipy.weave import converters
import numpy as np
from os.path import join, split

def ipopt(clsp):

    hor = clsp.time_hor
    obj = clsp.nb_obj
    alpha  = clsp.alpha
    beta = clsp.beta
    prod = clsp.cost_prod
    stor = clsp.cost_stor
    consprod = clsp.cons_prod
    constraint = clsp.constraint
    results = np.zeros(((obj+1)*hor+1), float)
    print results
    extra_code = open(join(split(__file__)[0],'LLBP/main.cpp')).read()
    verbose_level = 5
    code="""
        float r[hor];
        float** al = new float*[hor];//coeff demand function
        float** b = new float* [hor];//coeff demand function
        float** v = new float*[hor];//production cost
        float** h = new float*[hor];//storage cost
        float** a = new float*[hor];//consumption
        int status;
        for (int i = 0; i < hor; i ++)
        {
            r[i]=(float)constraint(0,i);
            al[i] = new float[obj];
            b[i] = new float[obj];
            v[i] = new float[obj];
            h[i] = new float[obj];
            a[i] = new float[obj];
            printf("r[%d]= %f\\n",i,r[i]);
            for(int j = 0; j < obj; j ++)
            {
                al[i][j] = alpha(j,i);
                b[i][j] = beta(j,i);
                v[i][j] = prod(j,i);
                h[i][j] = stor(j,i);
                a[i][j] = consprod(j,i);
                printf("al[%d][%d]= %f\\n",i,j,al[i][j]);
                printf("b[%d][%d]= %f\\n",i,j,b[i][j]);
                printf("v[%d][%d]= %f\\n",i,j,v[i][j]);
                printf("h[%d][%d]= %f\\n",i,j,h[i][j]);
                printf("a[%d][%d]= %f\\n",i,j,a[i][j]);
            }
        }
        // Create an instance of your nlp...
        status = ipopt(al,b,v,h,a,r,results,hor,obj,verbose_level);
        printf("END\\n");
    """
    wv.inline(code,["alpha","beta","prod","stor","consprod","constraint","hor","obj","results",'verbose_level'],\
    include_dirs=["LLBP/","/usr/include/coin/"],\
    support_code=extra_code,\
    libraries=['ipopt','lapack','pthread'],\
    sources =['LLBP/LbIpopt.cpp'],\
    type_converters=converters.blitz)
    return results[0],results[1:]

def heuristic(clsp):
    return [0],[100]
