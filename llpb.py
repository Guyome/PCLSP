#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       Copyright 2009 Guillaume Lanquepin <guillaume@himolde.no>
#       version 0.1

import scipy.weave as wv
from scipy.weave import converters
import numpy as np

def ipopt(clsp):

    hor = clsp.time_hor
    obj = clsp.nb_obj
    alpha  = clsp.alpha
    beta = clsp.beta
    prod = clsp.cost_prod
    stor = clsp.cost_stor
    consprod = clsp.cons_prod
    constraint = clsp.constraint
    ret = np.zeros(hor,float)
    code="""
        float r[hor];
        float** al = new float*[hor];//coeff demand function
        float** b = new float* [hor];//coeff demand function
        float** v = new float*[hor];//production cost
        float** h = new float*[hor];//storage cost
        float** a = new float*[hor];//consumption
        for (int i = 0; i < hor; i ++)
        {
            r[i]=(float)constraint(0,i);
            al[i] = new float[obj];
            b[i] = new float[obj];
            v[i] = new float[obj];
            h[i] = new float[obj];
            a[i] = new float[obj];
            for(int j = 0; j < obj; j ++)
            {
                al[i][j] = alpha(j,i);
                b[i][j] = beta(j,i);
                v[i][j] = prod(j,i);
                h[i][j] = stor(j,i);
                a[i][j] = consprod(j,i);
            }
        }
        // Create an instance of your nlp...
        SmartPtr<TNLP> mynlp = new FreeSolver(al,b,v,h,a,r,hor,obj);
        
        // Create an instance of the IpoptApplication
        SmartPtr<IpoptApplication> app = new IpoptApplication();
        
        // Initialize the IpoptApplication and process the options
        ApplicationReturnStatus status;
        
        status = app->Initialize();
        if (status != Solve_Succeeded) 
        {
            printf("\\n\\n*** Error during initialization!\\n");
        }

        status = app->OptimizeTNLP(mynlp);

        if (status == Solve_Succeeded) 
        {
            // Retrieve some statistics about the solve
            Index iter_count = app->Statistics()->IterationCount();
            printf("\\n\\n*** The problem solved in %d iterations!\\n", iter_count);
            Number final_obj = app->Statistics()->FinalObjective();
            printf("\\n\\n*** The final value of the objective function is %e.\\n", -final_obj);
        }
    """
    wv.inline(code,["alpha","beta","prod","stor","consprod","constraint","hor","obj"],\
    include_dirs=["LLBP/","/usr/include/coin/"],\
    library_dirs=["/usr/lib/"],\
    libraries=['ipopt','lapack','pthread','dl'],\
    sources =['LLBP/FreeSolver.cpp'],\
    headers=['"FreeSolver.hpp"','"IpSolveStatistics.hpp"','"IpIpoptApplication.hpp"'],\
    type_converters=converters.blitz,verbose=3,compiler='gcc')


