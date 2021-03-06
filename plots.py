#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       Copyright 2009 Guillaume Lanquepin <guillaume@himolde.no>
#       version 0.1
import matplotlib.pyplot as plt
import numpy as np

def print_tab(tabular):
    """
    This function plot an tabular
    """
    try:
        row,col = tabular.shape
        tab = ""
        for i in xrange(row):
            tab += "\n"
            line = ""
            for j in xrange(col):
                line += "--------"
                tab += str(tabular[i][j])+"\t"
    except:
        tab = "\n"
        line = ""
        for i in xrange(len(tabular)):
            line += "--------"
            tab += str(tabular[i])+"\t"
    print line,
    print tab
    print line

def show_data(clsp):
    print "### Data"
    print "Period:"+str(clsp.time_hor)+"\tProduct:"+str(clsp.nb_obj)
    print "Demand function (slope):"
    print_tab(clsp.alpha)
    print "Demand function (y-intercept,):"
    print_tab(clsp.beta)
    print "Setup Cost:"
    print_tab(clsp.cost_setup)
    print "Holding cost:"
    print_tab(clsp.cost_stor)
    print "Prodcution cost:"
    print_tab(clsp.cost_prod)
    print "Product consumption:"
    print_tab(clsp.cons_prod)
    print "Production constraint:"
    print_tab(clsp.constraint)
    
def show_varaibles(clsp):
    demand = np.zeros((clsp.nb_obj,clsp.time_hor))
    for j in xrange(clsp.nb_obj):
            for t in xrange(clsp.time_hor):
                demand[j, t] = clsp.alpha[j, t]-clsp.beta[j, t]*clsp.price[j, t]
    print"### Variable"
    print "Demand:"
    print_tab(demand)
    print "Storage:"
    print_tab(clsp.storage)
    print "Production:"
    print_tab(clsp.production)
    print "Setup structure:"
    print_tab(clsp.setup>0)
    
def graphic(lower, upper):
    ind = np.arange(len(lower))
    plt.plot(ind, lower, 'r--',label='Lower bound')
    plt.plot(ind, upper, 'b--',label='Upper bound')
    plt.xlabel('Number of iteration')
    plt.ylabel('Profit')
    plt.legend()
    plt.title('Convergence of upper and lower bounds')
    plt.show()

        
