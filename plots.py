#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       Copyright 2009 Guillaume Lanquepin <guillaume@himolde.no>
#       version 0.1

def print_tab(tabular):
    row,col = tabular.shape
    tab = ""
    for i in xrange(row):
        tab += "\n\t|"
        for j in xrange(col):
            tab += str(tabular[i][j])+"|"
    print tab

def show_data(clsp):
    print "### Data"
    print "Period:"+str(clsp.time_hor)+"\tProduct:"+str(clsp.nb_obj)
    print "Alpha:",
    print_tab(clsp.alpha)
    print "Beta:",
    print_tab(clsp.beta)
    print "Setup Cost:",
    print_tab(clsp.cost_setup)
    print "Storage cost:",
    print_tab(clsp.cost_stor)
    print "Prodcution cost:",
    print_tab(clsp.cost_prod)
    print "Product consumption:",
    print_tab(clsp.cons_prod)
    print "Production constraint",
    print_tab(clsp.constraint)
