#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       Copyright 2009 Guillaume Lanquepin <guillaume@himolde.no>
#       version 0.1
from tools import *

def show_data(clsp):
    print "### Data"
    print "Period:"+str(clsp.t)+"\tProduct:"+str(clsp.J)
    print "Alpha:",
    print_tab(clsp.alpha)
    print "Beta:",
    print_tab(clsp.beta)
    print "Setup Cost:",
    print_tab(clsp.costSetup)
    print "Storage cost:",
    print_tab(clsp.costStor)
    print "Prodcution cost:",
    print_tab(clsp.costProd)
    print "Product consumption:",
    print_tab(clsp.constraint)
    print "Production constraint",
    print_tab(clsp.constraint)
