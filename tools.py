#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       Copyright 2009 Guillaume Lanquepin <guillaume@himolde.no>
#       version 0.1
#
# This file containts all tools needed by heuristic program
import csv
from sys import exit
import numpy as np

def import_data(address):
    """
    This function import data form cvs file define as 
    in the following example:
    
    "T";4
    "J";2
    "Alpha";100;100;100;100
    ;100;100;100;100
    "Beta";1;1;1;1
    ;1;1;1;1
    "Setup cost";3;3;3;3
    ;3;3;3;3
    "Storage cost";2;2;2;2
    ;2;2;2;2
    "Production cost";20;20;20;20
    ;2;2;2;2
    "Consuption";0.3;0.3;0.3;0.3
    ;2;2;2;2
    "Constraint";35;40;42;24
    """
    try:
        inputcsv = csv.reader(open(address, "r"), delimiter=";", lineterminator="\n")
    except IOError:
        print "File not exists or is unreadable, please check it."
        exit(1)

    data = list() # all data
    item = list() # each tabular
    count = 0
    subcount = 0
    try:
        for row in inputcsv:
            if count < 2 : # read Time period and number of product
                data.append(int(row[1]))
            else :
                item.append(row[1:])
                subcount +=1 
                if subcount == data[1]:
                    data.append(np.array(item, dtype=float))
                    item = list()
                    subcount = 0
            count += 1
        if (data[1] > 1):
            data.append(np.array(item, dtype=float)) # manage the last tabular
    except:
        print "File is not well formated, please correct it."
        exit(1)
    return data

if __name__ == '__main__':
    TEST = import_data("data_test.csv")
    print TEST

