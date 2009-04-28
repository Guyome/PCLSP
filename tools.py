#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       Copyright 2009 Guillaume Lanquepin <guillaume@himolde.no>
#       version 0.1
#
# This file containts all tools needed by heuristic program
import csv
import numpy as np

def import_data(address):
    try:
        inputcsv = csv.reader(open(address, "r"), delimiter=";", lineterminator="\n")
    except IOError:
        print "I/O error with the specified file"

    data = list() # all data
    item = list() # each tabular
    count = 0
    for row in inputcsv:
        if count < 2 : # read Time period and number of product
            data.append(int(row[1]))
        else :
            item.append(row[1:])
            if ((count-2)%data[1] == 1):
                data.append(np.array(item, dtype=float))
                item = list()
        count += 1
    data.append(np.array(item, dtype=float)) # manage the last tabular
    return data

def print_tab(tabular):
    row,col = tabular.shape
    tab = ""
    for i in xrange(row):
        tab += "\n\t|"
        for j in xrange(col):
            tab += str(tabular[i][j])+"|"
    print tab


if __name__ == '__main__':
    test = import_data("data_test.csv")
    print test

