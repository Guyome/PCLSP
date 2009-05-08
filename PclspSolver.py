#!/usr/bin/env python
# -*- coding: utf-8 -*-
#       Copyright 2009 Guillaume Lanquepin <guillaume@himolde.no>
#       version 0.1

import sys
from optparse import OptionParser
import heuristic as he

def main():
    """
    Main function, deals with arguments and launch program
    """
    # Usual verifications and warnings
    if not sys.argv[1:]:
        print "You must specify at least input file\n"
        print "More help avalaible with -h or --help option"
        sys.exit(0)

    parser = OptionParser()
    parser.add_option("-v", "--verbose",
        action="store",type="int",dest="verbose",
        help="Write output informations (not only errors).",
        default=1)
    parser.add_option("-p", "--smooth-parameter",
         action="store", type="int", dest="theta",
        help="Coefficient for smoothing lagrangian coefficients.",
        default=0.5)
    parser.add_option("-e", "--epsilon",
         action="store", type="int", dest="eps",
        help="Minimal difference between upper and lower bounds",
        default=10)
    parser.add_option("-c", "--cycle",
        action="store", type="int", dest="cycle",
        help="Maximal number of cycle to compute optimum ",
        default=100)
    parser.add_option("-g", "--graphic",
        action="store_true", dest="graphic",
        help="Plot convergence of bounds ",
        default=False)
    parser.add_option("-f", "--file", 
        action="store", type="string", dest="file",
        help="File who containts data of the problem")
    (options, args) = parser.parse_args() 

    problem = he.HEURISTIC(options.theta,options.eps, options.cycle, options.verbose)
    problem.import_data(options.file)
    problem.solve()
    if options.graphic:
        problem.show_convergence()    
    sys.exit(0)
    
if __name__ == '__main__':
    main()

