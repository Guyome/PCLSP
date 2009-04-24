#!/usr/bin/env python
#       
#       Copyright 2009 Guillaume Lanquepin <guillaume@himolde.no>
#		version 0.1

import numpy as np

class ClspHeuristic:
	
	#_alpha,_beta,_prod,_stor,_cost,_constraint,_item
	def __init__(self,_dataFolder="data_test.csv",_eps=0.1,_cycle=500):
		#problem paramters
		self.T = _T 
		self.J = _J
		#problem variables
		coef = np.zeros
		#algorithm parameters
		self.eps = _eps #min difference between upper and lower bound (stop condition)
		self.cycle = _cycle #max number of iteration (stop condition)
		slef.optUpper = np.zeros(self.T) #vector fo optimun return by Upper
		slef.optLower = np.zeros(self.T) #vector fo optimun return by Lower
		

	def Upper(self):
		
	def Lower(self):
		
	def main(self):
		
		return 0

if __name__ == '__main__': main()
