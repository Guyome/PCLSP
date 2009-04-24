#!/usr/bin/env python
#
#       Copyright 2009 Guillaume Lanquepin <guillaume@himolde.no>
#		version 0.1

import numpy as np
from tools import importData

class ClspHeuristic:

	def __init__(self,dataFile="data_test.csv",_eps=0.1,_cycle=500):
		#problem paramters
		self.T,\ # time period
		self.J,\ # number of products
		self.alpha,\ # demand function parameter
		self.beta,\ # demand function parameter
		self.costSetup,\ # setup cost
		self.costStor,\ # storage cost
		self.costProd,\ # production cost
		self.consProd,\ # product consumption
		self.constraint = importData(dataFile) # production constraint
		#problem variables
		self.coef = np.zeros(self.T) # lagrangian multiplier
		self.production = np.zeros((self.J,self.T))
		self.price = np.zeros((self.J,self.T))
		self.setup = np.zeros((self.J,self.T))
		self.storage = np.zeros((self.J,self.T))
		#algorithm parameters
		self.eps = _eps #min difference between upper and lower bound (stop condition)
		self.cycle = _cycle #max number of iteration (stop condition)
		slef.optUpper = np.zeros(self.T) #vector fo optimun return by Upper
		slef.optLower = np.zeros(self.T) #vector fo optimun return by Lower

	def __verbose(self,start):
		if start:
			print "### Data"
			print "Period:"+self.T+"\tProduct:"+self.J
			print "Alpha:"
			print self.alpha
			print "Beta:"
			print self.beta
			print "Setup Cost:"
			print self.costS

	def __Upper(self):
		return 0
	def __Lower(self):
		return 0
	def main(self):

		return 0
