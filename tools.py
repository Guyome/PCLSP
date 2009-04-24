#!/usr/bin/env python
#
#		Copyright 2009 Guillaume Lanquepin <guillaume@himolde.no>
#		version 0.1
#
# This file containts all tools needed by heuristic program
import csv
import numpy as np

def importData(address):
		try:
			input= open(address, "r")
			inputCsv=csv.reader(input,delimiter=";",lineterminator="\n")
		except IOError:
			print "I/O error with the specified file"

		data = list() # all data
		item = list() # each tabular
		count = 0
		for row in inputCsv:
			if count < 2 : # read Time period and number of product
				data.append(int(row[1]))
			else :
				item.append(row[1:])
				if ((count-2)%data[1] == 1):
					data.append(np.array(item,dtype=float))
					item = list()
			count += 1
		data.append(np.array(item,dtype=float)) # manage the last tabular
		return data

if __name__ == '__main__':
	test = importData("data_test.csv")
	print test
