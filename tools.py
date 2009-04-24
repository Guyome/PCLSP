#!/usr/bin/env python
#       
#       Copyright 2009 Guillaume Lanquepin <guillaume@himolde.no>
#		version 0.1
import csv

def importData(adress):
		try:
			input= open(address, "r")
			inputCsv=csv.reader(input,delimiter=";",lineterminator="\n")
		except IOError:
			print "I/O error with the specified file"
			
		data = list()
		item = list()
		count = 0
		for row in inputCsv:
			if count < 2 :
				data.append(row[1])
				count =+ 1
			else
				item.append(row[1:]])
				if ((count-2)%data[2] == 0):
					data.append(item)
					item = list()
				
