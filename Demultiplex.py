#!/usr/bin/env python
"""
input a FASTQ file and return position and average read scores
"""


####################
# Import libraries #
####################

import argparse as ap
import re
import gzip
import numpy as np

#####################
# Argument Parser   #
#####################
def get_args():
	parser = ap.ArgumentParser()
	parser.add_argument("-f","--file", help="put the path the FASTQ file after '-f'"
	,type=str)
	parser.add_argument("-o","--out", help="put the path the output file after '-o'"
	,type=str)
	return(parser.parse_args())

####################
# GLOBAL VARIABLES #
####################

args = get_args()
file = args.file
outfile = args.out


#############
# FUNCTIONS #
#############

def convert_phred(letter):
    """Converts a single character into a phred score (integer)"""
    return(ord(letter)-33)

def qual_Length(file):
	"""
	Returns the length of quality
	Will be upper bound of 101 for read files
	and a lower bound of 8 for index files
	"""
	with gzip.open(file,"rt") as fh:
		for i,line in enumerate(fh): #loop through line in file
			line = line.rstrip() #remove new line char
			if i % 4 == 3: #quality scores 
				return(len(line))
	
def populate_mean(file):	
	"""
	This function takes the quality score of the FASTQ file
	and returns the mean for each position
	"""
	qLength = qual_Length(file)
	LNcount = 0
	all_qscores = [0.0 for x in range(qLength)]
	with gzip.open(file,"rt") as fh:
		for i,line in enumerate(fh): #loop through line in file
			LNcount = i + 1
			if i % 4 == 3: #quality scores
				line = line.rstrip() #remove new line char
				for j,score in enumerate(line): #loop through quality score in line
					all_qscores[j] += (convert_phred(score)) 
						
	return([(x/(LNcount/4)) for x in all_qscores])

def outwrite(string):
	"""
	input a string to be writen to an output file
	
	writes to file name
	returns None
	"""
	with open(outfile,"w") as fwrite:	
		fwrite.write(string)
		
	return None
	
def writeString(results):
	"""
	takes the results (means) and returns a string
	"""
	string = "# Base Pair\tMean Quality Score\n"

	for i,x in enumerate(results):
		string += str(i) + "\t" + str(x) + "\n"
	return(string)
	
def main():
	"""
	Calls functions to execute program
	"""
	means = populate_mean(file)
	printIt = writeString(means)
	outwrite(printIt)
	
########
# MAIN #
########
main()