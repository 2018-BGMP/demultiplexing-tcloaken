#!/usr/bin/env python3
"""
input gzipped FASTQ files and return demultiplexed 
Fastq files and a report file and two files for undetermined


Includes N recovery.


"""


####################
# Import libraries #
####################

import argparse as ap
#import re
import gzip
#import numpy as np
from itertools import izip

#####################
# Argument Parser   #
#####################
def get_args():
	parser = ap.ArgumentParser()
	parser.add_argument("-R1","--Readfile1", help="put the path the FASTQ file after the flag followed by a space for the forward sequence read"
	,type=str)
	parser.add_argument("-R2","--Readfile2", help="put the path the FASTQ file after the flag followed by a space for the forward index read"
	,type=str)
	parser.add_argument("-R3","--Readfile3", help="put the path the FASTQ file after the flag followed by a space for the reverse index read"
	,type=str)
	parser.add_argument("-R4","--Readfile4", help="put the path the FASTQ file after the flag followed by a space for the reverse sequence read"
	,type=str)
	parser.add_argument("-i","--index", help="put the path the index file after '-i', The index file is list with a name for the barcode followed by a space followed by the barcode followed by a newline char"
	,type=str)
	parser.add_argument("-N","--Nrecovery", help="N recovery number, -1 no N recovery, 0 is based on minimum hamming distance of the barcodes, otherwise number of n's to recover, default is 0"
	,type=int,default=0)
	parser.add_argument("-s","--QualityScoreThreshold",help="The phred score that the average read quality score must be over to be accepted, default is 30"
	,type=int,default=30)
	
	return(parser.parse_args())

####################
# GLOBAL VARIABLES #
####################

args = get_args()
R1 = args.Readfile1 #forward sequence read
R2 = args.Readfile2 #forward index
R3 = args.Readfile3 #reverse index
R4 = args.Readfile4 #reverse sequence read
index = args.index  #index list, example: "A5 ACGGATAC"
Nreco = args.Nrecovery #number of N's to recover
phredThres = args.QualityScoreThreshold #theshold phred score


#############
# FUNCTIONS #
#############

def hamming_distance(s1, s2):
    #Return the Hamming distance between equal-length sequences
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

def findHamMin(list):
	"""
	takes a list of the same length strings
	and returns the minimum hamming distance
	"""
	totalStrings = len(list)
	check = 0
	listRes = []
	for elm in list:
		for elm2 in list:
			if elm != elm2:
				listRes.append(hamming_distance(elm,elm2))
			else:
				check += 1
				
	if check > totalStrings:
		raise StandardError("Two of the same barcodes")
	return min(listRes)
	
def makeIndex(index):
	"""propogate index dictionary"""
	inDextionary = {}
	with open(index,"r") as fh:
		for line in fh:
			line = line.rstrip() #remove \n newline char
			line = line.split(" ") #parse space between barcode and name for barcode
			inDextionary[line[1]] = line[0]
	return inDextionary

def convert_phred(letter):
    """Converts a single character into a phred score (integer)"""
    return(ord(letter)-33)

	
def RevComp(nucleotides):
	"""
	Given a DNA sequence (nucleotides) return the reverse compliment
	"""
	swap = {"A":"T","T":"A","G":"C","C":"G","N":"N"}#dictionary
	return ''.join([ swap[x] for x in nucleotides[::-1]])
	

	


###########################
#FUNCTIONS FOR N RECOVERY #
###########################

def permutides(n):
	
	""" permutes Nucleotides and returns a list of the columns
	n = 1   n=2   n=3  
	A		A A  A A A   ...
	C		A C  A A C
	T		A T  A A T 
	G		A G  A A G
			C A  A C A
			C C  A C C 
			C T  A C T
			C G  A C G 
			T A  A T .
			T C  A T .
			T T  A T .
			T G  A T
			G A  A G
			G C  A G
			G T  A G
			G G  A G	
	"""	
	nucl = ["A","C","T","G"]
	bigSet = []
	for i in range(n,0,-1):
		bigSet.append("".join([x*4**(i-1) for x in nucl]*(4**(n-i))))
	return bigSet
	
def buildSet(barcode,NrecMin):
	"""given a barcode with at least one 'N' then we build a
	set	with different nucleotides to replace the 'Ns'"""
	#set of possibilities that the string could be
	#if it were actually nucleotides instead of "N"s
	setN = set()
	#get positions in barcode where "N"s are 
	pos = [i for i, ltr in enumerate(barcode) if ltr == "N"] 
	n = len(pos) #the number of "N"s
	if n >= NrecMin:
		#too many N's, abort
		return setN #return an empty set
	bar = list(barcode) #make it mutable
	listN = permutides(n) #get permutation list
	#loop through number of permutations (4 ^ # of N's in barcode)
	for i in range(4**n):
		#loop through number of positions of "N"
		for j,x in enumerate(pos):
			#Change each "N" to a possibility of nucleotide
			bar[x] = listN[j][i] 
		#add it to the set of possibilities
		setN.add("".join(bar))
		
	return setN #return set of possibilities

def Nrecovery(setPoss,indexes,Nrec):
	"""given a set of possibilities and a dictionary of 
	indexes then return the intersection and a bool of 
	whether or not the intersection as one element
	"""
	setV = set(indexes.keys()) #make indexes into a set
	#intersection of both sets
	intersect = setPoss.intersection(setV)
	if len(intersect) == 1:
		return intersect, True
	else:
		return intersect, False
		

		
##########################################################
# END of N recovery ######################################
##########################################################

def qualCutOff(seqRead,QC):
	"""given a sequence calulate the mean quality score,
	if it is above a given quality cutoff then return True
	else return False"""
	return (sum([convert_phred(x) for x in seqRead])/len(seqRead) > QC)
	

def indexMatching(R2,R3,index):
	"""
	given a couple of barcodes, return whether or not they
	match up.  If there's "N"s they won't match since
	index should not have any "N"s 
	"""
	return ((R2 == RevComp(R3)) and (R2 in index))
	
	

def Demult(index, Ngo,QScore):
	"""
	Given a index dictionary
	        Ngo : number for N recover (negative means no N recovery)
			QScore : number for Lowerbound threshold average quality score
			
	Sort fastq files by index dictionary, and apply
	quality assement,  if below quality cutoff then
	the reads will be in "undetermined" based on from either
	R1 or R4
	
	returns 4 dictionaries: quality assessed R1 reads where barcodes R2 and R3 match
							quality assessed R4 reads where barcodes R2 and R3 match
							undetermined R1 and R2 reads where either R2 and R3 did not match or average quality score of the read did not exceed threshold
							report file, csv
							
	"""
	QSlowerBound = QScore# threshold quality score (strictly greater than)
	files = {}  	  	 # read 1 
	files2 = {} 	  	 # read 2
	files3 = {} 	  	 # undetermined reads 1 and 2
	file4 = {}  	 	 # report file
	Ns1 = set() 		 # Read count if an "N" is present in barcode1
	Ns2 = set()     	 # ""  ""  barcode1
	NsRecovered = set()  # N's recovery succesfully from barcode1 or barcode2
	badQual1 = set()     # average quality score below threshold for R1
	badQual2 = set()	 # ""   ""  R4
	indexHopping = set() # barcodes do not match, putative index hopping
	EmployNRecovery = (Ngo >= 0) #bool if true, will employ N recovery
	LNcount     = 0 #line counter
	ReadCounter = 0 #updated every 4 line counts
	
	# instantiate dictionary with strings according to key
	# and do it for both read 1 and read 2
	# and two separate undetermined read 1 and read 2
	for key in index.values():
		files[key+"_R1"] = ""
	for key in index.values():
		files2[key+"_R2"] = ""	
	countData = {}
	for key in files:
		countData[key] = 0
	for key in files2:
		countData[key] = 0

	files3["Undetermined_R1"] = ""
	files3["Undetermined_R2"] = ""
	countData["Undetermined_R1"] = 0
	countData["Undetermined_R2"] = 0
	#open up all the files
	with gzip.open(R1,"rt") as fR1, gzip.open(R2,"rt") as fR2, gzip.open(R3,"rt") as fR3, gzip.open(R4,"rt") as fR4:
		
		#keey track of these variables through the loop to come
		qualScores1 = ""  #R1
		qualScores2 = ""  #R4
		barcode1    = ""  #R1
		barcode2    = ""  #R4
		header1     = ""  #R1
		header2     = ""  #R4
		message1    = ""  #R1
		message2    = ""  #R4
		Seq1        = ""  #R1
		Seq2        = ""  #R4
		BarcodesMatch = False #barcodes for this read match bool
		currentKey  = ""  #current key associated with barcode
		
		#loop through each line in all the files 
		for lR1,lR2,lR3,lR4 in izip(fR1,fR2,fR3,fR4): 
			
			if LNcount % 4 == 0: #header
				ReadCounter += 1 #increment read counter
				#update current headers
				header1 = lR1.rstrip()
				header2 = lR4.rstrip()
				
				
			if LNcount % 4 == 1: #sequence
				#update current sequence
				Seq1 = lR1
				Seq2 = lR4
				#get barcodes
				barcode1 = lR2.rstrip()
				barcode2 = lR3.rstrip()
				#are there "N"s in the barcodes?
				Nbar1 = ("N" in barcode1) 
				Nbar2 = ("N" in barcode2)
				if Nbar1: #if theres N's
					Ns1.add(ReadCounter)
					if EmployNRecovery: 
						
						# employ recovery
						intersection1,bool1 = Nrecovery(buildSet(barcode1,Ngo),index,Ngo)
						if bool1: #if succesful
							#update barcode1
							barcode1 = list(intersection1)[0]
							
						#else:
							 #"Can't unravel N's in barcode1"
				
				if Nbar2: #if theres N's
					Ns2.add(ReadCounter)
					if EmployNRecovery:
						
						# employ recovery
						intersection2,bool2 = Nrecovery(buildSet(RevComp(barcode2),Ngo),index,Ngo)
						if bool2: #if succesful
							#update barcode2
							barcode2 = RevComp(list(intersection2)[0])
							
						#else:
							#"Can't unravel N's in barcode2"
				#check to see if the index matches
				if not indexMatching(barcode1,barcode2,index):
					#  NO MATCH: putative index hopping
					indexHopping.add(ReadCounter)
					BarcodesMatch = False
					
				else:
					#they MATCH! good stuff
					BarcodesMatch = True
					currentKey = index[barcode1]
					if Nbar1 or Nbar2:
						#recoverend an N
						NsRecovered.add(ReadCounter)
			
			if LNcount % 4 == 2: #message
				#update current messages
				message1 = lR1
				message2 = lR4
				
			if LNcount % 4 == 3: #quality scores
				qualScores1 = lR1.rstrip() #remove new line char
				qualScores2 = lR4.rstrip()
				if not (qualCutOff(qualScores1,QSlowerBound)):
					#if it did not pass average quality score test
					#add the read counter to the set of bad quality for
					#R1
					badQual1.add(ReadCounter)
					#doesn't meet quality
					files3["Undetermined_R1"] += header1+":"+barcode1+"\n"+ Seq1 + message1 + lR1
					countData["Undetermined_R1"] += 1
				elif BarcodesMatch:
					#meets quality, and barcodes match
					files[currentKey+"_R1"] += header1+":"+barcode1+"\n" + Seq1 + message1 + lR1
					countData[currentKey+"_R1"] += 1
				else:
					#quality scores ok but not matching barcode
					files3["Undetermined_R1"] += header1+":"+barcode1+"\n" + Seq1 + message1 + lR1
					countData["Undetermined_R1"] += 1
				if not (qualCutOff(qualScores2,QSlowerBound)):
					#similiarly for R4
					badQual2.add(ReadCounter)
					files3["Undetermined_R2"] += header2+":"+RevComp(barcode2)+"\n" + Seq2 + message2 + lR4
					countData["Undetermined_R2"] += 1
				elif BarcodesMatch:
					#quality scores ok and matches not
					files2[currentKey+"_R2"] += header2+":"+barcode1+"\n" + Seq2 + message2 + lR4
					countData[currentKey+"_R2"] += 1
				else:
					#quality scores but not matching barcode
					files3["Undetermined_R2"] += header2 +":"+RevComp(barcode2)+"\n"+ Seq2 + message2 + lR4
					countData["Undetermined_R2"] += 1
			LNcount += 1
	
	#count data	
	#and prepare report file
	for key,value in countData.items():
		file4[key+" % reads"] = (float(value)/float(ReadCounter))*100 #percent of total reads
	file4["_Min_Hamming_Distance_of_Barcodes"] = findHamMin(list(index.keys()))
	Ns1 = len(list(Ns1))
	file4["_Ns in R1 barcode"] = Ns1
	Ns2 = len(list(Ns2))
	file4["_Ns in R2 barcode"] = Ns2
	NsRecovered = len(list(NsRecovered))
	file4["_Ns_Recovered"] = NsRecovered
	badQ1 = len(list(badQual1))
	file4["_bad_ave_quality_R1"] = badQ1
	badQ2 = len(list(badQual2))
	file4["_bad_ave_quality_R2"] = badQ2
	indexHopping = len(list(indexHopping))
	file4["_Putative_Index_hopping"] = indexHopping
	file4["_Ns in R1 barcode %"] = (float(Ns1)/float(ReadCounter))*100
	file4["_Ns in R2 barcode %"] = (float(Ns2)/float(ReadCounter))*100	
	file4["_Putative_Index_hopping %"] = (float(indexHopping)/float(ReadCounter))*100
	file4["_bad_ave_quality_R1 %"] = (float(badQ1)/float(ReadCounter))*100
	file4["_bad_ave_quality_R2 %"] = (float(badQ2)/float(ReadCounter))*100
	return files,files2,files3,file4

def writeFiles(dictionary1,dictionary2,dictionary3,dictionary4):
	"""
	Given four dictionarys of strings and numbers
	write files, (number of indexes twice + 2)
	and write a report file in CSV format
	return None
	"""
	
	for key,value in dictionary1.items():
		with gzip.open(key,"wt") as fh:
			fh.write(value)
	for key,value in dictionary2.items():
		with gzip.open(key,"wt") as fh:
			fh.write(value)
	for key, value in dictionary3.items():
		with gzip.open(key,"wt") as fh:
			fh.write(value)
	#write report
	string = ""
	
	for key in sorted(dictionary4.iterkeys()):
		string += key + "\t" + str(dictionary4[key]) + "\n"

	with open("report","w") as fh:
		fh.write(string)
	
	return None
	
def Nrecos(number,index):
	"""
	give number from the args parse
	return number for the "N" recovery
	"""
	if number < 0:
		return -1
	elif number == 0:
		return findHamMin(list(index.keys()))
	else:
		return number
		
def main():
	"""
	Calls functions to execute program
	"""
	Index = makeIndex(index) #get Indexes together
	Ngo = Nrecos(Nreco,Index) #get "N" recovery information
	fileR1,fileR2,fileUndeterminedR1R2,fileReport = Demult(Index,Ngo,phredThres)
	writeFiles(fileR1,fileR2,fileUndeterminedR1R2,fileReport)
	
########
# MAIN #
########
main()
