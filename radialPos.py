""" 
Author: Jason Wherry	Date: 2/1/2020		To Run: python3 radial.py

OBJ:
	Write a program to read the file and compute the following:
	1.	Number of genomic windows --> 90,877
	2.	Number of NPs --> 408
	3.	On average, how many windows are present in an NP?
			AVG Windows Present --> 5,482.81 --> 5469.32
	4.	What is the smallest number of windows present in any NP? The largest? 
		  Min Windows Present --> 31 --> 29   Max Windows Present --> 21,249 --> 21,230
	5.	On average, what is the number of NPs in which a window is detected? The smallest? The largest?
			AVG NPs -->  24.56    MAX NPs --> 94    MIN NPs --> 0 
	6.	What is the number of GWs <= 100 NPs? --> 90852 #excluded 25 windows

  Radial Position - check
  Nuclear Organization
  NP similarity
1/28  the Jaccard Index  J(A, B) = |A INTERSECT B|
									-------------

"""

from operator import itemgetter
import constant
import csv
import time

start = time.time()

NP = []
GW = [] # len(GW) = 90,877    len(GW[0]) = 408
avg = [] #avg NP
avgG = [] # valid GWs
rank = [] # ranking NPs for Radial position
valid_GW = [] # track GW > 0 NPs and GW <= 100 NPs
nonvalid_GW = []


def tests():
	# print('# of Valid Genomic Windows:', len(avgG), '   Difference in Length = ', 90877 - len(avgG))
	# print('AVG Profiles detected in window:', round(sum(avgG) / len(avgG), 2) )
	# print('MIN Profiles detected in window:', min(avgG))
	# print('MAX Profiles detected in window:', max(avgG))

	# print('# of Invalid Genomic Windows: ', len(nonvalid_GW) )
	# print(*nonvalid_GW)
	
	print(*ranked_NP.items(), sep='\t')
	print('\n# of Ranked Nuclear Profiles', len(ranked_NP) )
	

	total = 0
	tempMin = 10000
	tempMax = 0
	for i in avg:
		total += i[1]

		if tempMin >= i[1]:
			tempMin = i[1]

		if tempMax <= i[1]:
			tempMax = i[1]	

	# print('\n')
	# print('# of Nuclear Profiles', len(avg))
	# print('avg Windows Present per NP:', round(total/len(avg), 2) )
	# print('min Windows Present in NP:', tempMin)
	# print('max Windows Present in NP:', tempMax)


filename = "input.csv"

with open(filename, 'r') as file:
	csvreader = csv.reader(file)
	next(csvreader) # skip NPs...F10A2, F10A3
	count = 0
	for row in csvreader: #read in the rows of data (Genomic Windows)
		row = list(map(int, row[3:]) )
		GW.append(row)

		if sum(row) <= 100: # filter the data: if a GW has more than 100 NPs, don't include it
			valid_GW.append(count)
			avgG.append(sum(row))
		else:
			nonvalid_GW.append(count)

		count += 1

file.close()

# disallow appendage to sum if the genomic window is not valid
countCol = 0
with open(filename, 'r') as file:
	csvreader = csv.reader(file)

	columns = list(zip(*csvreader)) #read in the columns of data (Nuclear Profiles)
	columns = list(columns)

	for column in columns[3:]: # chrom start stop are indecies 0 1 2
		columnName = column[0]
		column = list(map(int, column[1:]) ) # NP name excluded

		for i in nonvalid_GW:
			column[i] = 0
		
		avg.append((columnName, sum(column)) ) # column name ex. F10A2 && summation of 1s & 0s

	NP = columns[3:]

file.close()


# rank NPs based on their radial position
avg.sort(key=itemgetter(1)) # sort on 2nd index of each element

ranked_NP = {} # dictionary for storing ranked NPs
index = 0
rank = 1
for i in range(0, 408):
	if index == 82:
		rank += 1
		index = 0
	
	ranked_NP.update({ avg[i] : rank }) # NP name, # windows present, rank
	index += 1


tests()

end = time.time()
print('\nTime: ', round(end - start, 2), 'sec')
