""" 
Author: Jason Wherry	Date: 2/5/2020		To Run: python3 NuclOrg.py

Output:
	# of Nuclear Profiles in HIST1:  162
	AVG Profiles detected in window: 30.76
	MIN Profiles detected in window: 0
	MAX Profiles detected in window: 42

	# of Valid Genomic Windows: 80
	avg Windows Present per NP: 15.0
	min Windows Present in NP: 1
	max Windows Present in NP: 68

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
	print('\n')
	print('# of Nuclear Profiles in HIST1:  ', len(HIST1NP) )
	print(' AVG Profiles detected in window:', round(sum(avgG[69716:69796]) / (80), 2) )
	print(' MIN Profiles detected in window:', min(avgG[69716:69796]))
	print(' MAX Profiles detected in window:', max(avgG[69716:69796]))

	total = 0
	tempMin = 10000
	tempMax = 0

# Traverse dictionary HIST1: { ('F10A3', 41), ('F10A5', 1), ('F10B5', 67), ...}
	for i in HIST1NP: 
		total += i[1]

		if i[1] <= tempMin:
			tempMin = i[1]

		if i[1] >= tempMax:
			tempMax = i[1]	

	print('\n# of Valid Genomic Windows:', len(avgG[69716:69796]) )
	print('avg Windows Present per NP:', round(total/len(HIST1NP), 2) )
	print('min Windows Present in NP:', tempMin)
	print('max Windows Present in NP:', tempMax)


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
		column = list(map(int, column[69716:69796]) ) # NP name excluded
		avg.append((columnName, sum(column)) ) # column name ex. F10A2 && summation of 1s & 0s

	NP = columns[3:]

file.close()

HIST1NP = []

for item in avg:
	if item[1] > 0:
		HIST1NP.append(item) # list of 80 NPs - compare each element for Nuclear Org.

# rank NPs based on their radial position
avg.sort(key=itemgetter(1)) # sort on 2nd index of each element

ranked_NP = {} # dictionary for storing ranked NPs
index = 0
rank = 1
for i in range(0, 408):
	if index == 82:
		rank += 1
		index = 0
	
	ranked_NP.update({ avg[i] : rank })
	index += 1


tests()


end = time.time()
print('\nTime:', round(end - start, 2), 'sec')