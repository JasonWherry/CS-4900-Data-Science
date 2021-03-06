""" 
Author: Jason Wherry	Date: 2/5/2020		To Run: python3 NuclSim.py

Nuclear Profile similarity: 
	- The Jaccard Index  J(A, B) = |A INTERSECT B| for each pair of relevant NPs
  	- compare a column of 80 rows with another column of 80 rows
  		- find their interesction to compute the similarity
	- store the computed values in a matrix
		- the matrix should contain similarities for all pairs of relevant NPs

Notes:
	- There are 162 Relevant NPs 


Output:
	sim Matrix length: 162 
	sim Matrix length of element: 162 

"""

from operator import itemgetter
import matplotlib.pyplot as plt
import seaborn as sb
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
nonvalid_GW = [] # values not in valid_GW
HIST1NP = []
seg = [] # for similarity comparison of NPs
simMatrix = []

def intersection(lst1, lst2): 
    return list(set(lst1) & set(lst2))

def tests():
	print('\n')
	print('# of Nuclear Profiles in HIST1:  ', len(HIST1NP) )
	print(' AVG Profiles detected in window:', round(sum(avgG[69716:69796]) / (80), 2) )
	print(' MIN Profiles detected in window:', min(avgG[69716:69796]))
	print(' MAX Profiles detected in window:', max(avgG[69716:69796]))

	total = 0
	tempMin = 10000
	tempMax = 0

# Traverse array HIST1: ('F10A3', 41), ('F10A5', 1), ('F10B5', 67), ...
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


	# print('\n', *HIST1NP, sep='\n\n')
	# print('Length of HIST1NP', len(HIST1NP))
	# print(*seg, sep='\n\n')
	print('\nLength of seg', len(seg))
	print('\nLength of seg[1][1]', len(seg[1][1]))

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
		seg.append(column)
		avg.append( (columnName, sum(column)) ) # column name ex. F10A2 && summation of 1s & 0s

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
	
	ranked_NP.update({ avg[i] : rank })
	index += 1

for item in avg:
	if item[1] > 0:
		HIST1NP.append(item) # list of 80 NPs - compare each element for Nuclear Org.

# deletes irrelevant NPs from seg list
for item in seg:
	if sum(item) > 0:
		pass
	else:
		for item in seg:
			if sum(item) > 0:
				pass
			else:
				seg.remove(item)

for item in seg:
	if sum(item) > 0:
		pass
	else:
		for item in seg:
			if sum(item) > 0:
				pass
			else:
				seg.remove(item)


# for item in seg:
# 	print(item, sep='\n\n')
# print('\n')
# print(len(seg)) # 162
# print(len(seg[1])) ## 80
# print(seg[1])


# two nested loops to compute jaccard values
denominator = 0
for i in range(0, 162):
	mMatrix = []
	
	for j in range(0, 162):
		numerator = 0
		denominator = min( sum(seg[i]), sum(seg[j]) )# + list(m01) + list(m10)

		for k in range(0, 80): # COLUMNS ARE TRAVERSED
			numerator += min(seg[i][k], seg[j][k]) # m11

		mMatrix.append(numerator/denominator)

	simMatrix.append(mMatrix) # simMatrix[i][j] = numerator/denominator


# print('\n\n')
# print('seg MATRIX:', *simMatrix, sep='\n\n')
# print('\nsim Matrix length', len(simMatrix), '\n')
# print('\nsim Matrix element', len(simMatrix[0]), '\n')

plt.imshow(simMatrix, cmap='hot', interpolation='nearest')
plt.title('Similarity Matrix')
plt.show()

distMatrix = []

for item in simMatrix:
	dMatrix = []
	for i in range(0, len(item)):
		dMatrix.append(1 - item[i])
	distMatrix.append(dMatrix)

plt.imshow(distMatrix, cmap='hot', interpolation='nearest')
plt.title('Distance Matrix')
plt.show()


#tests()

end = time.time()
print('\nTime:', round(end - start, 2), 'sec')