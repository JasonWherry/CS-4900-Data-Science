""" Author: Jason Wherry	Date: 2/18/2020		To Run: python3 clusterPart21.py

Clustering Part 2:
	Tasks:
		- For each cluster, find the NP whose average dissimilarity to all the objects in the cluster is minimal.
		- These are the centers of the new clusters.
			- Found the center of each cluster
			- Create a dictionary for each cluster in order to retrieve the average NP's name 

		- Calculate the mean of each cluster
			- stored in the variable 'AVGmeanNP'
		- Repeat until clusters no longer change
		- Assess the quality of the clustering by adding up the variation within each cluster

	Notes:
		- I have stored the NPs by cluster (1, 2, and 3)
		- For every cluster, the NP that is the minimum difference from the AVG NP Objects is what I seek

	Output:
			clust1 length 76	clust2 length 31	clust3 length 55
			clust1 + clust2 + clust3 =  162		
			clust1 length 80	
			Time: 18.5 sec
				
"""

from operator import itemgetter
from statistics import mean
import matplotlib.pyplot as plt
import seaborn as sb
import constant
import csv
import time
import random
from statistics import variance

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
diffMatrix = []
clusters = []

def intersection(lst1, lst2): 
    return list(set(lst1) & set(lst2))

def plotClusters(clusters, clust1, clust10, clust2, clust3):
	plt.imshow(clusters, cmap='hot', interpolation='nearest')
	plt.title('Clusters')
	plt.ylabel('NPs')
	plt.xlabel('Windows')
	plt.show()

	plt.imshow(clust1, cmap='hot', interpolation='nearest')
	plt.title('Cluster 1')
	plt.ylabel('NPs')
	plt.xlabel('Windows')
	plt.show()

	plt.imshow(clust2, cmap='hot', interpolation='nearest')
	plt.title('Cluster 2')
	plt.ylabel('NPs')
	plt.xlabel('Windows')
	plt.show()

	plt.imshow(clust3, cmap='hot', interpolation='nearest')
	plt.title('Cluster 3')
	plt.ylabel('NPs')
	plt.xlabel('Windows')
	plt.show()


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
			if sum(item) > 0: # changed from 0 to 10
				pass
			else:
				seg.remove(item)

for item in seg:
	if sum(item) > 0:
		pass
	else:
		for item in seg:
			if sum(item) > 0: # changed from 0 to 10
				pass
			else:
				seg.remove(item)

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

# display similarity matrix
plt.imshow(simMatrix, cmap='hot', interpolation='nearest')
plt.title('Similarity Matrix')
# plt.show()

distMatrix = []

for item in simMatrix:
	dMatrix = []
	for i in range(0, len(item)):
		dMatrix.append(1 - item[i])
	distMatrix.append(dMatrix)

# display distance matrix
plt.imshow(distMatrix, cmap='hot', interpolation='nearest')
plt.title('Distance Matrix')
# plt.show()

for j in range(0, 3):
	r1 = []
	for i in range(0, 80):
		r1.append(random.randint(0, 1))
	seg.append(r1)	

# measure the distance between each point and each of the k clusters
# two nested loops to compute jaccard values
denominator = 0
for i in range(162, 165):
	mMatrix = []
	
	for j in range(0, 162):
		numerator = 0
		denominator = min( sum(seg[i]), sum(seg[j]) )# + list(m01) + list(m10)

		for k in range(0, 80): # COLUMNS ARE TRAVERSED
			numerator += min(seg[i][k], seg[j][k]) # m11

		mMatrix.append(numerator/denominator)

	diffMatrix.append(mMatrix) # simMatrix[i][j] = numerator/denominator

# Assign each relevant NP (point) to the nearest cluster
for i in range(0, 162):
	val1 = diffMatrix[0][i]
	val2 = diffMatrix[1][i]
	val3 = diffMatrix[2][i]

	if min(val1, val2, val3) == val1:
		clusters.append(1) # append k-value

	elif min(val1, val2, val3) == val2:
		clusters.append(2) # append k-value

	elif min(val1, val2, val3) == val3:
		clusters.append(3) # append k-value


# for i in range(0, len(clusters)):
# 	print(clusters[i])

print('length', len(clusters) )

# 3 clusters and 3 cluster dictionaries
clust1 = []; clust10 = {}; clust2 = [] ; clust20 = {}; clust3 = []; clust30 = {}

# loop through seg matrix: 162 x 80
# clusters is 162 x 1
for i in range(0, 162):
	if clusters[i] == 1:
		clust1.append(seg[i])
		clust10.update({avg[i][0]:seg[i]})

	elif clusters[i] == 2:
		clust2.append(seg[i])
		clust20.update({avg[i][0]:seg[i]})

	elif clusters[i] == 3:
		clust3.append(seg[i])
		clust30.update({avg[i][0]:seg[i]})	


# print('\nClust1', clust1)
# print('\n\nClust2', clust2)
# print('\n\nClust3', clust3)
# print('\n\nClust3 has', len(clust3), 'elements')
# print('\n\nClust3 element length is', len(clust3[0]) )

# print('\nCLUSTERS')
# print('\n\tClust1 + Clust2 + Clust3 = ',len(clust1),'+',len(clust2),'+',len(clust3),'=>',len(clust1)+len(clust2)+len(clust3))
# print('\tLength of each Cluster object (NP) = ', len(clust3[0]),'\n\n' )

# Finds the number in a list closest to variable k
def closest(lst, K): 
    return lst[min(range(len(lst)), key = lambda i: abs(lst[i]-K))] 

# determine the mediod of each cluster
def findMinAvgDiff(cluster, clustDict, clustNum):
	# Find the center of each cluster
	meanNP = map(mean, zip(*cluster))
	# print('\nmeanNP', *meanNP, sep='\n')

	# list of the means for every NP in a cluster
	meanNPList = list(meanNP)
	# print(*meanNPList, sep='\n')

	# take the NP mean for every NP in a cluster and find the overall average mean
	AVGmeanNP = sum(meanNPList) / len(meanNPList)
	indexOfNP = meanNPList.index(closest(meanNPList, AVGmeanNP))

	# print('index of NP', indexOfNP, '\n')

	keyCount = 0
	chosenOne = 0 
	for key, vaue in clustDict.items():
		if keyCount == indexOfNP:
			chosenOne = key # find the name of the NP
		keyCount += 1

	print('CLUSTER', clustNum)
	print('Length of cluster', len(cluster), '  length of dictionary', len(clustDict))
	print('\n\tAVG of the cluster means is', AVGmeanNP)
	print('\tClosest Value is',closest(meanNPList, AVGmeanNP), ' is at index', indexOfNP ) 
	print('\tThe center of cluster', clustNum, 'is', chosenOne, 'at index', indexOfNP)
	print('\n')

	chosenOne = 0

	return(AVGmeanNP)


def calcVariance(clust, indexOfCenter) :
	count = 0
	distList = []
	for sublist in clust:

		if count == indexOfCenter:
			pass

		else:
			distList.append( abs( (sum(clust[indexOfCenter])/len(clust[indexOfCenter])) - (sum(sublist)/len(sublist)) ) )

		count += 1

	# print(distList)

	return variance(distList)

# clustersAVG is an list of 3 items -> clust1, 2, & 3's closest value found
def clusterize(clustersAVG, seg) :
	diffMatrix = [] # previously global variables
	clusters = [] # previously global variables

	# measure the distance between each point and each of the k clusters
	denominator = 0
	for i in range(0, 3):
		mMatrix = []
		
		for j in range(0, 162): # two nested loops to compute jaccard values
			numerator = 0
			denominator = min( clustersAVG[i], sum(seg[j]) )# + list(m01) + list(m10)

			for k in range(0, 80): # COLUMNS ARE TRAVERSED
				numerator += min(seg[i][k], seg[j][k]) # m11

			mMatrix.append(numerator/denominator)

		diffMatrix.append(mMatrix) # simMatrix[i][j] = numerator/denominator

	# Assign each relevant NP (point) to the nearest cluster
	for i in range(0, 162):
		val1 = diffMatrix[0][i]
		val2 = diffMatrix[1][i]
		val3 = diffMatrix[2][i]

		if min(val1, val2, val3) == val1:
			clusters.append(1) # append k-value

		elif min(val1, val2, val3) == val2:
			clusters.append(2) # append k-value

		elif min(val1, val2, val3) == val3:
			clusters.append(3) # append k-value

	# 3 clusters and 3 cluster dictionaries
	clust1 = []; clust10 = {}; clust2 = [] ; clust20 = {}; clust3 = []; clust30 = {}

	# loop through seg matrix: 162 x 80
	# clusters is 162 x 1
	for i in range(0, 162):
		if clusters[i] == 1:
			clust1.append(seg[i])
			clust10.update({avg[i][0]:seg[i]})

		elif clusters[i] == 2:
			clust2.append(seg[i])
			clust20.update({avg[i][0]:seg[i]})

		elif clusters[i] == 3:
			clust3.append(seg[i])
			clust30.update({avg[i][0]:seg[i]})	


	# print('\nClust1', clust1)
	# print('\n\nClust2', clust2)
	# print('\n\nClust3', clust3)
	# print('\n\nClust3 has', len(clust3), 'elements')
	# print('\n\nClust3 element length is', len(clust3[0]) )

	# print('\nCLUSTERS')
	# print('\n\tClust1 + Clust2 + Clust3 = ',len(clust1),'+',len(clust2),'+',len(clust3),'=>',len(clust1)+len(clust2)+len(clust3))
	# print('\tLength of each Cluster object (NP) = ', len(clust3[0]),'\n\n' )

	tempMeanList = []
	tempMeanList.append( findMinAvgDiff(clust1, clust10, 1) )
	tempMeanList.append( findMinAvgDiff(clust2, clust20, 2) )
	tempMeanList.append( findMinAvgDiff(clust3, clust30, 3) )


	clusters = clust1 + clust2 + clust3

	# print(clusters, sep='\n\n')
	# print('\nlength of clusters', len(clusters))
	# print('\nlength of cluster element', len(clusters[0]))
	# print('\nlength of clust 1', len(clust1))
	# print('\nlength of clust 1', len(clust1[0]))

	# plotClusters(clusters, clust1, clust10, clust2, clust3)

	# print('\nVariance clust1', variance(clust1), '\n')
	# print('clust1', clust1)
	# print('clust1[0] length', len(clust1[0]))

	# print('Clust1 Variance', calcVariance(clust1, 0))
	# print('Clust2 Variance', calcVariance(clust2, 6))
	# print('Clust3 Variance', calcVariance(clust3, 9))

	file = open('clust1.txt', 'a')
	file2 = open('clust1Names.txt', 'a')
	file3 = open('clust1Ranks.txt', 'a')
	# print cluster dictionaries
	

	for key, val in clust10.items():

		file.write(str(val)) # change list to str to write to file
		file.write(',')
		file2.write(key)
		file2.write(',')
		
	file.close()
	file2.close()
	file3.close()

	file = open('clust2.txt', 'a')
	file2 = open('clust2Names.txt', 'a')
	# print cluster dictionaries
	for key, val in clust20.items():
		file.write(str(val)) # change list to str to write to file
		file.write(',')
		file2.write(key)
		file2.write(',')

	file.close()
	file2.close()


	file = open('clust3.txt', 'a')
	file2 = open('clust3Names.txt', 'a')
	# print cluster dictionaries
	for key, val in clust30.items():
		# print(val)
		file = open('clust3.txt', 'a') # 'a' used to append to file
		# val2 = val + '\n'
		file.write(str(val)) # change list to str to write to file
		file.write(',')
		file2.write(key)
		file2.write(',')
		
	file.close()
	file2.close()

	print(clust20)
	print()
	print(ranked_NP)

	return tempMeanList


tempAVGList = []
tempAVGList.append( findMinAvgDiff(clust1, clust10, 1) )
tempAVGList.append( findMinAvgDiff(clust2, clust20, 2) )
tempAVGList.append( findMinAvgDiff(clust3, clust30, 3) )

clusters = clust1 + clust2 + clust3

# print(clusters, sep='\n\n')
# print('\nlength of clusters', len(clusters))
# print('\nlength of cluster element', len(clusters[0]))
# print('\nlength of clust 1', len(clust1))
# print('\nlength of clust 1', len(clust1[0]))

# plotClusters(clusters, clust1, clust10, clust2, clust3)

# print('Clust1 Variance', calcVariance(clust1, 10))
# print('Clust2 Variance', calcVariance(clust2, 3))
# print('Clust3 Variance', calcVariance(clust3, 15))

iteration3 = clusterize(tempAVGList, seg)
# iteration4 = clusterize(iteration3, seg)
# iteration5 = clusterize(iteration4, seg)
# clusterize(iteration5, seg)


end = time.time()
print('\nTime:', round(end - start, 2), 'sec')