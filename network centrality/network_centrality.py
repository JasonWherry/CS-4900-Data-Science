"""
Author: Jason Wherry	Date: 4/06/2020		To Run: python3 network_centrality.py
"""

import numpy as np
import pandas as pd
from pandas import DataFrame
import csv
import seaborn as sb
import matplotlib.pyplot as plt

GW = 81
NP = 165

fileName = "HIST1.txt"
dataFrame = pd.read_excel(fileName, index_col=0)
data = np.array(dataFrame)

dataFrame2 = dataFrame.iloc[2:, 2:]
frequency = []

# extract relevant Nuclear Profiles (NPs which contain at least one window in the region of interest)
def extract_relevant_NPs(dataFrame):
	# collect sum of each NP
	sums = dataFrame.sum(axis = 0)
	sumsWindows = dataFrame.sum(axis = 1)
	# print(sumsWindows)

	for i, row in sumsWindows.items():
		frequency.append(row / 165)


	indices = np.where(sums == 0)
	deletions, *y = indices # unpack tuple

	dataFrame2.drop(dataFrame2.columns[deletions], axis=1, inplace = True) # drop irrelevant NPs from dataFrame by column

	return frequency

frequency = extract_relevant_NPs(dataFrame2)


def co_segregation(df):
	co_seg = []
	final = []

	df = df.T

	for i, row in df.items():
		# print(row)
		temp = []
		for j, val in row.items():
			# print(val)
			temp.append(val)
		co_seg.append(temp)

	# print(co_seg) # 2D array of NPs
	
	for j in range(0, len(co_seg)): # 0 to 81: windows
		a = co_seg[j] # window a
		finalTemp = []
		for d in range(0, len(co_seg)): # 0 to 81: NPs
			b = co_seg[d] # window b
			count = 0

			for z in range(0, len(co_seg[0])):
				if a[z] and b[z] == 1:
					count += 1

			finalTemp.append(count / 165)
		final.append(finalTemp)

	final = pd.DataFrame(final)
	
	return final

co_segregation_Table = co_segregation(dataFrame2)
# print('\n co_segregation_table')
# print(co_segregation_Table)

# print(frequency)

def linkage(coseg_table):
	for i in range(0, 81):
		for j in range(0, 81):
			coseg_table.iloc[i, j] = coseg_table.iloc[i, j] - (frequency[i] * frequency[j])

	# print(coseg_table)

	return coseg_table

linkage_table = linkage(co_segregation_Table)
# print('\n linkage_table')
# print(linkage_table)

def normalized_linkage(linkage_table):
	for i in range(0, 81):
		for j in range(0, 81):
			if linkage_table.iloc[i,j] < 0:
				linkage_table.iloc[i, j] = linkage_table.iloc[i, j] / min( (frequency[i] * frequency[j]), (1 - frequency[i]) * (1 - frequency[j]) )
			elif linkage_table.iloc[i, j] > 0:
				linkage_table.iloc[i, j] = linkage_table.iloc[i, j] / min( (frequency[j] * (1 - frequency[i]) ), (frequency[i] * (1 - frequency[j])) )

	return linkage_table

	
normalized_linkage_table = normalized_linkage(linkage_table)
print('\n normalized_linkage_table')
print(normalized_linkage_table)
print('\n')
""" 
heatmap = sb.heatmap(normalized_linkage_table, vmin = 0, vmax = 1)
plt.title('normalized_linkage_table')
plt.show()
"""


""" ------------------------------------ ------------------------------------ -----------------------------------
 L_avg(DataFrame)
	- find L-avg of the normalized_linkage_table
	- parameter: 	pandas DataFrame
	- return:		pandas avg (includes diagonal), my avg (ignores diagonal)
"""
def L_avg(normalized_linkage_table):
	summation = 0
	average = normalized_linkage_table.mean(axis = 0)

	for i in range(0, 81):
		if average.iloc[i] != 1:
			summation += average.iloc[i]

	summation = summation / 80			# ignore reflexive edges
	average = average.mean(axis = 0)	# Does not ignore reflexive edges

	return average, summation

pandas_L_avg, my_L_avg = L_avg(normalized_linkage_table)
print('pandas L-avg:\t', pandas_L_avg, '\nmy L-avg:\t', my_L_avg)
print('\n')

"""
 determine_edges(normalized_linkage_table, L_avg)
 	- determine edges of the network
	- parameter types: 	pandas DataFrame, numpy.float64
	- return:			adjacency matrix (Edges)
"""
def determine_edges(normalized_linkage_table, L_avg):
	temp_data_frame = normalized_linkage_table
	for i in range(0, 81):
		for j in range(0, 81):
			if i == j:
				temp_data_frame.iloc[i,j] = 0 # mark diagonal as 0

			elif normalized_linkage_table.iloc[i,j] > L_avg:
				temp_data_frame.iloc[i,j] = 1

			elif normalized_linkage_table.iloc[i, j] <= L_avg:
				temp_data_frame.iloc[i,j] = 0

	return temp_data_frame

adjacency_matrix = determine_edges(normalized_linkage_table, my_L_avg)
print('adjacency_matrix (Edges)')
print(adjacency_matrix)
print('\n')


"""
 degree_centrality(adjacency_matrix)
 	- calculates the degree centrality of a window
	- parameter type: 	pandas DataFrame (81 x 81)
	- return:			pandas series
"""
def degree_centrality(adjacency_matrix):
	summation = adjacency_matrix.sum(axis=1) # sum each row
	for i in range(0, 81):
		summation[i] = summation[i] / 80

	return(summation)

degree_centrality_matrix = degree_centrality(adjacency_matrix)
# print('degree centrality\n', degree_centrality_matrix)
# print('\n')


"""
 degree_centrality_stats(degree_centrality_matrix, my_L_avg)
 	- finds avg, min and max degree centrality
 	- rank the windows in ascending order based on their degree centrality
	- parameter type: 	pandas series, numpy.float64
	- return:			Average, min and max degree centrality
"""
def degree_centrality_stats(degree_centrality_matrix, my_L_avg):
	print('L-avg:', my_L_avg, '\n')
	# find some stats of degree centrality
	average = degree_centrality_matrix.mean()
	minimum = degree_centrality_matrix.min()
	maximum = degree_centrality_matrix.max()

	# rank windows (ascending order)
	ranked_windows = degree_centrality_matrix.sort_values(axis=0, ascending=True)

	# a ranked list of windows and their centrality values (ascending order) 
	counter = 1
	print('rank\t', 'window\t', 'degree centrality')
	for i, j in ranked_windows.items():
		print(counter, '\t', i, '\t', j)
		counter += 1

	return average, minimum, maximum

dc_avg, dc_min, dc_max = degree_centrality_stats(degree_centrality_matrix, my_L_avg)
print('\n\n')
print('Average, min, and max degree centrality')
print('avg:', dc_avg, '\tmin:', dc_min, '\tmax:', dc_max)
print('\n')


