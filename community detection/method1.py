"""
Author: Jason Wherry	Date: 4/18/2020		To Run: python3 method1.py
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
# print('\n normalized_linkage_table')
# print(normalized_linkage_table)
# nlt = normalized_linkage_table.copy()
# print('\n')
 
# heatmap = sb.heatmap(normalized_linkage_table, vmin = 0, vmax = 1)
# plt.title('normalized_linkage_table')
# plt.show()



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
# print('pandas L-avg:\t', pandas_L_avg, '\nmy L-avg:\t', my_L_avg)
# print('\n')

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
# print('adjacency_matrix (Edges)')
# print(adjacency_matrix)
# print('\n')

edge_copy = adjacency_matrix.copy()

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
	# print('L-avg:', my_L_avg, '\n')
	# find some stats of degree centrality
	average = degree_centrality_matrix.mean()
	minimum = degree_centrality_matrix.min()
	maximum = degree_centrality_matrix.max()

	# rank windows (ascending order)
	ranked_windows = degree_centrality_matrix.sort_values(axis=0, ascending=True)

	# a ranked list of windows and their centrality values (ascending order) 
	# counter = 1
	# print('rank\t', 'window\t', 'degree centrality')
	# for i, j in ranked_windows.items():
	# 	print(counter, '\t', i, '\t', j)
	# 	counter += 1

	return average, minimum, maximum, ranked_windows

dc_avg, dc_min, dc_max, ranked_windows = degree_centrality_stats(degree_centrality_matrix, my_L_avg)
# print('\n\n')
# print('Average, min, and max degree centrality')
# print('avg:', dc_avg, '\tmin:', dc_min, '\tmax:', dc_max)
# print('\n')

""" ------------------------------------ ------------------------------------ -----------------------------------"""
LAD =   [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
Hist1 = [0,0,0,1,1,1,1,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,1,1,1,1,0,0,0,0,1,0,0,0]

"""
 community_hubs(ranked_windows, edges)
 	- finds the five nodes with the largest degree centrality
	- parameter type: 	pandas series, pandas DataFrame (copy edge matrix)
	- return:			five pandas series (hubs)
"""
def community_hubs(ranked_windows, edges):
	print('\n')
	print('5 Nodes with the Highest Degree Centrality ')
	print('window 53', ranked_windows.loc[53])
	print('window 34', ranked_windows.loc[34])
	print('window 29', ranked_windows.loc[29])
	print('window 36', ranked_windows.loc[36])
	print('window 16', ranked_windows.loc[16])

	hub1 = edges.iloc[53,:]
	hub2 = edges.iloc[34,:]
	hub3 = edges.iloc[29,:]
	hub4 = edges.iloc[36,:]
	hub5 = edges.iloc[16,:]

	return(hub1, hub2, hub3, hub4, hub5)

hub1, hub2, hub3, hub4, hub5 = community_hubs(ranked_windows, edge_copy)

"""
 develop_communities(edges, hub1, hub2, hub3, hub4, hub5)
 	- finds the neighbors for each community
	- parameter type: 	pandas DataFrame, pandas series
	- return:			none
"""
def develop_communities(edges, hub1, hub2, hub3, hub4, hub5):
	comm1 = pd.DataFrame(index=range(81),columns=range(81))
	comm1 = comm1[comm1.columns].astype(float)
	comm2 = pd.DataFrame(index=range(81),columns=range(81))
	comm2 = comm2[comm2.columns].astype(float)
	comm3 = pd.DataFrame(index=range(81),columns=range(81))
	comm3 = comm3[comm3.columns].astype(float)
	comm4 = pd.DataFrame(index=range(81),columns=range(81))
	comm4 = comm4[comm4.columns].astype(float)
	comm5 = pd.DataFrame(index=range(81),columns=range(81))
	comm5 = comm5[comm5.columns].astype(float)

	for i in range(0, 81):
		curr_row = edge_copy.iloc[i, :]
		counts = []
		keys = []
		count1 = 0
		count2 = 0
		count3 = 0
		count4 = 0
		count5 = 0
		for j in range(0, 81):
			if curr_row[j] == hub1[j]:
				count1 += 1
			if curr_row[j] == hub2[j]:
				count2 += 1
			if curr_row[j] == hub3[j]:
				count3 += 1
			if curr_row[j] == hub4[j]:
				count4 += 1
			if curr_row[j] == hub5[j]:
				count5 += 1
		counts.append(count1)
		counts.append(count2)
		counts.append(count3)
		counts.append(count4)
		counts.append(count5)

		max_val = max(counts)	# max value for a node	

		for n in range(0, 5):
			if counts[n] == max_val:
				keys.append(n)

		# community 1
		if 0 in keys:
			for k in range(0, 81):
				comm1.iloc[i, k] = 1
		else:
			for k in range(0, 81):
				comm1.iloc[i, k] = 0

		# community 2
		if 1 in keys:
			for k in range(0, 81):
				comm2.iloc[i, k] = 1
		else:
			for k in range(0, 81):
				comm2.iloc[i, k] = 0

		# community 3
		if 2 in keys:
			for k in range(0, 81):
				comm3.iloc[i, k] = 1
		else:
			for k in range(0, 81):
				comm3.iloc[i, k] = 0

		# community 4
		if 3 in keys:
			for k in range(0, 81):
				comm4.iloc[i, k] = 1
		else:
			for k in range(0, 81):
				comm4.iloc[i, k] = 0

		# community 5
		if 4 in keys:
			for k in range(0, 81):
				comm5.iloc[i, k] = 1
		else:
			for k in range(0, 81):
				comm5.iloc[i, k] = 0

	# find the nodes (windows with 1 values) in the community 1
	nodes1 = []
	LAD_gene = 0
	Hist1_gene = 0
	for n in range(0, 81):
		if comm1.iloc[n, 0] == 1:
			nodes1.append(n)
		if LAD[n] == comm1.iloc[n, 0]:
			LAD_gene += 1
		if Hist1[n] == comm1.iloc[n, 0]:
			Hist1_gene += 1
	print('\nnodes in community 1, hub = window 53: ', len(nodes1), 'total')
	print(nodes1)
	print('\n')
	print('% of nodes in community 1 that contain a LAD:', (LAD_gene/81)*100, '%')
	print('% of nodes in community 1 that contain a Hist1:', (Hist1_gene/81)*100, '%')
	print('\n')

	# community 2
	nodes2 = []
	LAD_gene2 = 0
	Hist1_gene2 = 0
	for n in range(0, 81):
		if comm2.iloc[n, 0] == 1:
			nodes2.append(n)
		if LAD[n] == comm2.iloc[n, 0]:
			LAD_gene2 += 1
		if Hist1[n] == comm2.iloc[n, 0]:
			Hist1_gene2 += 1
	print('\nnodes in community 2, hub = window 34: ', len(nodes2), 'total')
	print(nodes2)
	print('\n')
	print('% of nodes in community 2 that contain a LAD:', (LAD_gene2/81)*100, '%')
	print('% of nodes in community 2 that contain a Hist1:', (Hist1_gene2/81)*100, '%')
	print('\n')

	# community 3
	nodes3 = []
	LAD_gene3 = 0
	Hist1_gene3 = 0
	for n in range(0, 81):
		if comm3.iloc[n, 0] == 1:
			nodes3.append(n)
		if LAD[n] == comm3.iloc[n, 0]:
			LAD_gene3 += 1
		if Hist1[n] == comm3.iloc[n, 0]:
			Hist1_gene3 += 1
	print('\nnodes in community 3, hub = window 29: ', len(nodes3), 'total')
	print(nodes3)
	print('\n')
	print('% of nodes in community 3 that contain a LAD:', (LAD_gene3/81)*100, '%')
	print('% of nodes in community 3 that contain a Hist1:', (Hist1_gene3/81)*100, '%')
	print('\n')

	# community 4
	nodes4 = []
	LAD_gene4 = 0
	Hist1_gene4 = 0
	for n in range(0, 81):
		if comm4.iloc[n, 0] == 1:
			nodes4.append(n)
		if LAD[n] == comm4.iloc[n, 0]:
			LAD_gene4 += 1
		if Hist1[n] == comm4.iloc[n, 0]:
			Hist1_gene4 += 1
	print('\nnodes in community 4, hub = window 36: ', len(nodes4), 'total')
	print(nodes4)
	print('\n')
	print('% of nodes in community 4 that contain a LAD:', (LAD_gene4/81)*100, '%')
	print('% of nodes in community 4 that contain a Hist1:', (Hist1_gene4/81)*100, '%')
	print('\n')

	# community 5
	nodes5 = []
	LAD_gene5 = 0
	Hist1_gene5 = 0
	for n in range(0, 81):
		if comm5.iloc[n, 0] == 1:
			nodes5.append(n)
		if LAD[n] == comm5.iloc[n, 0]:
			LAD_gene5 += 1
		if Hist1[n] == comm5.iloc[n, 0]:
			Hist1_gene5 += 1
	print('\nnodes in community 5, hub = window 16: ', len(nodes5), 'total')
	print(nodes5)
	print('\n')
	print('% of nodes in community 5 that contain a LAD:', (LAD_gene5/81)*100, '%')
	print('% of nodes in community 5 that contain a Hist1:', (Hist1_gene5/81)*100, '%')
	print('\n')
	

	# visualize community with heatmap 
	title = 'community1: hub is window 53' 
	heatmap = sb.heatmap(comm1, vmin = 0, vmax = 1)
	plt.title(title)
	plt.show()

	title = 'community2: hub is window 34' 
	heatmap = sb.heatmap(comm2, vmin = 0, vmax = 1)
	plt.title(title)
	plt.show()

	title = 'community3: hub is window 29' 
	heatmap = sb.heatmap(comm3, vmin = 0, vmax = 1)
	plt.title(title)
	plt.show()

	title = 'community4: hub is window 36' 
	heatmap = sb.heatmap(comm4, vmin = 0, vmax = 1)
	plt.title(title)
	plt.show()

	title = 'community5: hub is window 16' 
	heatmap = sb.heatmap(comm5, vmin = 0, vmax = 1)
	plt.title(title)
	plt.show()

	# return temp_matrix

develop_communities(edge_copy, hub1, hub2, hub3, hub4, hub5)


