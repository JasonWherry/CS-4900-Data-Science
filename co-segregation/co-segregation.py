"""
Author: Jason Wherry	Date: 3/30/2020		To Run: python3 co-segregation.py
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
	print(sumsWindows)

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
print('co_segregation_Table')
print(co_segregation_Table)
print(frequency)

def linkage(coseg_table):
	for i in range(0, 81):
		for j in range(0, 81):
			coseg_table.iloc[i, j] = coseg_table.iloc[i, j] - (frequency[i] * frequency[j])

	print(coseg_table)

	return coseg_table

linkage_table = linkage(co_segregation_Table)

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


heatmap = sb.heatmap(normalized_linkage_table, vmin = 0, vmax = 1)
plt.title('normalized_linkage_table')
plt.show()


