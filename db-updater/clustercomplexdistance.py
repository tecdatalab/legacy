import sys
import numpy as np
from itertools import combinations

def calculate(data, residue, chain):
	origin = atomInfo[np.intersect1d(np.where(atomInfo[:, 16] == residue), np.where(atomInfo[:, 18] == chain))]
	target = atomInfo[np.where(atomInfo[:, 18] != chain)]
	target_list = []
	for each_chain in chainIdArray: 
		sub_list = target[np.where(target[:, 18] == each_chain)]
		if len(sub_list) != 0: 
			target_list.append(sub_list)
	for each_origin in origin: 
		flag = False
		for each_target_list in target_list: 
			for each_target in each_target_list: 
				distance = (float(each_origin[10]) - float(each_target[10])) ** 2 + (float(each_origin[11]) - float(each_target[11])) ** 2 + (float(each_origin[12]) - float(each_target[12])) ** 2
				if (distance < 4.5 * 4.5): 
					if not (chain, each_target[18]) in result: 
						result[chain, each_target[18]] = 1
					else: 
						result[chain, each_target[18]] += 1
					flag = True
					break
			if flag: break
		if flag: break
	return

filename = sys.argv[1]
with open('./chain/' + filename + '.cif') as f: 
	content = f.readlines()

atomInfo = []

for line in content: 
	if line.split()[0] == 'ATOM': 
		atomInfo.append(line.split())
atomInfo = np.array(atomInfo)

chainIdArray = np.unique(atomInfo[:, 18])
result = {}

if len(chainIdArray) < 2: 
	exit()

for each_chain in chainIdArray: 
	residueIndex = np.unique(np.array(atomInfo[np.where(atomInfo[:, 18] == each_chain)])[:, 16])
	for each_index in residueIndex: 
		print(each_index)
		calculate(atomInfo, each_index, each_chain)

f = open("distance/" + filename + ".distance", "w+")
_buffer = ""
for each_combination in list(combinations(chainIdArray, 2)): 
	#print(each_combination[0], end = '\t')
	_buffer = _buffer + str(each_combination[0]) + '\t'
	if each_combination in result: 
		#print(result[each_combination], end = '\t')
		_buffer = _buffer + str(result[each_combination]) + '\t'
	else: 
		#print('0', end = '\t')
		_buffer = _buffer + '0' + '\t'
	#print(each_combination[1], end = '\t')
	_buffer = _buffer + str(each_combination[1]) + '\t'
	if each_combination[::-1] in result: 
		#print(result[each_combination[::-1]])
		_buffer = _buffer + str(result[each_combination[::-1]]) + '\n'
	else: 
		#print('0')
		_buffer = _buffer + '0\n'
f.write(_buffer)
f.close()
