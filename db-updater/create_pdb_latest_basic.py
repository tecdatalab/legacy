import sys
import numpy as np

path = sys.argv[1]
pdb_list_file = sys.argv[2]
output_filename = sys.argv[3]

protein_list = np.array(open('./pdbidlistupdate/added.protandpna.pdb','r').read().split('\n'))
source_list = np.array(open(pdb_list_file, 'r').read().split('\n'))
destination_list = np.array(open('./pdbidlistupdate/pdb_latest_' + pdb_list_file[18:-10], 'r').read().split('\n'))
zernike_count = destination_list[0].split()[1]

for i in range(2, len(destination_list) - 1): 
	space_location = [j for j, c in enumerate(destination_list[i]) if c ==' ']
	destination_list[i] = destination_list[i][0:space_location[-5]]

for item in protein_list: 
	source_match = [i for i, entries in enumerate(source_list) if item in entries]

for index in source_match: 
	destination_match = [i for i, entries in enumerate(destination_list) if source_list[index] in entries]
	data = source_list[index] + '.cif.inv '
	with open(path + '/' + source_list[index] + '.cif.inv', 'r') as file: 
		data = data + file.read().replace('\n', ' ')
		data = data[:-1]
	if (len(destination_match) != 0): 
		destination_list[destination_match[0]] = data
	else: 
		zernike_count = int(zernike_count) + 1
		destination_list = np.append(destination_list, data)

destination_list[0] = 'ZERNIKE ' + str(zernike_count)

f = open(output_filename, 'w+')
for i in range(0, len(destination_list)): 
	f.write(destination_list[i] + '\n')
f.close()
