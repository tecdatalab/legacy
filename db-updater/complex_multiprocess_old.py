# use 5 cpus to identify complexes

import os
from multiprocessing import Pool

def complex_distance(pdb_id):
	perl_command = 'perl clustercomplexdistance.pl ' + pdb_id
	os.system(perl_command)

	proc = os.getpid()
	print pdb_id, proc

with open('./pdbidlistupdate/complexpdblist.txt') as complex_list:
	pdb_list = complex_list.read().split()

	pool = Pool(processes=5)
	pool.map(complex_distance, pdb_list)
	pool.close()
	pool.join()
	
