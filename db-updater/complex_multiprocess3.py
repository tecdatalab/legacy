# use 5 cpus to identify complexes

import os
from multiprocessing import Pool

def complex_distance(pdb_id):
	complex_output = './distance/' + pdb_id + '.distance'
	if not os.path.isfile(complex_output):
		perl_command = 'perl clustercomplexdistance.pl ' + pdb_id
		os.system(perl_command)

		proc = os.getpid()
		print pdb_id, proc

with open('/net/kihara/home/han273/test/complexpdblist_part3') as complex_list:
	pdb_list = complex_list.read().split()

	pool = Pool(processes=5)
	pool.map(complex_distance, pdb_list)
	pool.close()
	pool.join()
	
