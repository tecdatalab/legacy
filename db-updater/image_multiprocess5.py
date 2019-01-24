# use 4 cpus to generate img and anim

import os 
from multiprocessing import Pool
 
def image(struc_id):

	chain_file = './chain/' + struc_id + '.pdb'

	if not (os.path.isfile('./anim/'+struc_id+'.gif')):

		cp_command = 'cp ' + chain_file + ' ./'
		os.system(cp_command)

		perl_command = 'perl gifandpng.pl ' + struc_id 
		os.system(perl_command)

		mv_command_1 = 'mv ' + struc_id + '.gif ./anim/' 
		mv_command_2 = 'mv ' + struc_id + '.png ./img/'
		os.system(mv_command_1)
		os.system(mv_command_2) 

		rm_command = 'rm ' + struc_id + '.pdb'
		os.system(rm_command)

		proc = os.getpid()
		print struc_id, proc
 
if __name__ == '__main__':
	group_list = ['domain'];
	
	for element in group_list:
		group_type = element
		added_idlist = "./pdbidlistupdate/added." + group_type + "idlist.txt"

		with open(added_idlist) as struc_file:
			struc_list = struc_file.read().split()

			pool = Pool(processes=25)
   	 		pool.map(image, struc_list)

	pool.close()
	pool.join()

	
