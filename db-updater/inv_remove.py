import glob
import os
import os.path

def removefile(name, directory):
	if os.path.exists(name):
		command = 'mv '+name+' ./short_seq/'+directory
		os.system(command)
		print 'successfully remove '+directory+' file'

for filename in glob.glob('./short_seq/chain/*.pdb'):
	print filename
	pdb_id = filename.split('short_seq/chain/')[1].split('.pdb')[0]

	inv_name = './inv/'+pdb_id+'.pdb.inv'
	removefile(inv_name,'inv')

	inv_cacn_name = './inv_cacn/'+pdb_id+'.pdb.inv'
	removefile(inv_cacn_name,'inv_cacn')

	print '\n'
