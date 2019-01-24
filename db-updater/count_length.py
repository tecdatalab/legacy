import glob
import subprocess
import os.path
import os

def removefile(name, directory):
	if os.path.exists(name):
		command = 'mv '+name+' ./short_seq/'+directory
		os.system(command)
		print 'successfully remove '+directory+' file' 

for filename in glob.glob('./chain_cacn/*.pdb'):
	command = 'grep ^ATOM.........CA '+filename+' | wc -l'
	length = subprocess.check_output(command, shell = True)

	if int(length) <= 10:
		print filename
		pdb_id = filename.split('chain_cacn/')[1].split('.pdb')[0]
		
		chain_name = './chain/'+pdb_id+'.pdb'
		removefile(chain_name,'chain')
		
		chain_cacn_name = filename
		removefile(chain_cacn_name,'chain_cacn')

		inv_name = './inv/'+pdb_id+'.pdb.inv'
		removefile(inv_name,'inv')

		inv_cacn_name = './inv_cacn/'+pdb_id+'.pdb.inv'
		removefile(inv_cacn_name,'inv_cacn')

		gif_name = './anim/'+pdb_id+'.gif'
		removefile(gif_name,'anim')
		
		png_name = './img/'+pdb_id+'.png'
		removefile(png_name,'img')

		print '\n'
