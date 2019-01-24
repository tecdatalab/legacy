import glob
import os

#for file_name in glob.glob("/bio/kihara-scratch4/han273/allpdb_rerun_add/chain/*.pdb"):
	#command = 'cp '+file_name+' ./chain/'
	#print file_name
	#os.system(command)

#for file_name in glob.glob("/bio/kihara-scratch4/han273/allpdb_rerun_add/chain_cacn/*.pdb"):
	#command = 'cp '+file_name+' ./chain_cacn/'
	#print file_name
	#os.system(command)


for file_name in glob.glob("/bio/kihara-scratch4/han273/allpdb_rerun_add/inv/*.inv"):
	command = 'cp '+file_name+' ./inv/'
	print file_name
	os.system(command)

for file_name in glob.glob("/bio/kihara-scratch4/han273/allpdb_rerun_add/inv_cacn/*.inv"):
	command = 'cp '+file_name+' ./inv_cacn/'
	print file_name
	os.system(command)

for file_name in glob.glob("/bio/kihara-scratch4/han273/allpdb_rerun_add/anim/*.gif"):
	command = 'cp '+file_name+' ./anim/'
	print file_name
	os.system(command)

for file_name in glob.glob("/bio/kihara-scratch4/han273/allpdb_rerun_add/img/*.png"):
	command = 'cp '+file_name+' ./img/'
	print file_name
	os.system(command)