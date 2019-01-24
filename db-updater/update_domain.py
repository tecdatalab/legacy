# compare new domain boundary list with old one and existed domain chains

with open('./pdbidlistupdate/cath-domain-boundaries-seqreschopping-new.txt') as domain_new, open('./pdbidlistupdate/cath-domain-boundaries-seqreschopping.txt') as domain_old, open('./pdbidlistupdate/domainidlist.txt') as domain_id:
	domain_id_list = domain_id.read().split('\n')[:-1]
	domain_old_list = domain_old.read().split('\n')[:-1]
	domain_new_list = domain_new.read().split('\n')[:-1]

	#compare new domain boundary list with old one
	diff_list1 = list(set(domain_new_list) - set(domain_old_list)) 

	#extract domains that are not new but not processed before
	#read in domain ids in new domain boundary list
	new_domain_id_list = []
	for domain in domain_new_list:
		domain_id = domain.split()[0]
		if domain_id[5:] == '00':
			continue;

		pdb_id = domain_id[0:4]
		domain_chain_id = pdb_id + '-' + domain_id[4] + '-' + domain_id[5:]
		new_domain_id_list.append(domain_chain_id)

	#compare new domain id list with old one
	diff_list2 = list(set(new_domain_id_list) - set(domain_id_list))

	diff_list = diff_list1 + diff_list2
	new_pdb_list = [] #pdb ids of domains that need to be processed
	for element in diff_list:
		pdb_id = element[0:4]
		if not (pdb_id in new_pdb_list):
			new_pdb_list.append(pdb_id)

	output = open('./pdbidlistupdate/domain_new.pdb', 'w')
	output.write('\n'.join(new_pdb_list) + '\n')
	output.close()
			
