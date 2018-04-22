import Bio.PDB as pdb
import gzip
import pandas as pd
import os
import json

from angles import angles_tetra
from utils import printProgressBar

parser = pdb.PDBParser(QUIET=True)
		
		
def process_pdb_file(filename, output_f):
	args = filename.split('.')
	
	if not os.path.isfile(filename):
		print('file %s not found' % filename)
		# raise ValueError
		return None
	
	if len(args) != 0 and args[-1] == 'gz':
		f = gzip.open(filename, 'r')
	else:
		# ? should close before
		f = open(filename, 'r')

	structure = parser.get_structure('1', f)

	for model in structure:
		for chain in model:

			tetraArr, anglesArr = angles_tetra(chain)
			
			if tetraArr is None: continue
			
			for i in range(len(tetraArr)):
				arr = [tetraArr[i]] + anglesArr[i]
				output_f.write(' '.join(map(str,arr)) + '\n')
			
def process_all():
	output_file = 'stat.txt'
	output_f = open(output_file, 'w')
	
	# get processed
	tab = pd.read_table('processed.txt', header=None, names=['id'])
	processed = {k : 1 for k in tab['id'].values}
	g = open('processed.txt', 'a')
	
	f = open('../../../chain_seq.json')
	s = f.read()
	dump = json.loads(s)
	f.close()
	
	keys = dump.keys()
	for i in range(len(keys)):
		pdbID = keys[i]
		
		if pdbID in processed.keys(): continue
		# pdbID = '4FF1'
		print(pdbID)
		# pdbID = '4BB9'
		# if pdbID in ['3SX7']: continue
		
		input = '../../../pdb/%s.pdb.gz' % pdbID
		process_pdb_file(input, output_f)
		
		g.write(pdbID + '\n')
		
		printProgressBar(i, len(keys), length=50, fill=chr(219))
		# break
	output_f.close()
	g.close()
	
def process_filtered():
	filtered = 'lables.txt'
	tab = pd.read_table(filtered, header=None, names=['pdbID', 'chainID'])
	tab.sort_values('pdbID', inplace=True)
	
if __name__ == '__main__':
	process_all()
	
	