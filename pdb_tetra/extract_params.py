import Bio.PDB as pdb
import gzip
import pandas as pd
import os
import json

from angles import angles_tetra
from utils import printProgressBar

parser = pdb.PDBParser(QUIET=True)
				
def process_filtered():
        output_f = open('stat.txt', 'w')
	# dir with pdb files
	start_dir = '../pdb/'

	filtered = 'lables.txt'
	tab = pd.read_table(filtered, header=None, names=['pdbID', 'chainID'])
	tab.sort_values('pdbID', inplace=True)
	
	currpdbID = None
	for i in range(tab['pdbID'].values.size):
		pdbID = tab['pdbID'].values[i]
		chainID = tab['chainID'].values[i]
		filename = pdbID+'.pdb.gz'

                print(pdbID,chainID)

                #structure = list()
		if currpdbID is None or currpdbID != pdbID:
			currpdbID = pdbID
			
                        print(filename)
			if not os.path.isfile(start_dir + filename):
				print('file %s not found' % start_dir + filename)
				continue
			
			# get structure
			args = filename.split('.')
			if len(args) != 0 and args[-1] == 'gz':
				f = gzip.open(start_dir + filename, 'r')
			else:
				f = open(start_dir + filename, 'r')
			structure = parser.get_structure('1', f)
				
		# process chain
                chain = None
		for model in structure:
			try:
				chain = model[chainID]
				break
			except KeyError:
				continue
                if chain is None: continue
		tetraArr, anglesArr = angles_tetra(chain)
		
		if tetraArr is None: continue
		
		for i in range(len(tetraArr)):
			arr = [tetraArr[i]] + anglesArr[i]
			output_f.write(' '.join(map(str,arr)) + '\n')
		
	
if __name__ == '__main__':
	process_filtered()
