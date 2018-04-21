import Bio.PDB as pdb
import gzip
import pandas as pd
import os

from .utils import printProgressBar
from .aminoacids import aminoAcid, d1
	
def extract_seq(chain):
	""" returns seq of aminoacids in chain if there are any """
	seq = ""
	for residue in chain:
		sym = str(residue.get_resname())
		if sym == 'HOH':
			continue
		if sym in d1.keys():
			seq += d1[sym]
		elif len(seq) != 0:
			# TODO: convert to warning
			# h.write('ID: %s chain %s: unknown residue %s\n' % (pdbID, chainName, sym))
			seq += '-'
		continue
	return seq
	
if __name__ == '__main__':
	parser = pdb.PDBParser(QUIET=True)
	
	output = 'chain_seq.txt'
	g = open(output, 'w')
	
	start_dir = './pdb/'
	files = os.listdir(start_dir)
	# print(len(files))
	
	for ifile in range(len(files)):
		filename = files[ifile]
		if filename == '.' or filename == '..': continue
		
		pdbID = filename.split('.')[0]
		
		g.write(str(pdbID) + '\n')
		
		args = filename.split('.')
		if len(args) != 0 and args[-1] == 'gz':
			f = gzip.open(start_dir + filename, 'r')
		else:
			# ? should close before
			f = open(start_dir + filename, 'r')
		
		structure = parser.get_structure('1', f)
			
		for model in structure:
			for chain in model:
				chainName = chain.get_id()
				seq = extract_seq(chain)
				g.write('%s %s\n' % (chainName, seq))
			# process only first model
			break
		
		printProgressBar(ifile, len(files), length=50, fill=chr(219))
	
	g.close()
