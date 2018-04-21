import Bio.PDB as pdb
import gzip
import pandas as pd
import os

from .angles import angles_tetra

parser = pdb.PDBParser(QUIET=True)

f = open('135L.pdb', 'r')

structure = parser.get_structure('1', f)

for model in structure:
	for chain in model:
		# print(chain.get_id())
		# atoms = chain.get_atoms()
		# l = chain.get_list()

		res = angles_tetra(chain)
		print(res)
		
		# print(atoms.next().get_coord())
		#for r in list(res):
		#	print(res.get_resname())