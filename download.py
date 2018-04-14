import urllib
from .util import printProgressBar
import os.path
import pandas
import re

def downloadID(pdbID, output_dir='./pdb/', check=True):
	url = 'https://files.rcsb.org/download/%s.pdb.gz' % pdbID
	output = output_dir + '%s.pdb.gz' % pdbID
	if check and os.path.isfile(output):
		return
	urllib.urlretrieve(url, output)
	# error = re.match('404 Not Found')


	
	
if __name__ == '__main__':
	# all DNA-protein (and not only protein ...) complexes from NDB
	# input = 'ndb-dna-protein-complexes-all.txt'

	# all proteins (and short peptides) with resolution <= 1.5 A
	input = 'pdb-protein-resolution-less-1.5A.txt'

	tab = pandas.read_table(input, usecols=[0], header=None)
	nrows = tab.values[:,0].size


	for i in range(nrows):
		pdbID = tab.values[i,0]
		
		downloadID(pdbID)
		
		printProgressBar(i, nrows, length=50, fill=chr(219))


