from math import sqrt, acos, pi

import numpy as np

from .extract_seq import extract_seq


def angles_tetra(chain):
	""" return torsional angles phi, psi for each tetrapeptide in chain
			in this version we use six internal torsional angles:
			psi1, phi2, psi2, phi3, psi3, phi4
			(can be easily extended to eight angles)
	"""
	
	seq = extract_seq(chain)
	
	tetraArr = []
	anglesArr = []
	
	for i in range(len(seq) - 4):
		tetra = seq[i:i+4]
		if '-' in tetra: continue
		
		angles = []
		for k in range(4):
			resID = i+k+1
			
			# phi
			if k > 0:
				v1 = chain[resID]['C'].get_coord()
				v2 = chain[resID]['CA'].get_coord()
				v3 = chain[resID]['N'].get_coord()
				v4 = chain[resID-1]['C'].get_coord()
				angle = calculate_dihedral(v1, v2, v3, v4)
				angles.append(angle)
				
			# psi
			if k < 3:
				v1 = chain[resID+1]['N'].get_coord()
				v2 = chain[resID]['C'].get_coord()
				v3 = chain[resID]['CA'].get_coord()
				v4 = chain[resID]['N'].get_coord()
				angle = calculate_dihedral(v1, v2, v3, v4)
				angles.append(angle)
		
		tetraArr.append(tetra)
		anglesArr.append(angles)
		
	return [tetraArr, anglesArr]
	
def calculate_dihedral(r1, r2, r3, r4):
	v1 = r1-r2
	v2 = r3-r2
	n1 = np.cross(v1,v2)
	v3 = r4-r3
	n2 = np.cross(v2,v3)
	n1 /= norm(n1)
	n2 /= norm(n2)
	theta = acos(np.dot(n1, n2)) / np.pi * 180.0
	n3 = np.cross(v2, n1)
	if np.dot(n3, n2) < 0: theta = 360.0 - theta
	return theta
	
def norm(vec):
	return sqrt(sum(vec * vec))
