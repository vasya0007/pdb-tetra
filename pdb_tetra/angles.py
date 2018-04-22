from math import sqrt, acos, pi

import numpy as np

from extract_seq import extract_seq

EPS = 1.0e-6

def angles_tetra(chain):
	""" return torsional angles phi, psi for each tetrapeptide in chain
			in this version we use six internal torsional angles:
			psi1, phi2, psi2, phi3, psi3, phi4
			(can be easily extended to eight angles)
	"""
	
	seq = extract_seq(chain, full=True)
	
	resIDList = [residue.get_id() for residue in chain]
	
	tetraArr = []
	anglesArr = []
	
	for i in range(len(seq) - 4):
		tetra = seq[i:i+4]
		if '-' in tetra: continue
		
		angles = []
		missingAtoms = False
		for k in range(4):
			resIter = i+k
			
			resID = resIDList[resIter]
			# print(resIter, resID)
			
			# phi
			if k > 0:
				resIDPrev = resIDList[resIter-1]
				
				try:
					v1 = chain[resID]['C'].get_coord()
					v2 = chain[resID]['CA'].get_coord()
					v3 = chain[resID]['N'].get_coord()
					v4 = chain[resIDPrev]['C'].get_coord()
				except KeyError:
					# it is possible (see 4BB9, resID=327)
					# that there are missing atoms
					missingAtoms = True
					break
				
				angle = calculate_dihedral(v1, v2, v3, v4)
				angles.append(angle)
				
			# psi
			if k < 3:
				resIDNext = resIDList[resIter+1]
				
				try:
					v1 = chain[resIDNext]['N'].get_coord()
					v2 = chain[resID]['C'].get_coord()
					v3 = chain[resID]['CA'].get_coord()
					v4 = chain[resID]['N'].get_coord()
				except KeyError:
					missingAtoms = True
					break
				
				angle = calculate_dihedral(v1, v2, v3, v4)
				angles.append(angle)
		
		if missingAtoms: continue
		
		tetraArr.append(tetra)
		anglesArr.append(angles)
		
	return [tetraArr, anglesArr]

def calculate_angle(r1, r2, r3):
	v1 = r1-r2
	v2 = r3-r2
	v1 /= norm(v1)
	v2 /= norm(v2)
	theta = acos(np.dot(v1, v2)) / np.pi * 180.0
	return theta	
	
def calculate_dihedral(r1, r2, r3, r4):
	# print(r1)
	# print(r2)
	# print(r3)
	# print(r4)
	# print('angles')
	# print(calculate_angle(r1,r2,r3))
	# print(calculate_angle(r2,r3,r4))
	v1 = r1-r2
	v2 = r3-r2
	n1 = np.cross(v1,v2)
	v3 = r4-r3
	n2 = np.cross(v2,v3)
	n1 /= norm(n1)
	n2 /= norm(n2)
	# print(n1)
	# print(n2)
	# print(np.dot(n1,n2))
	
	try:
		theta = acos(np.dot(n1, n2)) / np.pi * 180.0
	except ValueError:
		# case: acos(+/-1.0) undefined
		if abs(np.dot(n1,n2)+1.0) < EPS:
			theta = 180.0
			return theta
		elif abs(np.dot(n1,n2)-1.0) < EPS:
			theta = 0.0
			return theta
		else:
			print(n1)
			print(n2)
			print(np.dot(n1,n2))
			raise ValueError
	
	n3 = np.cross(v2, n1)
	if np.dot(n3, n2) < 0: theta = 360.0 - theta
	return theta
	
def norm(vec):
	return sqrt(sum(vec * vec))
