
d1 = {
	"Ala" : "A",
	"Arg" : "R",
	"Asn" : "N",
	"Asp" : "D",
	"Asx" : "B",
	"Cys" : "C",
	"Glu" : "E",
	"Gln" : "Q",
	"Glx" : "Z",
	"Gly" : "G",
	"His" : "H",
	"Ile" : "I",
	"Leu" : "L",
	"Lys" : "K",
	"Met" : "M",
	"Phe" : "F",
	"Pro" : "P",
	"Ser" : "S",
	"Thr" : "T",
	"Trp" : "W",
	"Tyr" : "Y",
	"Val" : "V"
}
d1 =  {k.upper() : v for k,v in d1.items()}
d1r = {v.upper() : k for k,v in d1.items()}


aminoAcidListLong  = 'ARNDBCEQZGHILKMFPSTWYV'
# B=D, Q=Z
aminoAcidListShort = 'ARNDCEQGHILKMFPSTWYVBZ'
aminoAcid = {aminoAcidListShort[iter] : iter for iter in range(len(aminoAcidListShort)-2)}
aminoAcid['B'] = aminoAcid['D']
aminoAcid['Z'] = aminoAcid['Q']