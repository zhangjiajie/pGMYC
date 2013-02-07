#! /usr/bin/env python

def extract_ref_alignment(nfin, nfout):
	fin=open(nfin,"r")
	fout=open(nfout,"w")
	l=fin.readline()
	seqlen = l.split()[1]
	seq_list=[]
	l=fin.readline().strip()
	while l!="":
		#print(l)
		ls=l.split()
		tname=ls[0]
		if tname.endswith("r"):
			seq_list.append(l)
		l=fin.readline().strip()
	fout.write(repr(len(seq_list)) + "  " + seqlen +"\n")

	for seq in seq_list:
		fout.write(seq+"\n")

	fout.close()
	fin.close()
