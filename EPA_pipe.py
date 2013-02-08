#! /usr/bin/env python
import json
import os 
import glob
import EXP
import sys
from ete2 import Tree
from ete2 import SeqGroup
from subprocess import call
import subprocess

def gen_alignment(seq_names = [], alignment = SeqGroup(), outputfile = "epa_parser_alignments.out"):
	"""generate alignment from the input taxa name list - seq_name, and SeqGroup - alignment"""
	newalign = SeqGroup(format='phylip')
	for taxa in seq_names:
		seq = alignment.get_seq(taxa)
		newalign.set_seq(taxa, seq)
	newalign.write(format='iphylip_relaxed', outfile = outputfile)
	return outputfile

#step0: pre-process simulated alignment:
def pre_pro_aln(nfin, nfout):
	fin=open(nfin,"r")
	fout=open(nfout,"w")
	l=fin.readline().strip()
	while l!="":
		fout.write(l + "\n")
		l=fin.readline().strip()
	fin.close()
	fout.close()
	os.remove(nfin)
	os.rename(nfout, nfin)
	return nfin

#step1: extract ref alignment
def extract_ref_alignment(nfin, nfout, num_prune = 0):
	fin=open(nfin,"r")
	fout=open(nfout,"w")
	l=fin.readline()
	seqlen = l.split()[1]
	seq_list=[]
	l=fin.readline().strip()
	while l!="":
		ls=l.split()
		tname=ls[0]
		if tname.endswith("r"):
			seq_list.append(l)
		l=fin.readline().strip()
	fout.write(repr(len(seq_list)) + "  " + seqlen +"\n")

	if len(seq_list) - num_prune >3:
		for seq in seq_list[num_prune:]:
			fout.write(seq+"\n")
	else:
		print("refence number taxa < 3")

	fout.close()
	fin.close()
	return nfout


#step2: build ref tree
def build_ref_tree(nfin, nfout):
	call(["raxmlHPC-SSE3","-m","GTRGAMMA","-s",nfin,"-n",nfout,"-p", "1234"], stdout=open(os.devnull, "w"), stderr=subprocess.STDOUT)
	os.rename("RAxML_bestTree."+nfout, nfout + ".tre")
	os.remove("RAxML_info." + nfout)
	os.remove("RAxML_log." + nfout)
	os.remove("RAxML_parsimonyTree." + nfout)
	os.remove("RAxML_result." + nfout)
	return nfout + ".tre"


#step3: run EPA
def run_EPA(nfin_aln, nfin_tree, nfout):
	call(["raxmlHPC-SSE3","-m","GTRGAMMA","-s",nfin_aln,"-n",nfout,"-p", "1234", "-f", "v", "-r", nfin_tree], stdout=open(os.devnull, "w"), stderr=subprocess.STDOUT)
	os.rename("RAxML_portableTree." + nfout + ".jplace", nfout + ".jplace")
	os.remove("RAxML_classification." + nfout)
	os.remove("RAxML_classificationLikelihoodWeights." + nfout)
	os.remove("RAxML_entropy." + nfout)
	os.remove("RAxML_info." + nfout)
	os.remove("RAxML_labelledTree." + nfout)
	os.remove("RAxML_originalLabelledTree." + nfout)
	return nfout + ".jplace"
	

#step4: extrac EPA placement & build the phylogenetic tree
def extract_placement(nfin_place, nfin_aln, nfout):
	jsondata = open (nfin_place)
	align_orgin = SeqGroup(sequences = nfin_aln, format="phylip_relaxed")
	data = json.load(jsondata)
	placements = data["placements"]
	placemap = {}
	"""find how many edges are used for placement"""
	for placement in placements:
		edges = placement["p"]
		curredge = edges[0][0]
		placemap[curredge] = placemap.get(curredge, [])
	
	num_epa_edges = len(placemap)
	
	"""group taxa to edges"""
	for placement in placements:
		edges = placement["p"]
		taxa_names = placement["n"]
		curredge = edges[0][0]
		a = placemap[curredge] 
		a.extend(taxa_names)
		placemap[curredge]  = a
	
	"""output alignment"""
	groups = placemap.items()
	aln_fnames = []
	tre_fnames = []
	for i,item in enumerate(groups):
		seqset_name = item[0]
		seqset = item[1]
		if len(seqset)>3:
			alnname = gen_alignment(seq_names = seqset, alignment = align_orgin, outputfile = nfout + repr(i) + ".phy")
			trename = build_ref_tree(alnname, nfout + repr(i))
			aln_fnames.append(alnname)
			tre_fnames.append(trename)
	
	return num_epa_edges, aln_fnames, tre_fnames


def estimate_ref_exp_rate(nfin):
	ref_model = EXP.mix_exp(tree = nfin)
	spe_rate = ref_model.init_tree()
	return spe_rate
	

#step5: test mix_exp model
def mix_exp_model(nfin_tree_list, nfin_ref_tree = None):
	spe_rate = -1
	spe_num_list = []
	if nfin_ref_tree!= None:
		spe_rate = estimate_ref_exp_rate(nfin_ref_tree)
	
	for ntre in nfin_tree_list:
		epa_exp = None
		if spe_rate > 0:
			epa_exp = EXP.mix_exp(ntre, sp_rate = spe_rate, fix_sp_rate = True)
		else:
			epa_exp = EXP.mix_exp(ntre)
		epa_exp.search(reroot = True)
		num_spe = epa_exp.count_species(print_log = False)
		spe_num_list.append(num_spe)
	
	return spe_num_list
	

#step6: batch test
def batch_test(folder="./", suf = "phy", num_spe_tree = 10):
	phyl = glob.glob(folder + "*." + suf)
	num_tre = len(phyl)
	num_err_exp = 0
	num_err_epa = 0
	num_correct = 0
	for phy in phyl:
		fin1 = pre_pro_aln(nfin=phy, nfout="temp1")
		fin2 = extract_ref_alignment(nfin = fin1 , nfout = "temp2.phy")
		fin3 = build_ref_tree(nfin = fin2, nfout = "temp3")
		fin4 = run_EPA(nfin_aln =fin1 , nfin_tree=fin3 , nfout="p3")
		num_e, alns, tres = extract_placement(nfin_place = fin4, nfin_aln = fin1, nfout = "p4")
		sp_num_l = mix_exp_model(tres, nfin_ref_tree = fin3)
		for aln in alns:
			os.remove(aln)
		for tre in tres:
			os.remove(tre)
		jk1 = glob.glob("*.reduced")
		for jkfile in jk1:
			os.remove(jkfile)
		num_diff_epa = num_e - num_spe_tree
		if num_diff_epa > 0:
			num_err_epa = num_err_epa + num_diff_epa
		
		for spn in sp_num_l:
			num_err_exp = num_err_exp + spn - 1
			if spn == 1:
				num_correct = num_correct + 1
		
		os.remove("temp2.phy")
		os.remove("temp3.tre")
		os.remove("p3.jplace")
		print("Testing file: " + phy + ", accu num errors exp:" + repr(num_err_exp))
	
	err_rate_epa = float(num_err_epa) / float(num_spe_tree * num_tre)
	err_rate_exp = float(num_err_exp) / float(num_spe_tree * num_tre)
	correct_rate = float(num_correct) / float(num_spe_tree * num_tre)
	
	print("Errors made by EPA: " + repr(num_err_epa) + "	" + repr(err_rate_epa))
	print("Errors made by EXP: " + repr(num_err_exp) + "	" + repr(err_rate_exp))
	print("Correct EXP estimate: " + repr(num_correct) + "	" + repr(correct_rate))
		 

if __name__ == "__main__":
	
	batch_test(folder="./", suf = "phy", num_spe_tree = 10)
	#pre_pro_aln(nfin="pv.phy", nfout="jz.phy")
	#sys.exit()
	
	"""
	fin1 = "jz.phy"
	fin2 = extract_ref_alignment(nfin = fin1 , nfout = "p1.phy")
	fin3 = build_ref_tree(nfin = fin2, nfout = "p2")
	#print(trees)
	fin4 = run_EPA(nfin_aln =fin1 , nfin_tree=fin3 , nfout="p3")
	num_e, alns, tres = extract_placement(nfin_place = fin4, nfin_aln = fin1, nfout = "p4")
	sp_num_l = mix_exp_model(tres, nfin_ref_tree = fin3)
	#sp_num_l = mix_exp_model(tres)
	print("Number clusters id by EPA:" + repr(num_e))
	cnt = 0
	for spn in sp_num_l:
		print("EPA cluster " + repr(cnt) + ", num spe = " + repr(spn))
		cnt = cnt + 1
	
	for aln in alns:
		os.remove(aln)
	for tre in tres:
		os.remove(tre)
	jk1 = glob.glob("*.reduced")
	for jkfile in jk1:
		os.remove(jkfile)
	"""
	
