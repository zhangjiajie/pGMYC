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
import re
from collections import deque


def token(tree):
	#print(tree)
	tree = tree[:-1]
	tks = []
	tree_len = len(tree)
	curridx = 0
	stat = 0
	buf = ""
	while curridx < tree_len:
		if stat == 0:
			if tree[curridx] == "(" or tree[curridx] == ")" or tree[curridx] == ":" or tree[curridx] == ",":
				tks.append(tree[curridx])
			else:
				stat = 1
				buf = tree[curridx]
			curridx = curridx + 1
			 
		else:
			if tree[curridx] == "(" or tree[curridx] == ")" or tree[curridx] == ":" or tree[curridx] == ",":
				stat = 0
				tks.append(buf)
				buf = ""
				tks.append(tree[curridx])
			else:
				buf = buf + tree[curridx]
			curridx = curridx + 1
	return tks

def extract_sub_tree(tree, leafs):
	tks = token(tree)
	stak = deque()
	stak.append(tks[0])
	start_taxa = -1
	end_taxa = -1
	for i in range(len(tks))[1:]:
		if (len(leafs)!=0) and (tks[i] in leafs) and (start_p == -1):
			leafs.remove(tks[i])
			start_taxa = i-1 
		
		

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
	fout.write(repr(len(seq_list[num_prune:])) + "  " + seqlen +"\n")

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

#build tree with -g
def build_constrain_tree(nsfin, ntfin, nfout):
	call(["raxmlHPC-SSE3","-m","GTRGAMMA","-s",nsfin, "-g", ntfin, "-n",nfout,"-p", "1234"], stdout=open(os.devnull, "w"), stderr=subprocess.STDOUT)
	#call(["raxmlHPC-SSE3","-m","GTRGAMMA","-s",nsfin, "-g", ntfin, "-n", nfout, "-p", "123"])
	os.rename("RAxML_bestTree."+nfout, nfout + ".tre")
	os.remove("RAxML_info." + nfout)
	os.remove("RAxML_log." + nfout)
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
	tree = data["tree"]
	placemap = {}
	"""find how many edges are used for placement"""
	for placement in placements:
		edges = placement["p"]
		curredge = edges[0][0]
		placemap[curredge] = placemap.get(curredge, [])
	
	num_epa_edges = 0 #len(placemap)
	
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
			num_epa_edges = num_epa_edges + 1
			alnname = gen_alignment(seq_names = seqset, alignment = align_orgin, outputfile = nfout + repr(i) + ".phy")
			trename = build_ref_tree(alnname, nfout + repr(i))
			aln_fnames.append(alnname)
			tre_fnames.append(trename)
	
	return num_epa_edges, aln_fnames, tre_fnames

#step4: extrac EPA placement & build the phylogenetic tree
def extract_placement2(nfin_place, nfin_aln, nfout):
	jsondata = open (nfin_place)
	align_orgin = SeqGroup(sequences = nfin_aln, format="phylip_relaxed")
	data = json.load(jsondata)
	placements = data["placements"]
	tree = data["tree"]
	#print(tree)
	
	ete_tree = tree.replace("{", "[&&NHX:B=")
	ete_tree = ete_tree.replace("}", "]")
	root = Tree(ete_tree, format=1)
	leaves = root.get_leaves()
	print("This should be the ref tree")
	print(leaves)
	placemap = {}
	"""find how many edges are used for placement"""
	for placement in placements:
		edges = placement["p"]
		curredge = edges[0][0]
		placemap[curredge] = placemap.get(curredge, [])
	
	num_epa_edges = 0 #len(placemap)
	
	"""group taxa to edges"""
	for placement in placements:
		edges = placement["p"]
		taxa_names = placement["n"]
		curredge = edges[0][0]
		a = placemap[curredge] 
		a.extend(taxa_names)
		placemap[curredge]  = a
	
	
	rep = re.compile(r"\{[0-9]*\}")
	"""output alignment"""
	groups = placemap.items()
	aln_fnames = []
	tre_fnames = []
	for i,item in enumerate(groups):
		seqset_name = item[0]
		seqset = item[1]
		if len(seqset)>3:
			num_epa_edges = num_epa_edges + 1
			
			multi_fcating = "("
			for seqname in seqset:
				#multi_fcating = multi_fcating + "," + seqname
				multi_fcating = multi_fcating + seqname + ","
			
			multi_fcating = multi_fcating[:-1] 
			
			multi_fcating = "{" + repr(seqset_name) + "}," + multi_fcating + ")"
			
			mtfc_tree = tree.replace("{" + repr(seqset_name) + "}", multi_fcating)
			
			mtfc_tree = rep.sub("", mtfc_tree)
			#print(mtfc_tree)
			
			curr_placement = []
			for n in seqset:
				curr_placement.append(n)
			
			for leaf in leaves:
				print(leaf.name)
				seqset.append(leaf.name)
			#print(seqset)
			
			alnname = gen_alignment(seq_names = seqset, alignment = align_orgin, outputfile = nfout + repr(i) + ".phy")
			mtfc_out = open(nfout + repr(i) + ".mttree", "w")
			mtfc_out.write(mtfc_tree)
			mtfc_out.close()
			
			trename = build_constrain_tree(alnname, nfout + repr(i) + ".mttree", nfout+repr(i))
			
			full_tree = Tree(trename, format=1)
			print(curr_placement)
			full_tree.prune(curr_placement)
			full_tree.write(outfile= nfout + repr(i) + ".subtree", format=5)
			"""
			leafA = full_tree.get_leaves_by_name(curr_placement[0])[0]
			leaflist = []
			for n in curr_placement[1:]:
				leaflist.append(full_tree.get_leaves_by_name(n)[0])
			common = leafA.get_common_ancestor(leaflist)
			common.up = None
			common.write(outfile= nfout + repr(i) + ".subtree", format=5)
			"""
			#aln_fnames.append(alnname)
			#tre_fnames.append(trename)
	
	return num_epa_edges, aln_fnames, tre_fnames


def estimate_ref_exp_rate(nfin):
	ref_model = EXP.mix_exp(tree = nfin)
	spe_rate = ref_model.init_tree()
	return spe_rate
	

#step5: test mix_exp model
def mix_exp_model(nfin_tree_list, nfin_ref_tree = None):
	spe_rate = -1
	spe_num_list = []
	true_spe_num_list = []
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
		leaves = epa_exp.tree.get_leaves()
		
		true_spe = set([leaves[0].name.split(".")[0]])
		for leaf in leaves:
			true_spe.add(leaf.name.split(".")[0])
		true_spe_num = len(true_spe)
		print("True num spe for " + ntre + "  :" + repr(true_spe_num))
		true_spe_num_list.append(true_spe_num)
		
	
	return spe_num_list, true_spe_num_list
	

#step6: batch test
def batch_test(folder="./", suf = "phy", num_spe_tree = 10):
	phyl = glob.glob(folder + "*." + suf)
	num_tre = len(phyl)
	num_err_exp = 0
	num_err_epa = 0
	num_correct = 0
	num_placement = 0
	true_num_sp = 0 
	for phy in phyl:
		fin1 = pre_pro_aln(nfin=phy, nfout="temp1")
		fin2 = extract_ref_alignment(nfin = fin1 , nfout = "temp2.phy", num_prune = 3)
		fin3 = build_ref_tree(nfin = fin2, nfout = "temp3")
		fin4 = run_EPA(nfin_aln =fin1 , nfin_tree=fin3 , nfout="p3")
		num_e, alns, tres = extract_placement(nfin_place = fin4, nfin_aln = fin1, nfout = "p4")
		num_placement = num_placement + num_e
		sp_num_l, true_sp_num_l = mix_exp_model(tres, nfin_ref_tree = fin3)
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
		
		for i in range(len(sp_num_l)):
			num_err_exp = num_err_exp + sp_num_l[i] - true_sp_num_l[i]
			true_num_sp = true_num_sp + true_sp_num_l[i]
			if sp_num_l[i] == true_sp_num_l[i]:
				num_correct = num_correct + 1
				print("True: " + repr(true_sp_num_l[i]) + "  -- Predic: " + repr(sp_num_l[i]) )
			else:
				print("True: " + repr(true_sp_num_l[i]) + "  -- Predic: " + repr(sp_num_l[i]) )
		
		os.remove("temp2.phy")
		os.remove("temp3.tre")
		#os.remove("p3.jplace")
		print("Testing file: " + phy + ", accu num errors exp:" + repr(num_err_exp))
	
	err_rate_epa = float(num_err_epa) / float(num_spe_tree * num_tre)
	err_rate_exp = float(num_err_exp) / float(num_spe_tree * num_tre)
	correct_rate = float(num_correct) / float(num_placement)
	
	print("Errors made by EPA: " + repr(num_err_epa) + "	" + repr(err_rate_epa))
	print("Errors made by EXP: " + repr(num_err_exp) + "	" + repr(err_rate_exp))
	print("Correct EXP estimate: " + repr(num_correct) + "	" + repr(correct_rate))
		 

if __name__ == "__main__":
	
	#batch_test(folder="./", suf = "phy", num_spe_tree = 10)
	#pre_pro_aln(nfin="pv.phy", nfout="jz.phy")
	#sys.exit()
	
	#extract_placement2(nfin_place = "p3.jplace", nfin_aln = "pipe_test.phy" , nfout = "mf")
	token("((6.10.r:0.84314588630213416209,(6.1.q,6.6.q,6.4.q,6.7.q,6.9.q,6.3.q,6.2.q,6.8.q,6.5.q),(8.10.r:0.66420620868791402369,7.10.r:0.64658585728757944633):0.36731678624631497465):0.61996417294109085194,((2.10.r:1.70135900093206116068,3.10.r:1.58117056143981926652):0.26930523710313375441,1.10.r:1.68142566753935174262):0.00000067440133086238,5.10.r:4.21404673447556810117)")
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
	
