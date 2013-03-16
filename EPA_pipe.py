#! /usr/bin/env python
import json
import os 
import glob
import EXP
import GMYC
import sys
import subprocess
import re
#from EPA_ME_Test import epa_me_species_counting
from ete2 import Tree
from ete2 import SeqGroup
from subprocess import call
from collections import deque

#unused
def token(tree):
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


#unused
def grab_sub_tree(tree, leafs):
	tks = token(tree)
	#stak = deque()
	#stak.append(tks[0])
	start_taxa = -1
	end_taxa = -1
	for i in range(len(tks))[1:]:
		if (len(leafs)!=0) and (tks[i] in leafs) and (start_taxa == -1):
			leafs.remove(tks[i])
			start_taxa = i - 1 
		elif (len(leafs)==1) and (tks[i] in leafs) and (end_taxa == -1):
			end_taxa = i + 3
		else:
			if tks[i] in leafs:
				leafs.remove(tks[i])
		
	subs = deque(tks[start_taxa:(end_taxa + 1)])
	print (subs)
	left_cnt = 0
	right_cnt = 0
	coma_cnt = 0
	sumb = 0 
	for tk in subs:
		if "(" == tk:
			left_cnt = left_cnt + 1
		elif ")" == tk:
			right_cnt = right_cnt + 1
		elif "," == tk:
			coma_cnt = coma_cnt + 1
	if left_cnt > right_cnt:
		sumb = left_cnt * 2
		append_num = left_cnt - right_cnt
		append = [")"] * append_num
		subs.extend(append)
	elif right_cnt > left_cnt:
		sumb = right_cnt * 2
		append_num = right_cnt - left_cnt
		append = ["("] * append_num
		subs.extendleft(append)
	else:
		sumb = left_cnt + right_cnt
	
	if sumb!= (coma_cnt*2):
		subs.append(")")
		subs.appendleft("(")
	
	print("".join(subs))


def gen_alignment(seq_names = [], alignment = SeqGroup(), outputfile = "epa_parser_alignments.out"):
	"""generate alignment from the input taxa name list - seq_name, and SeqGroup - alignment"""
	newalign = SeqGroup(format='phylip')
	for taxa in seq_names:
		seq = alignment.get_seq(taxa)
		newalign.set_seq(taxa, seq)
	newalign.write(format='iphylip_relaxed', outfile = outputfile)
	return outputfile


def rm_redudent_seqs(aln_file_name):
	call(["bin/raxmlHPC-SSE3","-m","GTRGAMMA","-s",aln_file_name,"-f","c","-n","ck"], stdout=open(os.devnull, "w"), stderr=subprocess.STDOUT)
	if os.path.exists(aln_file_name+".reduced"):
		os.remove("RAxML_info." + "ck")
		return aln_file_name+".reduced"
	else:
		return aln_file_name
		
def rm_redudent_seqs_m(aln_file_name):
	call(["bin/raxmlHPC-PTHREADS-SSE3","-m","GTRGAMMA","-s",aln_file_name,"-f","c","-n","ck", "-T", "2"], stdout=open(os.devnull, "w"), stderr=subprocess.STDOUT)
	if os.path.exists(aln_file_name+".reduced"):
		os.remove("RAxML_info." + "ck")
		return aln_file_name+".reduced"
	else:
		return aln_file_name


#step0: pre-process simulated alignment:
def pre_pro_aln(nfin, nfout, numcpu = 1):
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
	faln = None 
	if numcpu == 1:
		faln = rm_redudent_seqs(nfin)
	else:
		faln = rm_redudent_seqs_m(nfin)
	return faln


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


#step1: extract ref alignment
def extract_ref_query_alignment(nfin, nfout):
	fin = open(nfin,"r")
	foutr = open(nfout+".ref.afa","w")
	foutq = open(nfout+".query.afa", "w")
	 
	l=fin.readline()
	seqlen = l.split()[1]
	
	rseq_list=[]
	qseq_list=[]
	
	l=fin.readline().strip()
	while l!="":
		ls=l.split()
		tname=ls[0]
		if tname.endswith("r"):
			rseq_list.append(l)
		else:
			qseq_list.append(l)
			
		l=fin.readline().strip()
	
	for seq in rseq_list:
		seqs = seq.split()
		name = seqs[0]
		s = seqs[1]
		foutr.write(">" + name + "\n")
		foutr.write(s + "\n")
	
	for seq in qseq_list:
		seqs = seq.split()
		name = seqs[0]
		s = seqs[1]
		foutq.write(">" + name + "\n")
		foutq.write(s + "\n")
		
	foutr.close()
	foutq.close()
	fin.close()
	return nfout + nfout+".ref.afa", nfout+".query.afa"


#step2: build ref tree
def build_ref_tree(nfin, nfout):
	call(["raxmlHPC-SSE3","-m","GTRGAMMA","-s",nfin,"-n",nfout,"-p", "1234"], stdout=open(os.devnull, "w"), stderr=subprocess.STDOUT)
	os.rename("RAxML_bestTree."+nfout, nfout + ".tre")
	os.remove("RAxML_info." + nfout)
	os.remove("RAxML_log." + nfout)
	os.remove("RAxML_parsimonyTree." + nfout)
	os.remove("RAxML_result." + nfout)
	return nfout + ".tre"


def build_ref_tree_m(nfin, nfout, numcpu = "2"):
	call(["raxmlHPC-PTHREADS-SSE3","-m","GTRGAMMA","-s",nfin,"-n",nfout,"-p", "1234", "-T", numcpu], stdout=open(os.devnull, "w"), stderr=subprocess.STDOUT)
	os.rename("RAxML_bestTree."+nfout, nfout + ".tre")
	os.remove("RAxML_info." + nfout)
	os.remove("RAxML_log." + nfout)
	os.remove("RAxML_parsimonyTree." + nfout)
	os.remove("RAxML_result." + nfout)
	return nfout + ".tre"

#build tree with -g
def build_constrain_tree(nsfin, ntfin, nfout):
	call(["raxmlHPC-SSE3","-m","GTRGAMMA","-s",nsfin, "-g", ntfin, "-n",nfout,"-p", "1234"], stdout=open(os.devnull, "w"), stderr=subprocess.STDOUT)
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


#step4: extrac EPA placement & build the phylogenetic tree -g option
def extract_placement2(nfin_place, nfin_aln, nfout):
	jsondata = open (nfin_place)
	align_orgin = SeqGroup(sequences = nfin_aln, format="phylip_relaxed")
	data = json.load(jsondata)
	placements = data["placements"]
	tree = data["tree"]
	
	ete_tree = tree.replace("{", "[&&NHX:B=")
	ete_tree = ete_tree.replace("}", "]")
	root = Tree(ete_tree, format=1)
	leaves = root.get_leaves()
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
			
			#genrate the newwick string to be inserted into the ref tree
			multi_fcating = "("
			for seqname in seqset:
				multi_fcating = multi_fcating + seqname + ","
			multi_fcating = multi_fcating[:-1] 
			multi_fcating = "{" + repr(seqset_name) + "}," + multi_fcating + ")"
			mtfc_tree = tree.replace("{" + repr(seqset_name) + "}", multi_fcating)
			mtfc_tree = rep.sub("", mtfc_tree)
			
			#The target taxa names
			curr_placement = []
			for n in seqset:
				curr_placement.append(n)
			
			#Now seqset contains refseq
			for leaf in leaves:
				seqset.append(leaf.name)
			
			#generate aligment with ref seqs
			alnname = gen_alignment(seq_names = seqset, alignment = align_orgin, outputfile = nfout + repr(i) + ".phy")
			
			#write multifurcating tree
			mtfc_out = open(nfout + repr(i) + ".mttree", "w")
			mtfc_out.write(mtfc_tree)
			mtfc_out.close()
			
			#raxml constrait search
			trename = build_constrain_tree(alnname, nfout + repr(i) + ".mttree", nfout+repr(i))
			
			#read in the fully resolved tree
			full_tree = Tree(trename, format=1)
			
			#the place where the tree can be safely rooted
			ref_node = full_tree.get_leaves_by_name(leaves[0].name)[0]
			
			#reroot 
			full_tree.set_outgroup(ref_node)
			
			#find the common ancestor of the target taxa
			leafA = full_tree.get_leaves_by_name(curr_placement[0])[0]
			leaflist = []
			for n in curr_placement[1:]:
				leaflist.append(full_tree.get_leaves_by_name(n)[0])
			common = leafA.get_common_ancestor(leaflist)
			common.up = None
			common.write(outfile= nfout + repr(i) + ".subtree", format=5)
			
			os.remove(nfout + repr(i) + ".mttree")
			os.remove(trename)
			
			aln_fnames.append(alnname)
			tre_fnames.append(nfout + repr(i) + ".subtree")
	
	return num_epa_edges, aln_fnames, tre_fnames


def estimate_ref_exp_rate(nfin):
	ref_model = EXP.exponential_mixture(tree = nfin)
	spe_rate = ref_model.null_model()
	return spe_rate


#step5: test mix_exp model
def mix_exp_model(nfin_tree_list, nfin_ref_tree = None, re_root = True):
	spe_rate = -1
	spe_num_list = []
	true_spe_num_list = []
	if nfin_ref_tree!= None:
		spe_rate = estimate_ref_exp_rate(nfin_ref_tree)
	
	for ntre in nfin_tree_list:
		epa_exp = None
		if spe_rate > 0:
			epa_exp = EXP.exponential_mixture(ntre, sp_rate = spe_rate, fix_sp_rate = True)
		else:
			epa_exp = EXP.exponential_mixture(ntre)
			
		epa_exp.search(reroot = re_root, strategy = "Brutal")
		num_spe = epa_exp.count_species(print_log = False)
		spe_num_list.append(num_spe)
		#species_list = epa_exp.species_list
		
		leaves = epa_exp.tree.get_leaves()
		true_spe = set([leaves[0].name.split(".")[0]])
		for leaf in leaves:
			true_spe.add(leaf.name.split(".")[0])
		true_spe_num = len(true_spe)
		true_spe_num_list.append(true_spe_num)
		
	
	return spe_num_list, true_spe_num_list


def logging(tsl, sl, sfout):
	fout = open(sfout, "a")
	for i in range(len(tsl)):
		tn = tsl[i]
		pn = sl[i]
		fout.write(repr(tn) + "	" + repr(pn) + "\n")
	fout.close()


def dppdiv():
	#"./dppdiv-pthreads-sse3 -T 4 -in testdata/sample.phy -tre testdata/sample.tree -cal testdata/sample.cal -n 1000 -out output/sample"
	pass


def r8s(sfin_tree, sfout):
	#r8s/r8s -b -f datafile
	t = Tree(sfin_tree, format = 1 )
	tree_str = t.write(format=5)
	leaf_names = t.get_leaf_names()
	r8s_data_file = sfin_tree + ".r8s.temp"
	fout = open(r8s_data_file, "w")
	fout.write("#nexus\n")
	fout.write("begin trees;\n")
	fout.write("tree single_Haplo = " + tree_str + "\n")
	fout.write("end;\n")
	fout.write("\n")
	fout.write("begin r8s;\n")
	fout.write("blformat lengths=persite nsites=1000 ultrametric=no;\n")
	fout.write("mrca root " + leaf_names[0] + " " + leaf_names[1] + ";\n")
	fout.write("fixage taxon=root age=100;\n")
	fout.write("divtime method=pl algorithm=tn;\n")
	fout.write("set smoothing=0.5 penalty=log checkGradient=yes;\n")
	fout.write("describe plot=chrono_description;\n")
	fout.write("end;\n")
	fout.close()
	call(["r8s/r8s","-b","-f",r8s_data_file], stdout=open("r8slog", "w"), stderr=subprocess.STDOUT)
	fin = open("r8slog")
	lines = fin.readlines()
	trees = lines[-1]
	trees = trees.split("=")[1].strip()
	fin.close()
	fout = open(sfout, "w")
	fout.write(trees + "\n")
	fout.close()
	os.remove(r8s_data_file)
	os.remove("r8slog")
	return sfout


def gmyc_model(nfin_tree_list, nfin_ref_tree = None, re_root = False):
	spe_num_list = []
	true_spe_num_list = []
	for ntre in nfin_tree_list:
		umtree = r8s(sfin_tree = ntre, sfout = "gmyc_temp" + ntre)
		num_spe = GMYC.gmyc(tree = umtree)
		spe_num_list.append(num_spe)
		
		t =Tree(ntre, format = 1)
		leaves = t.get_leaves()
		true_spe = set([leaves[0].name.split(".")[0]])
		for leaf in leaves:
			true_spe.add(leaf.name.split(".")[0])
		true_spe_num = len(true_spe)
		true_spe_num_list.append(true_spe_num)
		#os.remove(umtree)
		
	return spe_num_list, true_spe_num_list


#step6: batch test normal me
def batch_test_me(folder="./", suf = "phy", num_spe_tree = 10, sout = "log.txt"):
	phyl = glob.glob(folder + "*." + suf)
	num_tre = len(phyl)
	#num_err_exp = 0
	num_err_epa = 0
	num_correct = 0
	num_placement = 0
	#true_num_sp = 0
	cnt = 1 
	for phy in phyl:
		print("Testing " + repr(cnt) + " / " + repr(len(phyl))  )
		cnt = cnt + 1
		fin1 = pre_pro_aln(nfin=phy, nfout="temp1")
		fin2 = extract_ref_alignment(nfin = fin1 , nfout = "temp2.phy", num_prune = 2)
		fin3 = build_ref_tree(nfin = fin2, nfout = "temp3")
		fin4 = run_EPA(nfin_aln =fin1 , nfin_tree=fin3 , nfout="p3")
		num_e, alns, tres = extract_placement(nfin_place = fin4, nfin_aln = fin1, nfout = "p4")
		
		sp_num_l, true_sp_num_l = mix_exp_model(tres, nfin_ref_tree = fin3)
		num_placement = num_placement + len(true_sp_num_l)
		logging(true_sp_num_l, sp_num_l, sout)
		
		num_diff_epa = num_e - num_spe_tree
		if num_diff_epa > 0:
			num_err_epa = num_err_epa + num_diff_epa
		
		for i in range(len(sp_num_l)):
			print("True: " + repr(true_sp_num_l[i]) + "  -- Predic: " + repr(sp_num_l[i]) )
			if sp_num_l[i] == true_sp_num_l[i]:
				num_correct = num_correct + 1
		
		for aln in alns:
			os.remove(aln)
		for tre in tres:
			os.remove(tre)
		jk1 = glob.glob("*.reduced")
		for jkfile in jk1:
			os.remove(jkfile)
		os.remove("temp2.phy")
		os.remove("temp3.tre")
		os.remove("p3.jplace")
	
	err_rate_epa = float(num_err_epa) / float(num_spe_tree * num_tre)
	correct_rate = float(num_correct) / float(num_placement)
	
	print("Errors made by EPA: " + repr(num_err_epa) + "	" + repr(err_rate_epa))
	print("Correct EXP estimate: " + repr(num_correct) + "	" + repr(correct_rate))


#step6: batch test me using -g
def batch_test_me_g(folder="./", suf = "phy", num_spe_tree = 10, sout = "log.txt"):
	phyl = glob.glob(folder + "*." + suf)
	num_tre = len(phyl)
	#num_err_exp = 0
	num_err_epa = 0
	num_correct = 0
	num_placement = 0
	#true_num_sp = 0
	cnt = 1 
	for phy in phyl:
		print("Testing " + repr(cnt) + " / " + repr(len(phyl))  )
		cnt = cnt + 1
		fin1 = pre_pro_aln(nfin=phy, nfout="temp1")
		fin2 = extract_ref_alignment(nfin = fin1 , nfout = "temp2.phy", num_prune = 2)
		fin3 = build_ref_tree(nfin = fin2, nfout = "temp3")
		fin4 = run_EPA(nfin_aln =fin1 , nfin_tree=fin3 , nfout="p3")
		num_e, alns, tres = extract_placement2(nfin_place = fin4, nfin_aln = fin1, nfout = "p4")
		
		sp_num_l, true_sp_num_l = mix_exp_model(tres, nfin_ref_tree = fin3)
		num_placement = num_placement + len(true_sp_num_l)
		logging(true_sp_num_l, sp_num_l, sout)
		
		num_diff_epa = num_e - num_spe_tree
		if num_diff_epa > 0:
			num_err_epa = num_err_epa + num_diff_epa
		
		for i in range(len(sp_num_l)):
			print("True: " + repr(true_sp_num_l[i]) + "  -- Predic: " + repr(sp_num_l[i]) )
			if sp_num_l[i] == true_sp_num_l[i]:
				num_correct = num_correct + 1
		
		for aln in alns:
			os.remove(aln)
		for tre in tres:
			os.remove(tre)
		jk1 = glob.glob("*.reduced")
		for jkfile in jk1:
			os.remove(jkfile)
		os.remove("temp2.phy")
		os.remove("temp3.tre")
		os.remove("p3.jplace")
	
	err_rate_epa = float(num_err_epa) / float(num_spe_tree * num_tre)
	correct_rate = float(num_correct) / float(num_placement)
	
	print("Errors made by EPA: " + repr(num_err_epa) + "	" + repr(err_rate_epa))
	print("Correct EXP estimate: " + repr(num_correct) + "	" + repr(correct_rate))
	

def batch_test_gmyc(folder="./", suf = "phy", num_spe_tree = 10, sout = "log.txt"):
	phyl = glob.glob(folder + "*." + suf)
	num_correct = 0
	for phy in phyl:
		fin1 = pre_pro_aln(nfin=phy, nfout="temp1")
		fin2 = extract_ref_alignment(nfin = fin1 , nfout = "temp2.phy", num_prune = 2)
		fin3 = build_ref_tree(nfin = fin2, nfout = "temp3")
		fin4 = run_EPA(nfin_aln =fin1 , nfin_tree=fin3 , nfout="p3")
		num_e, alns, tres = extract_placement2(nfin_place = fin4, nfin_aln = fin1, nfout = "p4")
		sp_num_l, true_sp_num_l = gmyc_model(tres)
		for aln in alns:
			os.remove(aln)
		for tre in tres:
			os.remove(tre)
		jk1 = glob.glob("*.reduced")
		for jkfile in jk1:
			os.remove(jkfile)
		for i in range(len(sp_num_l)):
			if sp_num_l[i] == true_sp_num_l[i]:
				num_correct = num_correct + 1
				print("True: " + repr(true_sp_num_l[i]) + "  -- Predic: " + repr(sp_num_l[i]) )
			else:
				print("True: " + repr(true_sp_num_l[i]) + "  -- Predic: " + repr(sp_num_l[i]) )
		
		os.remove("temp2.phy")
		os.remove("temp3.tre")
		os.remove("p3.jplace")
		print("Correct EXP estimate: " + repr(num_correct))


#This will do batch test of mix exp model on original trees with 10 species
def batch_mix_exp(folder="./", suf = "phy", num_spe_tree = 10, sout = "log.txt", t = "1"):
	phyl = glob.glob(folder + "*." + suf)
	print(folder + "*." + suf) 
	rt_correct = 0
	rt_num = 0
	for phy in phyl:
		
		fin1 = pre_pro_aln(nfin=phy, nfout="temp1", numcpu = t)
		gt = ground_truth(fin1)
		
		fin2 = None
		if t == "1":
			fin2 = build_ref_tree(nfin = fin1, nfout = "temp2")
		else:
			fin2 = build_ref_tree_m(nfin = fin1, nfout = "temp2", numcpu = t)
		me = EXP.exponential_mixture(fin2)
		me.search(reroot = True, strategy = "H0")
		num_spe = me.count_species(print_log = False)
		splist = me.species_list
		
		num_correct = 0 
		for spe in splist:
			corr = gt.is_correct(spe)
			if corr:
				num_correct = num_correct + 1
		os.remove(fin2)
		print("correct: " + repr(num_correct))
		print("delimit: " + repr(len(splist)))
		print("overesti: " + repr((float(len(splist)) - float(gt.get_num_species()))/float(gt.get_num_species())))
		rt_num = rt_num + (float(len(splist)) - float(gt.get_num_species()))/float(gt.get_num_species())
		rt_correct = rt_correct + float(num_correct)/float(gt.get_num_species())
		
	print("Average correct ration: "  +  repr(rt_correct/float(len(phyl))))
	print("Average overestimate: "  +  repr(rt_num/float(len(phyl))))
	print("num_correct:" + repr(rt_correct))
	print("num_samples:" + repr(len(phyl)))


class ground_truth:
	def __init__(self, refaln, type = ""):
		if type == "fasta":
			self.aln = SeqGroup(sequences=refaln)
		else:
			self.aln = SeqGroup(sequences=refaln, format='phylip_relaxed')
		self.true_spe = {}
		self._get_truth()
		
	
	def _get_truth(self):
		for entr in self.aln.get_entries():
			name = entr[0]
			gid = name.split(".")[0]
			self.true_spe[gid] = []
		
		for entr in self.aln.get_entries():
			name = entr[0]
			gid = name.split(".")[0]
			group = self.true_spe[gid]
			group.append(name)
			self.true_spe[gid] = group
	
	def is_correct(self,names):
		#*R*
		newnames = []
		for name in names:
			if name.startswith("*R*"):
				pass
			else:
				newnames.append(name)
			
		names_set = set(newnames)
		for key in self.true_spe.keys():
			sps = self.true_spe[key]
			sps_set = set(sps)
			if names_set == sps_set:
				return True
		return False
	
	def get_num_species(self):
		return len(self.true_spe.keys())


def batch_gmyc_umtree(folder="./", suf = "simulate_tree", num_spe_tree = 10, sout = "log.txt"):
	phyl = glob.glob(folder + suf +"*")
	rt_correct = 0
	rt_num = 0
	for phy in phyl:
		treename = phy.split("/")[-1]
		alnname = treename.split("_")[-1]
		palnname = alnname.split(".")
		alnname = folder + "simulated_set_" + palnname[0] + "." + palnname[1] + "_" + palnname[2] + ".phy"
		print(alnname)
		c_alnname = alnname + ".c"
		fin=open(alnname,"r")
		fout=open(c_alnname,"w")
		l=fin.readline().strip()
		while l!="":
			fout.write(l + "\n")
			l=fin.readline().strip()
		fin.close()
		fout.close()
		os.remove(alnname)
		os.rename(c_alnname, alnname)
		
		
		gt = ground_truth(alnname)
		
		spes = GMYC.gmyc(tree = phy)
		
		num_correct = 0 
		for spe in spes:
			corr = gt.is_correct(spe)
			if corr:
				num_correct = num_correct + 1
		
		print("correct: " + repr(num_correct))
		print("delimit: " + repr(len(spes)))
		print("overesti: " + repr((float(len(spes)) - float(gt.get_num_species()))/float(gt.get_num_species())))
		rt_num = rt_num + (float(len(spes)) - float(gt.get_num_species()))/float(gt.get_num_species())
		rt_correct = rt_correct + float(num_correct)/float(gt.get_num_species())
		
	print("Average correct ration: "  +  repr(rt_correct/float(len(phyl))))
	print("Average overestimate: "  +  repr(rt_num/float(len(phyl))))
	print("num_correct: " + repr(rt_correct))
	print("num_samples: " + repr(len(phyl)))


def batch_gmyc_umtree2(folder="./", suf = "simulate_tree", num_spe_tree = 10, sout = "log.txt"):
	phyl = glob.glob(folder + suf +"*")
	rt_correct = 0
	rt_num = 0
	for phy in phyl:
		treename = phy.split("/")[-1]
		alnname = treename.split("_")[-1]
		palnname = alnname.split(".")
		alnname = folder + "simulated_set_" + palnname[0] + "_" + palnname[1]  + ".phy"
		print(alnname)
		c_alnname = alnname + ".c"
		fin=open(alnname,"r")
		fout=open(c_alnname,"w")
		l=fin.readline().strip()
		while l!="":
			fout.write(l + "\n")
			l=fin.readline().strip()
		fin.close()
		fout.close()
		os.remove(alnname)
		os.rename(c_alnname, alnname)
		
		
		gt = ground_truth(alnname)
		
		spes = GMYC.gmyc(tree = phy)
		
		num_correct = 0 
		for spe in spes:
			corr = gt.is_correct(spe)
			if corr:
				num_correct = num_correct + 1
		
		print("correct: " + repr(num_correct))
		print("delimit: " + repr(len(spes)))
		print("overesti: " + repr((float(len(spes)) - float(gt.get_num_species()))/float(gt.get_num_species())))
		rt_num = rt_num + (float(len(spes)) - float(gt.get_num_species()))/float(gt.get_num_species())
		rt_correct = rt_correct + float(num_correct)/float(gt.get_num_species())
		
	print("Average correct ration: "  +  repr(rt_correct/float(len(phyl))))
	print("Average overestimate: "  +  repr(rt_num/float(len(phyl))))
	print("num_correct: " + repr(rt_correct))
	print("num_samples: " + repr(len(phyl)))



if __name__ == "__main__":
	
	extract_ref_query_alignment(nfin = "/home/zhangje/Desktop/simulated_set_10_2.phy", nfout = "/home/zhangje/Desktop/cao")
	#epa_me_species_counting(refaln = "/home/zhangje/Desktop/test/ref.afa", queryaln = "/home/zhangje/Desktop/test/query.afa", folder="/home/zhangje/Desktop/test/" )
	
	
	if len(sys.argv) < 6: 
		print("usage: ./EXP_pipe.py -folder <folder contain test files> -suf <suffix of the test files> -num_spe <num of species per test file> -m <epa_me/epa_me_g/epa_gmyc/me/gmyc> -o <output file> -T <num cpus> -pvflag")
		sys.exit()
	
	sfolder = "./"
	ssuf = "phy"
	snum_spe = 10
	method = "me"
	sfout = "log.txt"
	sT = "1"
	pvflag = False
	
	for i in range(len(sys.argv)):
		if sys.argv[i] == "-folder":
			i = i + 1
			sfolder = sys.argv[i]
		elif sys.argv[i] == "-suf":
			i = i + 1
			ssuf = sys.argv[i]
		elif sys.argv[i] == "-num_spe":
			i = i + 1
			snum_spe = int(sys.argv[i])
		elif sys.argv[i] == "-m":
			i = i + 1
			method = sys.argv[i]
		elif sys.argv[i] == "-o":
			i = i + 1
			sfout = sys.argv[i]
		elif sys.argv[i] == "-T":
			i = i + 1
			sT = sys.argv[i]
		elif sys.argv[i] == "-pvflag":
			pvflag = True
	
	if method == "me":
		batch_mix_exp(folder = sfolder, suf = ssuf, num_spe_tree = snum_spe, sout = sfout, t = sT)
	elif method == "epa_me":
		batch_test_me(folder = sfolder, suf = ssuf, num_spe_tree = snum_spe, sout = sfout)
	elif method == "epa_me_g":
		batch_test_me_g(folder = sfolder, suf = ssuf, num_spe_tree = snum_spe, sout = sfout)
	elif method == "epa_gmyc":
		batch_test_gmyc(folder = sfolder, suf = ssuf, num_spe_tree = snum_spe, sout = sfout)
	elif method == "gmyc":
		if pvflag:
			batch_gmyc_umtree(folder = sfolder, suf = ssuf, num_spe_tree = snum_spe, sout = sfout)
		else:
			batch_gmyc_umtree2(folder = sfolder, suf = ssuf, num_spe_tree = snum_spe, sout = sfout)

