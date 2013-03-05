#! /usr/bin/env python
import sys
import math
import collections
import re
import json
import os 
import glob
import EXP
import subprocess
import random
from ete2 import Tree, TreeStyle, TextFace, SeqGroup, NodeStyle
from scipy.optimize import fmin
from collections import deque
from scipy import stats
from numpy import array
from subprocess import call

def gen_alignment(seq_names = [], alignment = SeqGroup(), outputfile = "epa_parser_alignments.out"):
	"""generate alignment from the input taxa name list - seq_name, and SeqGroup - alignment"""
	newalign = SeqGroup()
	for taxa in seq_names:
		seq = alignment.get_seq(taxa)
		newalign.set_seq(taxa, seq)
	newalign.write(outfile = outputfile)
	return outputfile, newalign


def gen_alignment2(seq_names = [], alignment = SeqGroup()):
	"""generate alignment from the input taxa name list - seq_name, and SeqGroup - alignment"""
	newalign = SeqGroup()
	for taxa in seq_names:
		if taxa.startswith("*R*"):
			seq = alignment.get_seq(taxa[3:])
		else:
			seq = alignment.get_seq(taxa)
		newalign.set_seq(taxa, seq)
	#newalign.write(outfile = outputfile)
	return newalign


def catalns(refs, alns, sfout):
	fout = open(sfout, "w")
	for aln in alns:
		fout.write(aln)
	fout.write("\n")
	for ref in refs:
		fout.write(ref)
	fout.close()
	return sfout


def count_and_pick_reads(align, outputfile):
	logss = ""
	numreads = 0
	l_reads = None
	max_non_gap = 0
	entries = align.get_entries()
	if len(entries) == 1:
		name = entries[0][0]
		if name.startswith("*R*"):
			return ""
	for entr in align.get_entries():
		name = entr[0]
		seq = entr[1]
		if not name.startswith("*R*"):
			numseqs = int(name.split("*")[-1])
			numreads = numreads + numseqs
		else:
			logss = "	find reference species: " + name
		
		seql = count_non_gap(seq)
		if seql > max_non_gap:
			l_reads = entr
			max_non_gap = seql
	if logss == "":
		logss = "	find new species:"
	logss = logss + " :reads number: " + repr(numreads) + "\n"
	fout = open(outputfile, "a")
	fout.write(">" + l_reads[0] + "*" + repr(numreads) + "\n")
	fout.write(l_reads[1] + "\n")
	fout.close()
	return logss


def remove_seq_len_smaller_than(f_fasta, min_l):
	fin = open(f_fasta)
	fout = open(f_fasta+".min"+repr(min_l)+".fasta", "a")
	line = fin.readline().strip()
	lastname = ""
	while line!="":
		if line.startswith(">"):
			line = line.split()[0]
			lastname = line
			#fout.write(line + "\n")
		else:
			if len(line) >= min_l:
				fout.write(lastname + "\n")
				fout.write(line + "\n")
		line = fin.readline().strip()
	fout.close()


def collapse_identical_seqs(f_fasta):
	fin = open(f_fasta)
	fout = open(f_fasta+".collapse.fasta", "a")
	seq_list=[]
	line = fin.readline().strip()
	while line!="":
		name = line
		seq = fin.readline().strip()
		flag = True
		for oneseq in seq_list:
			if oneseq[2] == seq:
				oneseq[1] = oneseq[1] + 1
				flag = False
				break
		if flag:
			seq_list.append([name, 1, seq])
		line = fin.readline().strip()
	for oneseq in seq_list:
		fout.write(oneseq[0]+"*"+repr(oneseq[1])+"\n")
		fout.write(oneseq[2]+"\n")
	
	fin.close()
	fout.close()


def processHMMseq(seqin):
	newseq = ""
	for s in seqin:
		if s == ".":
			pass
		elif s == "-":
			newseq = newseq + s
		elif s.isupper():
			newseq = newseq + s
	return newseq


def count_non_gap(seqin):
	cnt = 0
	for s in seqin:
		if s!="-":
			cnt = cnt + 1
	return cnt


def parse_HMM(f_stock):
	cnt = 0
	fin = open(f_stock)
	line = fin.readline()
	seqs = {}
	while line!="":
		if line.startswith("//"):
			break
		elif line.startswith("#"):
			pass
		elif line.startswith("\n"):
			cnt = cnt + 1
		else:
			line = line.strip()
			if cnt == 1:
				l2 = line.split()
				ss = processHMMseq(l2[1])
				seqs[l2[0]] = ss
			else:
				l2 = line.split()
				seq = seqs[l2[0]]
				ss = processHMMseq(l2[1])
				seqs[l2[0]] = seq + ss 
		line = fin.readline()
	fin.close()
	fout = open(f_stock+".fasta", "w")
	for key in seqs.keys():
		if count_non_gap(seqs[key]) >= 100:
			fout.write(">" + key + "\n")
			fout.write(seqs[key] + "\n")
	fout.close()


def chimera_removal(nuseach, nalign, nout, chimeraout):
	align = SeqGroup(nalign)
	newalign = open(nout, "w")
	chalign = open(chimeraout, "w")
	fus = open(nuseach)
	lines = fus.readlines()
	fus.close()
	for line in lines:
		its = line.split()
		c = its[-1]
		sname = its[1]
		if c == "Y" or c =="?":
			seq = align.get_seq(sname)
			chalign.write(">" + sname + "\n")
			chalign.write(seq + "\n")
		else:
			seq = align.get_seq(sname)
			newalign.write(">" + sname + "\n")
			newalign.write(seq + "\n")
	newalign.close()
	chalign.close()


#correct rooting method, this shound run after EPA
def extract_placement3(nfin_place, nfin_aln, nfout, min_lw = 0.9, logfile = "spcount.log"):
	jsondata = open (nfin_place)
	align_orgin = SeqGroup(sequences = nfin_aln)
	data = json.load(jsondata)
	placements = data["placements"]
	tree = data["tree"]
	
	ete_tree = tree.replace("{", "[&&NHX:B=")
	ete_tree = ete_tree.replace("}", "]")
	root = Tree(ete_tree, format=1)
	leaves = root.get_leaves()
	allnodes = root.get_descendants()
	allnodes.append(root)
	
	"""get refseq"""
	refseqset = []
	for leaf in leaves:
		refseqset.append(leaf.name)
	refali = gen_alignment2(seq_names = refseqset, alignment = align_orgin)
	
	placemap = {}
	"""find how many edges are used for placement"""
	for placement in placements:
		edges = placement["p"]
		curredge = edges[0][0]
		lw = edges[0][2] 
		if lw >= min_lw:
			placemap[curredge] = placemap.get(curredge, [])
	
	"""placement quality control"""
	discard_file = open(nfout+"discard.placement.txt", "w")
	"""group taxa to edges"""
	for placement in placements:
		edges = placement["p"]
		taxa_names = placement["n"]
		curredge = edges[0][0]
		lw = edges[0][2] 
		if lw >= min_lw:
			a = placemap[curredge] 
			a.extend(taxa_names)
			placemap[curredge]  = a
		else:
			discard_file.write(repr(taxa_names) + "\n")
	discard_file.close()
	
	groups = placemap.items()
	cnt_leaf = 0
	cnt_inode = 0
	
	"""check each edge""" 
	for i,item in enumerate(groups):
		seqset_name = item[0]
		seqset = item[1]
		
		"""check if placed on leaf node and find the node being placed on"""
		flag = False
		place_node = None
		for node in allnodes:
			if str(node.B) == str(seqset_name):
				place_node = node
				if node.is_leaf():
					flag = True 
				break
				
		"""find the furthest leaf of the placement node"""
		fnode = place_node.get_farthest_node()[0]
		outgroup_name = fnode.name
		
		"""find sister node"""
		snode = place_node.get_sisters()[0]
		if not snode.is_leaf():
			snode = snode.get_closest_leaf()[0]
		sister_name = snode.name
		
		"""generate aligment"""
		if flag:
			"""process leaf node placement"""
			cnt_leaf = cnt_leaf + 1
			newalign = SeqGroup()
			for taxa in seqset:
				seq = align_orgin.get_seq(taxa)
				newalign.set_seq(taxa, seq)
			if len(newalign.get_entries()) < 2:
				#count_and_pick_reads(align = newalign, outputfile = nfout + "_leaf_picked_otus.fasta")
				og_seq = align_orgin.get_seq(outgroup_name)
				sis_seq = align_orgin.get_seq(sister_name)
				newalign.set_seq("sister", sis_seq) #set the sister seqeunce to make 4 taxa
				newalign.set_seq("root_ref", og_seq) #set the outgroup name
				place_seq = align_orgin.get_seq(place_node.name)
				newalign.set_seq("*R*" + place_node.name, place_seq) #set the reference sequence name
				newalign.write(outfile = nfout + "_leaf_"+repr(cnt_leaf) + ".fasta")
			else:
				og_seq = align_orgin.get_seq(outgroup_name)
				newalign.set_seq("root_ref", og_seq) #set the outgroup name
				place_seq = align_orgin.get_seq(place_node.name)
				newalign.set_seq("*R*" + place_node.name, place_seq) #set the reference sequence name
				newalign.write(outfile = nfout + "_leaf_"+repr(cnt_leaf) + ".fasta")
		else:
			"""genrate the newwick string to be inserted into the ref tree"""
			rep = re.compile(r"\{[0-9]*\}")
			multi_fcating = "("
			for seqname in seqset:
				multi_fcating = multi_fcating + seqname + ","
			multi_fcating = multi_fcating[:-1] 
			multi_fcating = "{" + repr(seqset_name) + "}," + multi_fcating + ")"
			mtfc_tree = tree.replace("{" + repr(seqset_name) + "}", multi_fcating)
			mtfc_tree = rep.sub("", mtfc_tree)
			
			cnt_inode = cnt_inode + 1
			newalign = SeqGroup()
			for taxa in seqset:
				seq = align_orgin.get_seq(taxa)
				newalign.set_seq(taxa, seq)
			if len(newalign.get_entries()) < 2:
				count_and_pick_reads(align = newalign, outputfile = nfout + "_inode_picked_otus.fasta")
				sp_log(sfout = logfile, logs="	the palcement is on an internal node \n	find new species::reads number:1 \n")
			else:
				#og_seq = align_orgin.get_seq(outgroup_name)
				#newalign.set_seq("root_ref", og_seq)
				for entr in refali.get_entries():
					sname = entr[0]
					seqe = entr[1]
					newalign.set_seq(sname, seq)
				newalign.write(outfile = nfout + "_inode_"+repr(cnt_inode) + ".fasta")
				mtfc_out = open(nfout + "_inode_"+repr(cnt_inode) +  ".mttree", "w")
				mtfc_out.write(mtfc_tree)
				mtfc_out.close()


#extrac EPA placement & build the phylogenetic tree -g option
def extract_placement2(nfin_place, nfin_aln, nfout):
	jsondata = open (nfin_place)
	align_orgin = SeqGroup(sequences = nfin_aln)
	data = json.load(jsondata)
	placements = data["placements"]
	tree = data["tree"]
	
	ete_tree = tree.replace("{", "[&&NHX:B=")
	ete_tree = ete_tree.replace("}", "]")
	root = Tree(ete_tree, format=1)
	leaves = root.get_leaves()
	allnodes = root.get_descendants()
	allnodes.append(root)
	
	"""write refseq"""
	refseqset = []
	for leaf in leaves:
		refseqset.append(leaf.name)
	alnname, tali = gen_alignment(seq_names = refseqset, alignment = align_orgin, outputfile = nfout+".refalign.tmp.fasta")
	refalnfin = open(alnname)
	refaln = refalnfin.readlines()
	refalnfin.close()
	
	placemap = {}
	"""find how many edges are used for placement"""
	for placement in placements:
		edges = placement["p"]
		curredge = edges[0][0]
		placemap[curredge] = placemap.get(curredge, [])
	
	
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
	cnt_leaf = 0
	cnt_inode = 0 
	for i,item in enumerate(groups):
		seqset_name = item[0]
		seqset = item[1]
		
		#check if placed on leaf node
		flag = False
		for leaf in leaves:
			if str(leaf.B) == str(seqset_name):
				flag = True 
				break
		
		#genrate the newwick string to be inserted into the ref tree
		multi_fcating = "("
		for seqname in seqset:
			multi_fcating = multi_fcating + seqname + ","
		multi_fcating = multi_fcating[:-1] 
		multi_fcating = "{" + repr(seqset_name) + "}," + multi_fcating + ")"
		mtfc_tree = tree.replace("{" + repr(seqset_name) + "}", multi_fcating)
		mtfc_tree = rep.sub("", mtfc_tree)
			
		#generate aligment with ref seqs
		if flag:
			cnt_leaf = cnt_leaf + 1
			alnname, tali = gen_alignment(seq_names = seqset, alignment = align_orgin, outputfile = nfout + "_leaf_"+repr(cnt_leaf) + ".fasta")
			if len(tali.get_entries()) < 3:
				count_and_pick_reads(align = tali, outputfile = nfout + "_leaf_picked_otus.fasta")
				os.remove(alnname)
			else:
				alnfin = open(alnname)
				talns = alnfin.readlines()
				alnfin.close()
				catalns(refs = refaln, alns = talns , sfout = alnname)
				mtfc_out = open(nfout + "_leaf_"+repr(cnt_leaf) +  ".mttree", "w")
				mtfc_out.write(mtfc_tree)
				mtfc_out.close()
		else:
			cnt_inode = cnt_inode + 1
			alnname, tali = gen_alignment(seq_names = seqset, alignment = align_orgin, outputfile = nfout + "_inode_"+repr(cnt_inode) + ".fasta")
			if len(tali.get_entries()) < 3:
				count_and_pick_reads(align = tali, outputfile = nfout + "_inode_picked_otus.fasta")
				os.remove(alnname)
			else:
				alnfin = open(alnname)
				talns = alnfin.readlines()
				alnfin.close()
				catalns(refs = refaln, alns = talns , sfout = alnname)
				mtfc_out = open(nfout + "_inode_"+repr(cnt_inode) +  ".mttree", "w")
				mtfc_out.write(mtfc_tree)
				mtfc_out.close()


#build tree with -g
def build_constrain_tree(nsfin, ntfin, nfout):
	call(["/home/zhangje/bin/raxmlHPC-PTHREADS-SSE3","-m","GTRGAMMA","-s",nsfin, "-g", ntfin, "-n",nfout,"-p", "1234", "-T", "40"])#, stdout=open(os.devnull, "w"), stderr=subprocess.STDOUT)
	os.rename("RAxML_bestTree."+nfout, nfout + ".tre")
	os.remove("RAxML_info." + nfout)
	os.remove("RAxML_log." + nfout)
	os.remove("RAxML_result." + nfout)
	return nfout + ".tre"


#build tree with -g
def build_constrain_tree_l(nsfin, ntfin, nfout):
	call(["raxmlHPC-SSE3","-m","GTRGAMMA","-s",nsfin, "-g", ntfin, "-n",nfout,"-p", "1234"]) #, stdout=open(os.devnull, "w"), stderr=subprocess.STDOUT)
	os.rename("RAxML_bestTree."+nfout, nfout + ".tre")
	os.remove("RAxML_info." + nfout)
	os.remove("RAxML_log." + nfout)
	os.remove("RAxML_result." + nfout)
	return nfout + ".tre"


#build the phylogenetic tree -g option
def raxml_g(nfolder, nfout, nref_align, suf = "mttree"):
	align_orgin = SeqGroup(sequences = nref_align)
	ref_taxa = []
	for entr in align_orgin.get_entries():
		ref_taxa.append(entr[0])
	
	mttrees = glob.glob(nfolder + "*." + suf)
	cnt = 0
	for mtre in mttrees:
		cnt = cnt + 1
		align = mtre.split(".")[0] + ".fasta"
		print(align)
		#raxml constrait search
		trename = build_constrain_tree(align, mtre, "full"+repr(cnt))
		#read in the fully resolved tree
		full_tree = Tree(trename, format=1)
		all_taxa = full_tree.get_leaf_names()
		target_taxa = []
		for taxa in all_taxa:
			if taxa in ref_taxa:
				pass
			else:
				target_taxa.append(taxa)
		#the place where the tree can be safely rooted
		ref_node = full_tree.get_leaves_by_name(ref_taxa[0])[0]
		#reroot 
		full_tree.set_outgroup(ref_node)
		#find the common ancestor of the target taxa
		leafA = full_tree.get_leaves_by_name(target_taxa[0])[0]
		leaflist = []
		for n in target_taxa[1:]:
			leaflist.append(full_tree.get_leaves_by_name(n)[0])
		common = leafA.get_common_ancestor(leaflist)
		common.up = None
		common.write(outfile= mtre.split(".")[0] + ".subtree", format=5)
		os.remove(mtre)


#build a tree
def build_ref_tree_l(nfin, nfout):
	call(["raxmlHPC-SSE3","-m","GTRGAMMA","-s",nfin,"-n",nfout,"-p", "1234"], stdout=open(os.devnull, "w"), stderr=subprocess.STDOUT)
	os.rename("RAxML_bestTree."+nfout, nfout + ".tre")
	os.remove("RAxML_info." + nfout)
	os.remove("RAxML_log." + nfout)
	os.remove("RAxML_parsimonyTree." + nfout)
	os.remove("RAxML_result." + nfout)
	return nfout + ".tre"


#build a tree
def build_ref_tree(nfin, nfout):
	call(["/home/zhangje/bin/raxmlHPC-PTHREADS-SSE3","-m","GTRGAMMA","-s",nfin,"-n",nfout,"-p", "1234", "-T", "44"], stdout=open(os.devnull, "w"), stderr=subprocess.STDOUT)
	os.rename("RAxML_bestTree."+nfout, nfout + ".tre")
	os.remove("RAxML_info." + nfout)
	os.remove("RAxML_log." + nfout)
	os.remove("RAxML_parsimonyTree." + nfout)
	os.remove("RAxML_result." + nfout)
	return nfout + ".tre"


def compare_node(node):
	return node.dist


def find_lonest_br(tree):
	node_list = tree.get_descendants()
	node_list.sort(key=compare_node)
	node_list.reverse()
	rootnode = node_list[0]
	return rootnode


#build the phylogenetic tree
def raxml(nfolder, nfout, suf = "fasta"):
	#align_orgin = SeqGroup(sequences = nref_align)
	#ref_taxa = []
	#for entr in align_orgin.get_entries():
	#	ref_taxa.append(entr[0])
	
	naligns = glob.glob(nfolder + "*." + suf)
	cnt = 0
	for aln in naligns:
		cnt = cnt + 1
		trename = build_ref_tree(aln, "full"+repr(cnt))
		full_tree = Tree(trename, format=1)
		rootref = full_tree.get_leaves_by_name("root_ref")[0]
		if rootref.up.is_root():
			newrootnode = rootref.get_farthest_node()[0]
			full_tree.set_outgroup(newrootnode)
		
		rootref = full_tree.get_leaves_by_name("root_ref")[0]
		refroot_brl = rootref.dist
		full_tree.set_outgroup(rootref)
		real_tree = None
		
		for child in full_tree.get_children():
			if not child.is_leaf():
				real_tree = child
				real_tree.up = None
				real_tree.dist = 0.0
				break
		
		lnode = find_lonest_br(real_tree)
		if lnode.dist > refroot_brl:
			real_tree.set_outgroup(lnode)
			real_tree.dist = 0.0
		
		real_tree.write(outfile= aln.split(".")[0] + ".subtree", format=5)
		#os.remove(trename)


def estimate_ref_exp_rate(nfin):
	ref_model = EXP.exponential_mixture(tree = nfin)
	spe_rate = ref_model.null_model()
	return spe_rate


def sp_log(sfout, logs=""):
	f = open(sfout, "a")
	f.write(logs)
	f.write("\n")
	f.close()


def otu_picking(nfolder, nfout1, nfout2, nref_tree, n_align, suf = "subtree"):
	trees = glob.glob(nfolder + "*." + suf)
	spe_rate = estimate_ref_exp_rate(nref_tree)
	align = SeqGroup(sequences = n_align)
	logss = ""
	for tree in trees:
		logss = logss + "Searching species in tree: " + tree + ":\n"
		epa_exp = EXP.exponential_mixture(tree, sp_rate = spe_rate, fix_sp_rate = True)
		t = Tree(tree, format = 1)
		tsize = len(t.get_leaves())
		if tsize > 500:
			epa_exp.search(reroot = False, strategy = "H1")
			logss = logss + "H1\n"
		elif tsize < 20:
			epa_exp.search(reroot = False, strategy = "Brutal")
			logss = logss + "Brutal\n"
		else:
			epa_exp.search(reroot = False, strategy = "H0")
			logss = logss + "H0\n"
		num_spe = epa_exp.count_species(print_log = False)
		
		logss = logss + "	find number specices: " + repr(num_spe) + "\n"
		
		idx = tree.find("leaf")
		if idx >= 0:
			logss = logss + "	the palcement is on a leaf node" + "\n"
		else:
			logss = logss + "	the palcement is on an internal node" + "\n"
		for spe in epa_exp.species_list:
			newalign = gen_alignment2(seq_names = spe, alignment = align)
			
			if idx >= 0:
				morelog = count_and_pick_reads(newalign, nfout1)
				logss = logss + morelog
			else:
				morelog = count_and_pick_reads(newalign, nfout2)
				logss = logss + morelog
			
		sp_log(sfout = nfolder + "spcount.log", logs = logss)


def stas(sfin):
	cnt_read2 = 0 
	cnt_read3 = 0
	cnt_read4 = 0
	cnt_read5 = 0
	ref_map = 0
	ref_map2 = 0 
	ref_map5 = 0
	align = SeqGroup(sequences = sfin)
	for entr in align.get_entries():
		numreads = int(entr[0].split("*")[-1])
		if entr[0].startswith("*R*"):
			ref_map = ref_map + 1
			
		if entr[0].startswith("*R*"):
			if numreads > 4:
				ref_map5 = ref_map5 + 1
				
		if entr[0].startswith("*R*"):
			if numreads > 1:
				ref_map2 = ref_map2 + 1
		
		if numreads > 1:
			cnt_read2 = cnt_read2 + 1
		if numreads > 2:
			cnt_read3 = cnt_read3 + 1
		if numreads > 3:
			cnt_read4 = cnt_read4 + 1
		if numreads > 4:
			cnt_read5 = cnt_read5 + 1
	
	print(">=5 reads OTUs: " + repr(cnt_read5))
	print(">=4 reads OTUs: " + repr(cnt_read4))
	print(">=3 reads OTUs: " + repr(cnt_read3))
	print(">=2 reads OTUs: " + repr(cnt_read2))
	print(">=1 reads OTUs: " + repr(len(align.get_entries())))
	print("maped OTUs: " + repr(ref_map))
	print("maped OTUs >=5: " + repr(ref_map5))
	print("maped OTUs >=2: " + repr(ref_map2))


def random_remove_taxa(falign, num_remove, num_repeat = 10):
	align = SeqGroup(sequences = falign)
	entrs = align.get_entries()
	numseq = len(entrs)
	index = range(numseq)
	
	for i in range(num_repeat):
		newalign = SeqGroup()
		random.shuffle(index)
		idxs = index[num_remove:]
		for idx in idxs:
			newalign.set_seq(entrs[idx][0], entrs[idx][1])
		newalign.write(outfile = falign + "_" + repr(num_remove)+ "_" + repr(i + 1) + ".fasta")


def print_cluster_script(nfolder): #tree
	naligns = glob.glob(nfolder + "*" + ".fasta")
	appd = "/hits/sco/zhangje/biosoup/reduced_ref/"
	for aln in naligns:
		alnname = aln.split("/")[-1]
		print("/home/zhangje/bin/raxmlHPC-PTHREADS-SSE3 -m GTRGAMMA -p 1234 -T 48 -s " + appd + alnname + " -n " + alnname)


def print_cluster_script_EPA(nfolder): #EPA
	ntrees = glob.glob(nfolder + "RAxML_bestTree.*")
	appd = "/hits/sco/zhangje/biosoup/reduced_epa/"
	for tree in ntrees:
		treename = tree.split("/")[-1]
		alnname = str(treename[15:]) + ".combin.fasta"
		print("/home/zhangje/bin/raxmlHPC-PTHREADS-SSE3 -m GTRGAMMA -p 1234 -T 48 -f v -s " + appd + alnname + " -n " + str(treename[32:-6]) + " -r " + appd + treename)


def merge_align(nfolder, talign = "", pref = "origin_ref.fasta_"):
	#fout = open(sfout, "w")
	refalign = open(talign)
	lines = refalign.readlines()
	refalign.close()
	
	naligns = glob.glob(nfolder + pref + "*")
	for aln in naligns:
		newfout = open(aln+".combin.fasta", "w")
		newaln = SeqGroup(sequences = aln)
		for entr in newaln.get_entries():
			newfout.write(">" + entr[0] + "\n")
			newfout.write(entr[1] + "\n")
		
		for line in lines:
			newfout.write(line)
		#newfout.write("\n")
		newfout.close()


def count_reads(nfolder, pref = "bds_leaf_"):
	cnt = 0
	naligns = glob.glob(nfolder + pref + "*")
	for aln in naligns:
		a = SeqGroup(sequences = aln)
		for ent in a.get_entries():
			name = ent[0]
			if name == "root_ref":
				pass
			elif name.startswith("*R*"):
				pass
			else:
				numread = int(name.split("*")[-1])
				cnt = cnt + numread
	print cnt


if __name__ == "__main__":
	#remove_seq_len_smaller_than("/home/zhangje/GIT/gGMYC/biosoup/production/454input.fna", min_l = 100)
	#collapse_identical_seqs("/home/zhangje/GIT/gGMYC/biosoup/production/454.stock.fasta")
	#parse_HMM("/home/zhangje/GIT/gGMYC/biosoup/production/454.stock")
	
	#chimera_removal(nuseach = "/home/zhangje/GIT/gGMYC/biosoup/production/454.chimera.uchime", nalign = "/home/zhangje/GIT/gGMYC/biosoup/production/454.stock.fasta.collapse.fasta", nout = "/home/zhangje/GIT/gGMYC/biosoup/production/454.chimerafree2.fasta", chimeraout = "/home/zhangje/GIT/gGMYC/biosoup/production/454.chimera2.fasta")
	#extract_placement2(nfin_place = "/home/zhangje/GIT/gGMYC/biosoup/production/EPA/454.placement.jplace", nfin_aln = "/home/zhangje/GIT/gGMYC/biosoup/production/EPA/454.epainput.fasta", nfout = "/home/zhangje/GIT/gGMYC/biosoup/production/EPA/bds")
	#raxml_g(nfolder = "/home/zhangje/biosoup/production/EPA/", nfout = "/home/zhangje/biosoup/production/EPA/ttt", nref_align = "/home/zhangje/biosoup/production/EPA/ref_outs_547.fas", suf = "mttree")
	#print estimate_ref_exp_rate("/home/zhangje/GIT/gGMYC/biosoup/production/Test/ref_out547.tre")
	#otu_picking(nfolder = "/home/zhangje/biosoup/production/EPA3/OTU/", nfout1 = "/home/zhangje/biosoup/production/EPA3/OTU/bds_leaf_picked_otus.fasta", nfout2 = "/home/zhangje/biosoup/production/EPA3/OTU/bds_inode_picked_otus.fasta", nref_tree = "/home/zhangje/biosoup/production/chimera/ref_out547.tre", n_align = "/home/zhangje/biosoup/production/chimera/454.epainput.chimerafree.fasta", suf = "subtree")
	#stas(sfin = "/home/zhangje/GIT/gGMYC/bds_leaf_picked_otus.fasta")
	extract_placement3(nfin_place = "/home/zhangje/GIT/biosoup/production/EPA4/454.cfree.jpalce", nfin_aln = "/home/zhangje/GIT/biosoup/production/EPA4/454.epainput.chimerafree.fasta", nfout = "/home/zhangje/GIT/biosoup/production/EPA4/bds", logfile = "/home/zhangje/GIT/biosoup/production/EPA4/spcount.log")
	#raxml(nfolder = "/home/zhangje/biosoup/production/EPA3/", nfout = "/home/zhangje/biosoup/production/EPA3/ttt", suf = "fasta")
	#count_reads("/home/zhangje/GIT/gGMYC/biosoup/production/EPA3/")
	#random_remove_taxa(falign = "/home/zhangje/GIT/gGMYC/biosoup/production/reduce_reftree/origin_ref.fasta", num_remove = 270, num_repeat = 10)
	#print_cluster_script(nfolder = "/home/zhangje/GIT/gGMYC/biosoup/production/reduce_reftree/")
	#print_cluster_script_EPA(nfolder = "/home/zhangje/GIT/gGMYC/biosoup/production/reduce_reftree/")
	#merge_align(nfolder = "/home/zhangje/GIT/gGMYC/biosoup/production/reduce_reftree/", talign = "/home/zhangje/GIT/gGMYC/biosoup/production/reduce_reftree/454.chimerafree2.fasta", pref = "origin_ref.fasta_")
	
