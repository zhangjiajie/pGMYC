#! /usr/bin/env python
import json
import os 
import glob
from ete2 import Tree
from ete2 import SeqGroup
from subprocess import call
tempfiles=[]
alpha = 0.01

def genalignment(seq_names = [], alignment = SeqGroup(), outputfile = "epa_parser_alignments.out"):
	"""generate alignment from the input taxa name list - seq_name, and SeqGroup - alignment"""
	newalign = SeqGroup(format='phylip')
	for taxa in seq_names:
		seq = alignment.get_seq(taxa)
		newalign.set_seq(taxa, seq)
	newalign.write(format='iphylip_relaxed', outfile = outputfile)
	tempfiles.append(outputfile)
	tree_string = call_raxml_tree(alignment_file_name = outputfile, outputname = outputfile)
	alltrees = open(outputfile+"_combined.trees","w")
	#tempfiles.append(outputfile+"_combined.trees")
	alltrees.write(tree_string+"\n")
	return newalign 

	
def gen_random_tree_from_alignment(alignment, num_rand_trees=1, outputfile = "epa_parser_randtree.out", alignment_file = None):
	alltrees = open(outputfile+"_combined.trees","a")
	taxa_names = []	
	for seq in alignment:
		taxa_names.append(seq[0])
	for i in xrange(num_rand_trees):
		t = Tree()
		t.populate(len(taxa_names), names_library=taxa_names)
		t.write(outfile=outputfile+"_"+repr(i)+"_random.tre")
		tree_string = call_raxml_branchlength_optimization(alignment_file_name = alignment_file, tree_file_name = outputfile+"_"+repr(i)+"_random.tre", outputname = outputfile+"_"+repr(i)+"_raxml.tre")
		alltrees.write(tree_string+"\n")
	alltrees.close()
	
def gen_persite_likelihood(alignment_file, trees_file):
	call(["raxmlHPC-SSE3","-m","GTRGAMMA","-s",alignment_file,"-f","g","-n",trees_file,"-z", trees_file])
	tempfiles.append("RAxML_perSiteLLs."+trees_file)
		
    
def call_raxml_branchlength_optimization(alignment_file_name, tree_file_name, outputname):
	call(["raxmlHPC-SSE3","-m","GTRGAMMA","-s",alignment_file_name,"-f","e","-n",outputname,"-t", tree_file_name])
	outtree = open ("RAxML_result."+outputname,"r")
	tree_string = outtree.readline()
	outtree.close()
	return tree_string
	
def call_raxml_tree(alignment_file_name, outputname):
	call(["raxmlHPC-SSE3","-m","GTRGAMMA","-s",alignment_file_name,"-n",outputname,"-p", "1234"])
	outtree = open ("RAxML_bestTree."+outputname,"r")
	tree_string = outtree.readline()
	outtree.close()
	return tree_string

def call_consel(persite_file):
	call(["./makermt","-f","--puzzle",persite_file, "tmpconsel"])
	call(["./consel", "tmpconsel"])
	fnames = persite_file.split(".")
	call(["./catpv -s 6 "+"tmpconsel"+" > "+fnames[1]+".conselout"], shell = True)
	fin = open(fnames[1]+".conselout")
	fout = open("Result.txt","a")
	line_s = fin.readline()
	line_s = fin.readline()
	line_s = fin.readline()
	line_s = fin.readline()
	line_s = fin.readline()
	all_trees_sh = []
	rank = 0
	while line_s != "":
		line_s = line_s.strip()
		if line_s != "":
			items = line_s.split()
			sh = items[8]
			rank = items[1]
			#if sh < alpha:
			all_trees_sh.append(sh)
			
		line_s = fin.readline()
		
	ration_at_alpha = 1.0;
	
	for i,sh in enumerate(all_trees_sh):
		if float(sh) < alpha:
			ration_at_alpha = float(i)/float(len(all_trees_sh) - 1)
			break
	
	num_trees = len(all_trees_sh) - 1
	idx_minp = int(num_trees/2) + 1 
	fin.close()
	fout.write(fnames[1]+" "+repr(ration_at_alpha)+" "+repr(all_trees_sh[idx_minp])+"\n")
	fout.close()

		
def clean_up():
	tempfiles.append("tmpconsel.rmt")
	tempfiles.append("tmpconsel.vt")
	tempfiles.append("tmpconsel.pv")
	for jkfile in tempfiles:
		os.remove(jkfile)
	jk1 = glob.glob("RAxML_binaryModelParameters.*")
	jk2 = glob.glob("RAxML_info.*")
	jk3 = glob.glob("RAxML_log.*")
	jk4 = glob.glob("RAxML_parsimonyTree.*")
	jk5 = glob.glob("RAxML_result.*")
	jk6 = glob.glob("*random*")
	jk7 = glob.glob("RAxML_bestTree.*")

	#jk = []
	#jk1.extend(jk2).extend(jk3).extend(jk4).extend(jk5).extend(jk6)
	for jkfile in jk1:
		os.remove(jkfile)
	for jkfile in jk2:
		os.remove(jkfile)
	for jkfile in jk3:
		os.remove(jkfile)
	for jkfile in jk4:
		os.remove(jkfile)
	for jkfile in jk5:
		os.remove(jkfile)
	for jkfile in jk6:
		os.remove(jkfile)
	for jkfile in jk7:
		os.remove(jkfile)
		
	jk10 = glob.glob("*.reduced")	
	for jkfile in jk10:
		os.remove(jkfile)

def main():
	input_aln_file = "test.phy"
	placemap = {}
	jsondata = open ("RAxML_portableTree.epa.jplace")
	align_orgin = SeqGroup(sequences = "test.phy", format="phylip_relaxed")
	data = json.load(jsondata)
	placements = data["placements"]
	for placement in placements:
		edges = placement["p"]
		curredge = edges[0][0]
		placemap[curredge] = placemap.get(curredge, [])	
	for placement in placements:
		edges = placement["p"]
		taxa_names = placement["n"]
		curredge = edges[0][0]
		a = placemap[curredge] 
		a.extend(taxa_names)
		placemap[curredge]  = a
		#print (placemap[curredge])
	groups = placemap.items()
	for i,item in enumerate(groups):
		seqset_name = item[0]
		seqset = item[1]
		if len(seqset)>3:
			newalign = genalignment(seq_names = seqset, alignment = align_orgin, outputfile = repr(seqset_name)+"_"+repr(i)+".phy")
			gen_random_tree_from_alignment(alignment = newalign, num_rand_trees=10, outputfile = repr(seqset_name)+"_"+repr(i)+".phy", alignment_file = repr(seqset_name)+"_"+repr(i)+".phy")
			gen_persite_likelihood(alignment_file = repr(seqset_name)+"_"+repr(i)+".phy", trees_file = repr(seqset_name)+"_"+repr(i)+".phy_combined.trees")
			call_consel("RAxML_perSiteLLs."+repr(seqset_name)+"_"+repr(i)+".phy_combined.trees")
		#print(seqset)
	clean_up()
	


main()

