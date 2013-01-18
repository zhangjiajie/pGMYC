#! /usr/bin/env python
import sys
import os
import json
import operator
from ete2 import Tree, TreeStyle, TextFace, SeqGroup
from subprocess import call

class exp_distribution:
    def __init__(self):
		pass 
		
    def __str__(self):
        return "Nothing at all"
    
	def coalescent_prob(r, n , x):
		"""
		r: coalesce rate = 1/2N
		n: number of lineages 
		x: x
		"""
		prob = r * n * (n - 1.0) * math.exp( -1.0 * r * n * (n - 1.0) * x)
		return math.log(prob)
		
	def yule_prob(r, n , x):
		"""
		r: speciation rate
		n: number of lineages 
		x: x
		"""
		prob = r * n * math.exp(-1.0 * r * n * x)
		return math.log(prob)
		
	def gmyc_prob(b, x):
		prob = b * math.exp (-1.0 * b * x)
		return math.log(prob)
		
	
	def sum_log_l(bl, xl):
		s = 0.0
		for i in range(len(bl)):
			s = s + gmyc_prob(bl[i], xl[i])
		return s 
		
	
class waiting_time:
	def __init__(self):
		self.length = 1.0
		self.num_species = 0
		self.coalescents = []
		
		
		
		
		
		
	
		
		

    




if __name__ == "__main__":
        if len(sys.argv) != 6: 
            print("usage: ./ncbi_taxonomy.py <tree_of_life.tre> <id_name.txt> <id_rank.txt> <name_tax.txt> <outputfile>")
            sys.exit()
        t = ncbi_taxa()
        t.init_tax_tree(sys.argv[1], sys.argv[2], sys.argv[3])
        t.extract_sub_tax_tree(sys.argv[4], sys.argv[5])
