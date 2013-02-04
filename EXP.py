#! /usr/bin/env python
import sys
import os
import json
import operator
import Ultrametric_tree
import math
import random
from ete2 import Tree, TreeStyle, TextFace, SeqGroup
from subprocess import call
from scipy.optimize import fmin
import collections
from scipy import stats

class lh_ratio_test:
	def __init__(self, null_llh, llh, df):
		self.lr = 2.0 * (llh - null_llh)
		self.p = 1 - stats.chi2.cdf(self.lr, df)
	
	def get_p_value(self):
		return self.p

class exp_distribution:
	def __init__(self, data):
		self.data = data
		self.rate = 0.0
		self.estimate_rate()
		
	def __str__(self):
		return "Exponential distribution with rate = " + repr(self.rate)
	
	def estimate_rate(self):
		sumbr = 0.0
		numbr = len(self.data)
		for br in self.data:
			sumbr = sumbr + br
		if sumbr == 0:
			self.rate = 0
		else:
			self.rate = float(numbr) / sumbr
		
	def log_l(self, x):
		prob = self.rate * math.exp (-1.0 * self.rate * x)
		return math.log(prob)
		
	def sum_log_l(self):
		s = 0.0
		for br in self.data:
			s = s + self.log_l(br)
		return s 

class mix_exp:
	def __init__(self, tree):
		self.tree = Tree(tree, format = 1)
		self.spe_br = []
		self.coa_br = []
		self.null_logl = 0
		self.max_logl = float("-inf") 
		self.num_spe = -1
		self.q = collections.deque()
		self.rate1 = None 
		self.rate2 = None
		self.species_list = []
		self.min_brl = 0.0005
	
	def init_tree(self):
		self.coa_br = []
		for node in self.tree.get_descendants():
			node.add_feature("t", "coa")
			if node.dist > self.min_brl:
				self.coa_br.append(node.dist)
		e1 = exp_distribution(self.coa_br)
		self.null_logl = e1.sum_log_l()
	
	def set_model_data(self):
		self.spe_br = []
		self.coa_br = []
		for node in self.tree.get_descendants():
			if node.t == "coa":
				if node.dist > self.min_brl:
					self.coa_br.append(node.dist)
			else:
				if node.dist > self.min_brl:
					self.spe_br.append(node.dist)
		e1 = exp_distribution(self.coa_br)
		e2 = exp_distribution(self.spe_br)
		self.rate1 = e1.rate
		self.rate2 = e2.rate
		logl = e1.sum_log_l() + e2.sum_log_l()
		return logl
		
	
	def __compare_node(self, node):
		return node.dist
	
	def re_rooting(self):
		node_list = self.tree.get_descendants()
		node_list.sort(key=self.__compare_node)
		node_list.reverse()
		rootnode = node_list[0]
		self.tree.set_outgroup(rootnode)
	
	def search(self):
		self.re_rooting()
		
		self.init_tree()
		node_list = self.tree.get_descendants()
		node_list.sort(key=self.__compare_node)
		node_list.reverse()
		maxnode = None
		cnt = 0
		for node in node_list:
			node.add_feature("t", "spe")
			currnode = node
			while not node.up.is_root():
				node = node.up
				if cnt == 0:
					cnt = cnt + 1
					break
				if node.t == "spe":
					break
				else:
					node.add_feature("t", "spe")
			curr_logl = self.set_model_data()
			if curr_logl >= self.max_logl:
				self.max_logl = curr_logl
				maxnode = currnode
			
		self.init_tree()
		cnt = 0
		for node in node_list:
			if node == maxnode:
				node.add_feature("t", "spe")
				while not node.up.is_root():
					node = node.up
					if cnt == 0:
						cnt = cnt + 1
						break
					if node.t == "spe":
						break
					else:
						node.add_feature("t", "spe")
				break
			else:
				node.add_feature("t", "spe")
				while not node.up.is_root():
					node = node.up
					if cnt == 0:
						cnt = cnt + 1
						break
					if node.t == "spe":
						break
					else:
						node.add_feature("t", "spe")
		self.set_model_data()
		for t in self.tree.get_descendants():
			if t.t == "spe":
				t.add_face(TextFace("SPE"), column=0, position = "branch-right")
		ts = TreeStyle()
		ts.show_leaf_name = True
		ts.scale =  1000 # 1000 pixels per branch length unit
		self.tree.show(tree_style=ts)
			
	
	def count_species(self):
		node_list = self.tree.get_descendants()
		for node in node_list:
			if node.t == "spe":
				childs = node.get_children()
				flag = False
				for child in childs:
					if child.t == "spe":
						flag = True
						break
				if flag:
					for child in childs:
						if child.t == "spe":
							pass
						else:
							self.species_list.append(child.get_leaves())
				else:
					self.species_list.append(node.get_leaves())
		print(self.null_logl)
		print(self.max_logl)
		lhr = lh_ratio_test(self.null_logl, self.max_logl, 1)
		pvalue = lhr.get_p_value()
		print(pvalue)
		if pvalue >= 0.05:
			self.species_list = []
			self.species_list.append(self.tree.get_leaves()) 
		return len(self.species_list)
		
	
	def print_species(self):
		cnt = 1
		for sp in self.species_list:
			print("Species " + repr(cnt) + ":")
			for leaf in sp:
				print("          " + leaf.name)
			cnt = cnt + 1


if __name__ == "__main__":
	#me = mix_exp("2mtree.tre")
	#me = mix_exp("test.tree.tre")
	me = mix_exp("2m2.tre")
	#me = mix_exp("cox1.tre")
	#me = mix_exp("t1.tre")
	#me = mix_exp("1ms1.tre")
	me.search()
	print(me.count_species())
	me.print_species()
