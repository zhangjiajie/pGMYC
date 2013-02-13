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
from collections import deque
from scipy import stats

class lh_ratio_test:
	def __init__(self, null_llh, llh, df):
		self.lr = 2.0 * (llh - null_llh)
		self.p = 1 - stats.chi2.cdf(self.lr, df)
	
	def get_p_value(self):
		return self.p

class exp_distribution:
	def __init__(self, data, rate = -1):
		self.data = data
		self.rate = 0.0
		if rate < 0:
			self.estimate_rate()
		else:
			self.rate = rate
		
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

class species_setting:
	def __init__(self, spe_nodes, root, sp_rate = 0, fix_sp_rate = False):
		self.min_brl = 0.0001
		self.spe_rate = sp_rate
		self.fix_spe_rate = fix_sp_rate
		self.spe_nodes = spe_nodes
		self.root = root
		self.all_nodes = root.get_descendants()
		self.all_nodes.append(self.root)
		self.coa_nodes = []
		for node in self.all_nodes:
			if not (node in self.spe_nodes):
				self.coa_nodes.append(node) 
		
		self.active_nodes = []
		for node in self.spe_nodes:
			if node.is_leaf():
				self.active_nodes.append(node)
			else:
				childs = node.get_children()
				flag = False
				for child in childs:
					if not (child in self.spe_nodes):
						flag = True
						break
				if flag:
					self.active_nodes.append(node)
		
		
	def get_log_l(self):
		spe_br = []
		coa_br = []
		for node in self.spe_nodes:
			if node.dist > self.min_brl:
				spe_br.append(node.dist)
		
		for node in self.coa_nodes:
			if node.dist > self.min_brl:
				coa_br.append(node.dist)
				
		e1 = exp_distribution(coa_br)
		e2 = None
		if self.fix_spe_rate:
			e2 = exp_distribution(spe_br, rate = self.spe_rate)
		else:
			e2 = exp_distribution(spe_br)
		self.rate1 = e1.rate
		self.rate2 = e2.rate
		logl = e1.sum_log_l() + e2.sum_log_l()
		return logl
		
	def count_species(self):
		self.spe_list = []
		for node in self.active_nodes:
			one_spe = []
			if node.is_leaf():
				one_spe.append(node)
			else:
				one_spe.extend(node.get_leaves())
				"""
				childs = node.get_children()
				
				for child in childs:
					if not child in self.spe_nodes:
						if child.is_leaf():
							one_spe.append(child)
						else:
							one_spe.extend(child.get_leaves())
				"""
			self.spe_list.append(one_spe)
		return len(self.spe_list), self.spe_list

class exponential_mixture:
	def __init__(self, tree, sp_rate = 0, fix_sp_rate = False):
		self.min_brl = 0.0001
		self.tree = Tree(tree, format = 1)
		self.tree.dist = 0.0
		self.fix_spe_rate = fix_sp_rate
		self.fix_spe = sp_rate
		#self.sp_settings = deque([]) 
		self.max_logl = float("-inf") 
		self.max_setting = None
		self.null_logl = 0.0
		self.null_model()
		self.species_list = None
		self.counter = 0
		self.setting_set = set([])
	
	def null_model(self):
		coa_br = []
		all_nodes = self.tree.get_descendants()
		for node in all_nodes:
			if node.dist > self.min_brl:
				coa_br.append(node.dist)
		e1 = exp_distribution(coa_br)
		self.null_logl = e1.sum_log_l()
		return e1.rate

	def __compare_node(self, node):
		return node.dist
		
	def re_rooting(self):
		node_list = self.tree.get_descendants()
		node_list.sort(key=self.__compare_node)
		node_list.reverse()
		rootnode = node_list[0]
		self.tree.set_outgroup(rootnode)
		self.tree.dist = 0.0
	
	def comp_num_comb(self):
		for node in self.tree.traverse(strategy='postorder'):
			if node.is_leaf():
				node.add_feature("cnt", 1.0)
			else:
				acum = 1.0
				for child in node.get_children():
					acum = acum * child.cnt
				acum = acum + 1.0
				node.add_feature("cnt", acum)
		return self.tree.cnt
	
	def next(self, sp_setting):
		self.setting_set.add(frozenset(sp_setting.spe_nodes))
		self.counter = self.counter + 1
		#print("Search No. :" + repr(self.counter))
		logl = sp_setting.get_log_l()
		#print(logl)
		if logl > self.max_logl:
			self.max_logl = logl
			self.max_setting = sp_setting
			
		#print("Num active nodes:" + repr(len(sp_setting.active_nodes)))
		#print sp_setting.active_nodes
		cnt = 1
		for node in sp_setting.active_nodes:
			#print(cnt)
			cnt = cnt + 1
			if node.is_leaf():
				pass
			else:
				childs = node.get_children()
				sp_nodes = []
				for child in childs:
					sp_nodes.append(child)
				for nod in sp_setting.spe_nodes:
					sp_nodes.append(nod)
				new_sp_setting = species_setting(spe_nodes = sp_nodes, root = sp_setting.root, sp_rate = sp_setting.spe_rate, fix_sp_rate = sp_setting.fix_spe_rate)
				if frozenset(sp_nodes) in self.setting_set:
					pass
				else:
					self.next(new_sp_setting)
		
	def search(self, reroot = False):
		if reroot:
			self.re_rooting()
		first_node_list = []
		first_node_list.append(self.tree)
		first_childs = self.tree.get_children()
		for child in first_childs:
			first_node_list.append(child)
		
		#print(len(first_childs))  
		first_setting = species_setting(spe_nodes = first_node_list, root = self.tree, sp_rate = self.fix_spe, fix_sp_rate = self.fix_spe_rate)
		self.next(first_setting)
		
	
	def count_species(self, print_log = True):
		lhr = lh_ratio_test(self.null_logl, self.max_logl, 1)
		pvalue = lhr.get_p_value()
		if print_log:
			print("Speciation rate:" + repr(self.max_setting.rate2))
			print("Coalesecnt rate:" + repr(self.max_setting.rate1))
			print("Null logl:" + repr(self.null_logl))
			print("MAX logl:" + repr(self.max_logl))
			print("P-value:" + repr(pvalue))
		if pvalue < 0.001:
			num_sp, self.species_list = self.max_setting.count_species()
			self.fix_spe_rate = False
			self.max_logl = float("-inf") 
			self.search()
			return num_sp
		else:
			self.species_list = []
			self.species_list.append(self.tree.get_leaves()) 
			return 1
		
		
	def print_species(self):
		cnt = 1
		for sp in self.species_list:
			print("Species " + repr(cnt) + ":")
			for leaf in sp:
				print("          " + leaf.name)
			cnt = cnt + 1
	
	def showTree(self, scale = 1000):
		for t in self.max_setting.spe_nodes:
			t.add_face(TextFace("SPE"), column=0, position = "branch-right")
		ts = TreeStyle()
		ts.show_leaf_name = True
		ts.scale =  scale # scale pixels per branch length unit
		self.tree.show(tree_style=ts)
		

class mix_exp:
	def __init__(self, tree, sp_rate = 0, coa_rate = 0, fix_sp_rate = False):
		self.tree = Tree(tree, format = 1)
		self.tree.dist = 0
		self.spe_br = []
		self.coa_br = []
		self.null_logl = 0
		self.max_logl = float("-inf") 
		self.num_spe = -1
		self.q = collections.deque()
		self.rate1 = None 
		self.rate2 = None
		self.fix_spe_rate = fix_sp_rate
		self.fix_spe = sp_rate 
		self.fix_coa = coa_rate
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
		return e1.rate
	
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
		e2 = None
		if self.fix_spe_rate:
			e2 = exp_distribution(self.spe_br, rate = self.fix_spe)
		else:
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
	
	def search(self, reroot = True):
		if reroot:
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
		
		
		
	def showTree(self, scale = 1000):
		self.tree.add_face(TextFace("SPE"), column=0, position = "branch-right")
		for t in self.tree.get_descendants():
			if t.t == "spe":
				t.add_face(TextFace("SPE"), column=0, position = "branch-right")
		ts = TreeStyle()
		ts.show_leaf_name = True
		ts.scale =  scale # scale pixels per branch length unit
		self.tree.show(tree_style=ts)
	
	def count_species(self, print_log = True):
		node_list = self.tree.get_descendants()
		leaves = self.tree.get_leaves()
		#self.tree.add_feature("t", "spe")
		#node_list.append(self.tree)
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
					
		rest = []
		added = []
		for sp in self.species_list:
			for leaf in sp:
				added.append(leaf)
		
		for leaf in leaves:
			if not leaf in added:
				rest.append(leaf)
		
		if len(rest) != 0:
			self.species_list.append(rest)

		lhr = lh_ratio_test(self.null_logl, self.max_logl, 1)
		pvalue = lhr.get_p_value()
		if print_log:
			print("Speciation rate:" + repr(self.rate2))
			print("Coalesecnt rate:" + repr(self.rate1))
			print("Null logl:" + repr(self.null_logl))
			print("MAX logl:" + repr(self.max_logl))
			print("P-value:" + repr(pvalue))
		if pvalue >= 0.01:
			self.species_list = []
			self.species_list.append(self.tree.get_leaves()) 
			self.init_tree()
		return len(self.species_list)
		
	
	def print_species(self):
		cnt = 1
		for sp in self.species_list:
			print("Species " + repr(cnt) + ":")
			for leaf in sp:
				print("          " + leaf.name)
			cnt = cnt + 1


if __name__ == "__main__":
	if len(sys.argv) < 4: 
		print("usage: ./EXP.py  <tree_of_life.tre>  <-r/-n (reroot or not)>  <-s/-n (show or not)>  <scale>")
		sys.exit()
	#me = mix_exp(sys.argv[1])
	me = exponential_mixture(sys.argv[1])
	if sys.argv[2] == "-r":
		print(me.comp_num_comb())
		#me.search(reroot = True)
		#print("Number of species:" + repr(me.count_species()))
		#me.print_species()
	else:
		me.search(reroot = False)
		print("Number of species:" + repr(me.count_species()))
		me.print_species()
	
	if sys.argv[3] == "-s":
		if (len(sys.argv) > 4):
			me.showTree(scale = int(sys.argv[4]))
		else:
			me.showTree(scale = 500)
