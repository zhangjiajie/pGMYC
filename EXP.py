#! /usr/bin/env python
import sys
#import os
#import json
#import operator
#import Ultrametric_tree
import math
#import random
from ete2 import Tree, TreeStyle, TextFace, SeqGroup, NodeStyle
#from subprocess import call
from scipy.optimize import fmin
import collections
from collections import deque
from scipy import stats

#TODO: implement H3; root at each node; search with all 3 H; brutal search when search < max 

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
		self.logl = 0
		self.spe_list = []
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
		self.logl = logl
		return logl
		
	def count_species(self):
		if len(self.spe_list) != 0:
			return len(self.spe_list), self.spe_list
		else:
			for node in self.active_nodes:
				one_spe = []
				if node.is_leaf():
					one_spe.append(node.name)
				else:
					one_spe.extend(node.get_leaf_names())
				self.spe_list.append(one_spe)
			return len(self.spe_list), self.spe_list

class exponential_mixture:
	def __init__(self, tree, sp_rate = 0, fix_sp_rate = False):
		self.min_brl = 0.0001
		self.tree = Tree(tree, format = 1)
		self.tree.dist = 0.0
		self.fix_spe_rate = fix_sp_rate
		self.fix_spe = sp_rate
		self.max_logl = float("-inf") 
		self.max_setting = None
		self.null_logl = 0.0
		self.null_model()
		self.species_list = None
		self.counter = 0
		self.setting_set = set([])
		self.max_num_search = 10000


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
		logl = sp_setting.get_log_l()
		if logl > self.max_logl:
			self.max_logl = logl
			self.max_setting = sp_setting
		for node in sp_setting.active_nodes:
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


	def H1(self, reroot = True):
		if reroot:
			self.re_rooting()
			
		#self.init_tree()
		sorted_node_list = self.tree.get_descendants()
		sorted_node_list.sort(key=self.__compare_node)
		sorted_node_list.reverse()
		
		first_node_list = []
		first_node_list.append(self.tree)
		first_childs = self.tree.get_children()
		for child in first_childs:
			first_node_list.append(child)
		first_setting = species_setting(spe_nodes = first_node_list, root = self.tree, sp_rate = self.fix_spe, fix_sp_rate = self.fix_spe_rate)
		last_setting = first_setting
		max_logl = last_setting.get_log_l()
		max_setting = last_setting
		
		for node in sorted_node_list:
			if node not in last_setting.spe_nodes:
				curr_sp_nodes = []
				for nod in last_setting.spe_nodes:
					curr_sp_nodes.append(nod)
				
				chosen_branching_node = node.up #find the father of this new node
				if chosen_branching_node in last_setting.spe_nodes:
					for nod in chosen_branching_node.get_children():
						if nod not in curr_sp_nodes:
							curr_sp_nodes.append(nod)
				else:
					for nod in chosen_branching_node.get_children():
						if nod not in curr_sp_nodes:
							curr_sp_nodes.append(nod)
					while not chosen_branching_node.is_root():
						chosen_branching_node = chosen_branching_node.up
						for nod in chosen_branching_node.get_children():
							if nod not in curr_sp_nodes:
								curr_sp_nodes.append(nod)
						if chosen_branching_node in last_setting.spe_nodes:
							break
				new_setting = species_setting(spe_nodes = curr_sp_nodes, root = self.tree, sp_rate = self.fix_spe, fix_sp_rate = self.fix_spe_rate)
				new_logl = new_setting.get_log_l()
				if new_logl> max_logl:
					max_logl = new_logl
					max_setting = new_setting 
				last_setting = new_setting
				
			else:
				"""node already is a speciation node, do nothing"""
				pass
			
		self.max_logl = max_logl
		self.max_setting = max_setting


	def H2(self, reroot = True):
		"""Greedy"""
		if reroot:
			self.re_rooting()
			
		#self.init_tree()
		sorted_node_list = self.tree.get_descendants()
		sorted_node_list.sort(key=self.__compare_node)
		sorted_node_list.reverse()
		
		first_node_list = []
		first_node_list.append(self.tree)
		first_childs = self.tree.get_children()
		for child in first_childs:
			first_node_list.append(child)
		first_setting = species_setting(spe_nodes = first_node_list, root = self.tree, sp_rate = self.fix_spe, fix_sp_rate = self.fix_spe_rate)
		last_setting = first_setting
		max_logl = last_setting.get_log_l()
		max_setting = last_setting
		contin_flag = True 
		
		
		while contin_flag:
			curr_max_logl = float("-inf") 
			curr_max_setting = None
			contin_flag = False
			for node in last_setting.active_nodes:
				if node.is_leaf():
					pass
				else:
					contin_flag = True 
					childs = node.get_children()
					sp_nodes = []
					for child in childs:
						sp_nodes.append(child)
					for nod in last_setting.spe_nodes:
						sp_nodes.append(nod)
					new_sp_setting = species_setting(spe_nodes = sp_nodes, root = self.tree, sp_rate = self.fix_spe, fix_sp_rate = self.fix_spe_rate)
					logl = new_sp_setting.get_log_l()
					if logl > curr_max_logl:
						curr_max_logl = logl
						curr_max_setting = new_sp_setting
			
			if curr_max_logl > max_logl:
				max_setting = curr_max_setting
				max_logl = curr_max_logl
			
			last_setting = curr_max_setting
			
		self.max_logl = max_logl
		self.max_setting = max_setting



	def H3(self, reroot = True):
		pass


	def Brutal(self, reroot = False): #, reroot_node = None):
		if reroot:
			self.re_rooting()
		first_node_list = []
		first_node_list.append(self.tree)
		first_childs = self.tree.get_children()
		for child in first_childs:
			first_node_list.append(child)
		
		first_setting = species_setting(spe_nodes = first_node_list, root = self.tree, sp_rate = self.fix_spe, fix_sp_rate = self.fix_spe_rate)
		self.next(first_setting)


	def search(self, strategy = "H1", reroot = False):
		if strategy == "H1":
			self.H1(reroot)
		elif strategy == "H2":
			self.H2(reroot)
		elif strategy == "Brutal":
			self.Brutal(reroot)


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
			#when there are more than one speces, should do sth to make the results better?
			#self.fix_spe_rate = False
			#self.max_logl = float("-inf") 
			#self.search()
			return num_sp
		else:
			self.species_list = []
			self.species_list.append(self.tree.get_leaf_names()) 
			return 1


	def print_species(self):
		cnt = 1
		for sp in self.species_list:
			print("Species " + repr(cnt) + ":")
			for leaf in sp:
				print("          " + leaf)
			cnt = cnt + 1


	def showTree(self, scale = 500):
		style0 = NodeStyle()
		style0["fgcolor"] = "#000000"
		#style2["shape"] = "circle"
		style0["vt_line_color"] = "#0000aa"
		style0["hz_line_color"] = "#0000aa"
		style0["vt_line_width"] = 2
		style0["hz_line_width"] = 2
		style0["vt_line_type"] = 0 # 0 solid, 1 dashed, 2 dotted
		style0["hz_line_type"] = 0
		style0["size"] = 0
		
		for node in self.tree.get_descendants():
			node.set_style(style0)
			node.img_style["size"] = 0
		self.tree.set_style(style0)
		self.tree.img_style["size"] = 0
		
		
		style1 = NodeStyle()
		style1["fgcolor"] = "#000000"
		#style2["shape"] = "circle"
		style1["vt_line_color"] = "#ff0000"
		style1["hz_line_color"] = "#0000aa"
		style1["vt_line_width"] = 2
		style1["hz_line_width"] = 2
		style1["vt_line_type"] = 0 # 0 solid, 1 dashed, 2 dotted
		style1["hz_line_type"] = 0
		style1["size"] = 0
		
		style2 = NodeStyle()
		style2["fgcolor"] = "#0f0f0f"
		#style2["shape"] = "circle"
		style2["vt_line_color"] = "#ff0000"
		style2["hz_line_color"] = "#ff0000"
		style2["vt_line_width"] = 2
		style2["hz_line_width"] = 2
		style2["vt_line_type"] = 0 # 0 solid, 1 dashed, 2 dotted
		style2["hz_line_type"] = 0
		style2["size"] = 0
		
		for node in self.max_setting.active_nodes:
			node.set_style(style1)
			node.img_style["size"] = 0
			for des in node.get_descendants():
				des.set_style(style2)
				des.img_style["size"] = 0
		ts = TreeStyle()
		#ts.show_leaf_name = True
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
	if len(sys.argv) < 5: 
		print("usage: ./EXP.py  <tree_of_life.tre> <H1/Brutal> <-r/-n (reroot or not)>  <-s/-n (show or not)>  <scale>")
		sys.exit()
	me = exponential_mixture(sys.argv[1])
	
	if sys.argv[3] == "-r":
		me.search(reroot = True, strategy = sys.argv[2])
		print("Number of species:" + repr(me.count_species()))
		me.print_species()
	else:
		me.search(reroot = False, strategy = sys.argv[2])
		print("Number of species:" + repr(me.count_species()))
		me.print_species()
	
	if sys.argv[4] == "-s":
		if (len(sys.argv) > 5):
			me.showTree(scale = int(sys.argv[5]))
		else:
			me.showTree(scale = 500)
