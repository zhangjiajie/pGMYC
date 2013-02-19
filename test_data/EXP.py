#! /usr/bin/env python
import sys
import math
import collections
import Clustering
from ete2 import Tree, TreeStyle, TextFace, SeqGroup, NodeStyle
from scipy.optimize import fmin
from collections import deque
from scipy import stats
from numpy import array

#TODO: root at each node; branching and bound using broad first search; multiple rates and clustering 

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
		if prob == 0:
			return float("-inf") 
		else:
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
	"""init(), search() and count_species()"""
	def __init__(self, tree, sp_rate = 0, fix_sp_rate = False, max_iters = 20000, min_br = 0.0001):
		self.min_brl = min_br
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
		self.max_num_search = max_iters


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


	def H0(self, reroot = True):
		self.H1(reroot)
		self.H2(reroot = False)
		self.H3(reroot = False)


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
		
		if max_logl > self.max_logl:
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
			
		if max_logl > self.max_logl:
			self.max_logl = max_logl
			self.max_setting = max_setting


	def H3(self, reroot = True):
		if reroot:
			self.re_rooting()
		sorted_node_list = self.tree.get_descendants()
		sorted_node_list.sort(key=self.__compare_node)
		sorted_node_list.reverse()
		sorted_br = []
		for node in sorted_node_list:
			sorted_br.append(node.dist)
		maxlogl = float("-inf") 
		maxidx = -1
		for i in range(len(sorted_node_list))[1:]:
			l1 = sorted_br[0:i]
			l2 = sorted_br[i:]
			e1 = exp_distribution(l1)
			e2 = exp_distribution(l2)
			logl = e1.sum_log_l() + e2.sum_log_l()
			if logl > maxlogl:
				maxidx = i
				maxlogl = logl
		
		target_nodes = sorted_node_list[0:maxidx]
		
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
		target_node_cnt = 0
		while contin_flag:
			curr_max_logl = float("-inf") 
			curr_max_setting = None
			contin_flag = False
			unchanged_flag = True
			for node in last_setting.active_nodes:
				if node.is_leaf():
					pass
				else:
					contin_flag = True 
					childs = node.get_children()
					sp_nodes = []
					flag = False
					for child in childs:
						if child in target_nodes:
							flag = True
							#target_nodes.remove(child)
					if flag:
						unchanged_flag = False
						for child in childs:
							sp_nodes.append(child)
						for nod in last_setting.spe_nodes:
							sp_nodes.append(nod)
						new_sp_setting = species_setting(spe_nodes = sp_nodes, root = self.tree, sp_rate = self.fix_spe, fix_sp_rate = self.fix_spe_rate)
						logl = new_sp_setting.get_log_l()
						if logl > curr_max_logl:
							curr_max_logl = logl
							curr_max_setting = new_sp_setting
			if not unchanged_flag:
				target_node_cnt = target_node_cnt + 1
				if curr_max_logl > max_logl:
					max_setting = curr_max_setting
					max_logl = curr_max_logl
				last_setting = curr_max_setting
			
			if len(target_nodes) == target_node_cnt:
				contin_flag = False
			if contin_flag and unchanged_flag and last_setting!= None:
				for node in last_setting.active_nodes:
					if node.is_leaf():
						pass
					else:
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
				
		if max_logl > self.max_logl:
			self.max_logl = max_logl
			self.max_setting = max_setting


	def Brutal(self, reroot = False):
		if reroot:
			self.re_rooting()
		first_node_list = []
		first_node_list.append(self.tree)
		first_childs = self.tree.get_children()
		for child in first_childs:
			first_node_list.append(child)
		num_s = self.comp_num_comb()
		if num_s > self.max_num_search:
			print("Too many search iterations: " + repr(num_s) + ", using H0 instead!!!")
			self.H0(reroot = False)
		else:
			first_setting = species_setting(spe_nodes = first_node_list, root = self.tree, sp_rate = self.fix_spe, fix_sp_rate = self.fix_spe_rate)
			self.next(first_setting)


	def search(self, strategy = "H1", reroot = False):
		if strategy == "H1":
			self.H1(reroot)
		elif strategy == "H2":
			self.H2(reroot)
		elif strategy == "H3":
			self.H3(reroot)
		elif strategy == "Brutal":
			self.Brutal(reroot)
		else:# strategy == "H0":
			self.H0(reroot)


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



if __name__ == "__main__":
	
	if len(sys.argv) < 3: 
		print("usage: ./EXP.py -t <tree_of_life.tre> -m <H1/H2/H3/Brutal>  -r(reroot)  -s(show)  -c <scale> -maxiters <max num of brutal search(20000)> -minbr <minimal branch length accepted> -sprate <fix speciation rate>")
		sys.exit()
	
	stree = ""
	sreroot = False
	sstrategy = "H1"
	sprint = False 
	sshow = False
	sscale = 500
	max_iter = 20000
	min_brl = 0.0001
	spe_rate = -1.0
	
	for i in range(len(sys.argv)):
		if sys.argv[i] == "-t":
			i = i + 1
			stree = sys.argv[i]
		elif sys.argv[i] == "-m":
			i = i + 1
			sstrategy = sys.argv[i]
		elif sys.argv[i] == "-r":
			sreroot = True
		elif sys.argv[i] == "-s":
			sshow = True
		elif sys.argv[i] == "-c":
			i = i + 1
			sscale = int(sys.argv[i])
		elif sys.argv[i] == "-p":
			sprint = True 
		elif sys.argv[i] == "-maxiters":
			i = i + 1
			max_iter = int(sys.argv[i])
		elif sys.argv[i] == "-minbr":
			i = i + 1
			min_brl = float(sys.argv[i])
		elif sys.argv[i] == "-sprate":
			i = i + 1
			spe_rate = float(sys.argv[i])
		
	
	if stree == "":
		print("usage: ./EXP.py -t <tree_of_life.tre> -m <H1/H2/H3/Brutal>  -r(reroot)  -s(show)  -c <scale> -maxiters <max num of brutal search(20000)> -minbr <minimal branch length accepted> -sprate <fix speciation rate>")
		sys.exit()
	
	me = None 
	if spe_rate <= 0:
		me = exponential_mixture(tree= stree, max_iters = max_iter, min_br = min_brl )
	else:
		me = exponential_mixture(tree= stree, max_iters = max_iter, min_br = min_brl, sp_rate = spe_rate, fix_sp_rate = True)
	
	me.search(reroot = sreroot, strategy = sstrategy)
	
	if sprint:
		me.print_species()
	print("Number of species:" + repr(me.count_species()))
	
	if sshow:
		me.showTree(scale = sscale)

