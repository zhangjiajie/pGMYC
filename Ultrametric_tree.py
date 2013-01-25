#! /usr/bin/env python
import sys
import os
import json
import operator
import GMYC
from ete2 import Tree, TreeStyle, TextFace, SeqGroup


class um_tree:
	def __init__(self, tree):
		self.tree = Tree(tree, format = 1)
		self.tree.add_feature("age", 0)
		self.nodes = self.tree.get_descendants()
		internal_node = []
		cnt = 0
		for n in self.nodes:
			node_age = n.get_distance(self.tree)
			n.add_feature("age", node_age)
			if not n.is_leaf():
				n.add_feature("id", cnt)
				cnt = cnt + 1
				internal_node.append(n)
		self.nodes = internal_node
		one_leaf = self.tree.get_farthest_node()[0]
		one_leaf.add_feature("id", cnt+1)
		if one_leaf.is_leaf():
			self.nodes.append(one_leaf)
		self.nodes.sort(key=self.__compare_node)
		
	def __compare_node(self, node):
		return node.age
		
	def get_waiting_times(self, threshold_node = None, threshold_node_idx = 0):
		wt_list = []
		reach_t = False
		curr_age = 0.0
		curr_spe = 2
		curr_num_coa = 0
		coa_roots = []
		min_brl = 1000
		num_spe = -1
		
		if threshold_node == None:
			threshold_node = self.nodes[threshold_node_idx]
		
		for node in self.nodes:
			wt = None
			times = node.age - curr_age
			if times >= 0:
				if times < min_brl and times > 0:
					min_brl = times
				curr_age = node.age
				assert curr_spe >=0
				 
				if reach_t:
					fnode = node.up
					coa_root = None
					while not fnode.is_root():
						for coa_r in coa_roots:
							if coa_r.id == fnode.id:
								coa_root = coa_r
								break
						if coa_root!=None:
							break
						else:
							fnode = fnode.up
					wt = GMYC.waiting_time(length = times, num_coas =curr_num_coa, num_lines = curr_spe)
					for coa_r in coa_roots:
						coa = GMYC.coalescent(num_individual = coa_r.curr_n)
						wt.coas.add_coalescent(coa)
					if coa_root == None: #here can be modified to use multiple T
						curr_spe = curr_spe - 1
						curr_num_coa = curr_num_coa + 1
						node.add_feature("curr_n", 2)
						coa_roots.append(node)
					else:
						curr_n = coa_root.curr_n
						coa_root.add_feature("curr_n", curr_n + 1)
						
				else:
					if node.id == threshold_node.id:
						reach_t = True
						wt = GMYC.waiting_time(length = times, num_coas = 0, num_lines = curr_spe)
						num_spe = curr_spe
						curr_spe = curr_spe - 1
						curr_num_coa = curr_num_coa + 1
						node.add_feature("curr_n", 2)
						coa_roots.append(node)
					else:
						wt = GMYC.waiting_time(length = times, num_coas = 0, num_lines = curr_spe)
						curr_spe = curr_spe + 1
				if times > 0:
					wt_list.append(wt)
		
		if min_brl < 0.0001:
			min_brl = 0.0001
		
		for wt in wt_list:
			pass
			#wt.length = wt.length/min_brl
			#print(wt)
		
		return wt_list, num_spe


if __name__ == "__main__":
	#print("main function")
	#t = Tree("test.tree.tre", format = 1)
	#print(t)
	#t.show()
	#lvs = t.get_leaves()
	#for leaf in lvs:
		 #print (leaf.get_distance(t))
		 
	#nodes = t.get_descendants()
	#for n in nodes:
	#	print (n.dist)
	#print t.dist
	ut =um_tree(tree = "test.tree.tre")
	#ut =um_tree(tree = "2mtree.tre")
	wl, n = ut.get_waiting_times(threshold_node_idx = 0)
	print(n)
	print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
	wl, n = ut.get_waiting_times(threshold_node_idx = 1)
	print(n)
	print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
	wl, n = ut.get_waiting_times(threshold_node_idx = 2)
	print(n)
	print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
	wl, n = ut.get_waiting_times(threshold_node_idx = 3)
	print(n)
	print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
	wl, n = ut.get_waiting_times(threshold_node_idx = 4)
	print(n)
	print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
	wl, n = ut.get_waiting_times(threshold_node_idx = 5)
	print(n)
	print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
