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
		self.nodes.sort(key=self.__compare_node)
		#for n in self.nodes:
		#	print(repr(n.age) + repr(n.is_leaf()))
		
			
	def __compare_node(self, node):
		return node.age	
		
	def get_waiting_times(threshold_node):
		reach_t = False
		curr_age = 0.0
		curr_spe = 2
		curr_num_coa = 0
		num_spe = -1
		coa_roots = []
		
		for node in self.nodes:
			if reach_t:
				fnode = node.up
				coa_root = None
				while not fnode.is_root:
					for coa_r in coa_roots:
						if coa_r.id == fnode.id:
							coa_root = coa_r
							break
					if coa_root!=None:
						break
					else:
						fnode = fnode.up
				if coa_root == None: #here can be modified to use multiple T
					
				else:
					
					
			else:
				if node.id == threshold_node.id:
					reach_t = True
					num_spe = curr_spe
					times = node.age - curr_age
					curr_age = node.age
					wt = GMYC.waiting_time(length = time, num_coas = 0, num_lines = curr_spe)
					curr_spe = curr_spe - 1
					curr_num_coa = curr_num_coa + 1
					node.add_feature("curr_n", 2)
					coa_roots.append(node)
				else:
					times = node.age - curr_age
					curr_age = node.age
					wt = GMYC.waiting_time(length = time, num_coas = 0, num_lines = curr_spe)
					curr_spe = curr_spe + 1
					
				
		
			
			
	
	
	
if __name__ == "__main__":
	print("main function")
	#t = Tree("2mtree.tre", format = 1)
	#lvs = t.get_leaves()
	#for leaf in lvs:
		 #print (leaf.get_distance(t))
		 
	#nodes = t.get_descendants()
	#for n in nodes:
	#	print (n.dist)
	#print t.dist
	ut =um_tree(tree = "2mtree.tre")
