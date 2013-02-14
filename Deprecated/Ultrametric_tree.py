#! /usr/bin/env python
#import sys
#import os
#import json
#import operator
import GMYC
from ete2 import Tree, TreeStyle, TextFace, SeqGroup
import matplotlib.pyplot as plt


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
		self.species_list = []


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
		
		last_coa_num = 0
		tcnt = 0 
		for node in self.nodes:
			wt = None
			times = node.age - curr_age
			if times >= 0:
				if times < min_brl and times > 0:
					min_brl = times
				curr_age = node.age
				assert curr_spe >=0
				 
				if reach_t:
					if tcnt == 0:
						last_coa_num = 2
					fnode = node.up
					coa_root = None
					
					idx = 0
					while not fnode.is_root():
						idx = 0 
						for coa_r in coa_roots:
							if coa_r.id == fnode.id:
								coa_root = coa_r
								break
							idx = idx + 1
						
						if coa_root!=None:
							break
						else:
							fnode = fnode.up
							
					wt = GMYC.waiting_time(length = times, num_coas =curr_num_coa, num_lines = curr_spe)
					
					for coa_r in coa_roots:
						coa = GMYC.coalescent(num_individual = coa_r.curr_n)
						wt.coas.add_coalescent(coa)
					
					wt.coas.coas_idx = last_coa_num
					wt.num_curr_coa = last_coa_num
					if coa_root == None: #here can be modified to use multiple T
						curr_spe = curr_spe - 1
						curr_num_coa = curr_num_coa + 1
						node.add_feature("curr_n", 2)
						coa_roots.append(node)
						last_coa_num = 2
					else:
						curr_n = coa_root.curr_n
						coa_root.add_feature("curr_n", curr_n + 1)
						last_coa_num = curr_n + 1
					tcnt = tcnt + 1
				else:
					if node.id == threshold_node.id:
						reach_t = True
						tcnt = 0 
						wt = GMYC.waiting_time(length = times, num_coas = 0, num_lines = curr_spe)
						num_spe = curr_spe
						curr_spe = curr_spe - 1
						curr_num_coa = 2
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
			wt.count_num_lines()
		
		self.species_list = []
		all_coa_leaves = []
		for coa_r in coa_roots:
			leaves = coa_r.get_leaves()
			all_coa_leaves.extend(leaves)
			self.species_list.append(leaves)
		
		all_leaves = self.tree.get_leaves()
		for leaf in all_leaves:
			if leaf not in all_coa_leaves:
				self.species_list.append([leaf])
		
		#print(len(self.species_list) == num_spe)
		
		return wt_list, num_spe


	def show(self, wt_list):
		cnt = 1
		for wt in wt_list:
			print("Waitting interval "+ repr(cnt))
			print(wt)
			cnt = cnt + 1
			
	
	def get_species(self):
		sp_list = []
		for sp in self.species_list:
			spe = []
			for taxa in sp:
				spe.append(taxa.name)
			sp_list.append(spe)
		
		all_taxa_name = []
		for leaf in self.tree.get_leaves():
			all_taxa_name.append(leaf.name)
		
		return [all_taxa_name], sp_list
	
	def print_species(self):
		cnt = 1
		for sp in self.species_list:
			print("Species " + repr(cnt) + ":")
			cnt = cnt + 1
			taxas = ""
			for taxa in sp:
				taxas = taxas + taxa.name + ", "
			print("	" + taxas[:-1])
			
	def num_lineages(self, wt_list):
		nl_list = []
		times = []
		last_time = 0.0
		for wt in wt_list:
			nl_list.append(wt.get_num_branches())
			times.append(last_time)
			last_time = wt.length + last_time
		
		plt.plot(times, nl_list)
		plt.ylabel('Number of lineages')
		plt.xlabel('Time')
		plt.savefig("Time_Lines")
		plt.show()
				


if __name__ == "__main__":
	ut =um_tree(tree = "test.tree.tre", r_mode = False)
	wl, n = ut.get_waiting_times(threshold_node_idx = 0)
	print(n)
	ut.show(wl)
	print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
	wl, n = ut.get_waiting_times(threshold_node_idx = 1)
	print(n)
	ut.show(wl)
	print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
	wl, n = ut.get_waiting_times(threshold_node_idx = 2)
	print(n)
	ut.show(wl)
	print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
	wl, n = ut.get_waiting_times(threshold_node_idx = 3)
	print(n)
	ut.show(wl)
	print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
	wl, n = ut.get_waiting_times(threshold_node_idx = 4)
	print(n)
	ut.show(wl)
	print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
	wl, n = ut.get_waiting_times(threshold_node_idx = 5)
	print(n)
	ut.show(wl)
	print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
