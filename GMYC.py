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

class speciation:
	def __init__(self, num_lineage, rate = 0 , p = 1.0):
		self.num_lineages = num_lineage #n
		self.rate = rate #speciation rate 
		self.p = p
		
	def __str__(self):
		s = "spe_event with n=" + repr(self.num_lineages) +  ", rate= " + repr(self.rate) + "\n" 
		return s
	
	def update(self,spe_rate, spe_p):
		self.rate = spe_rate
		self.p = spe_p
	
	def update_p(self, spe_p):
		self.p = spe_p
	
	def update_rate(self, spe_rate):
		self.rate = spe_rate
	
	def getBirthRate(self): # this is b
		if self.num_lineages == 0:
			return 0
		else:
			return self.rate * math.pow(self.num_lineages, self.p)
		
	def getNumDivEvent(self):
		return self.num_lineages-1


class coalescent:
	def __init__(self, num_individual, rate = 0, p = 1.0):
		 self.num_individual = num_individual #n
		 self.rate = rate # 1/2N
		 self.p = p
	
	def getCoalesecntRate(self):
		return self.rate * math.pow(self.num_individual * (self.num_individual - 1), self.p)
		
	def getNumDivEvent(self):
		return self.num_individual - 1 


class coalescents:
	def __init__(self, num_coalescent, rate = 0, p = 1, const_rate = True, const_p = True):
		self.coa_list = []
		self.num_coa =  num_coalescent
		self.rate = rate 
		self.p = p
		self.const_rate = const_rate
		self.const_p = const_p
	
	def __str__(self):
		s = ""
		cnt = 1
		for coa in self.coa_list:
			s = s + "coa_event" + repr(cnt) + ", with n= "+repr(coa.num_individual) + ", rate="+ repr(coa.rate)+ "\n" 
			cnt = cnt + 1
		return s
	
	def add_coalescent(self,coa):
		if self.const_rate:
			coa.rate = self.rate 
		if self.const_p:
			coa.p = self.p
		self.coa_list.append(coa)
		
	def getSumCoaRate(self): #this is b
		sr = 0
		for coa in self.coa_list:
			sr = sr + coa.getCoalesecntRate()
		return sr
	
	def update(self, rate = 0, p = 1):
		for coa in self.coa_list:
			if self.const_rate:
				coa.rate = rate
			if self.const_p:
				coa.p = p
	
	def update_p(self, p):
		for coa in self.coa_list:
			if self.const_p:
				coa.p = p
	
	def update_rate(self, rate):
		for coa in self.coa_list:
			if self.const_rate:
				coa.rate = rate	
				
	def getNumDivEvent(self):
		ndiv = 0
		for coa in self.coa_list:
			ndiv = ndiv + coa.getNumDivEvent()
		return ndiv


class waiting_time:
	def __init__(self, length, num_coas = 0, num_lines = 0):
		self.length = length #x
		self.spe = speciation(num_lineage = num_lines)
		self.coas  =  coalescents(num_coalescent = num_coas)
		self.b = self.spe.getBirthRate() + self.coas.getSumCoaRate()
		self.spe_rate = 0
		self.spe_p = 1
		self.spe_n = num_lines
		self.coa_rate = 0
		self.coa_p = 1
		self.coa_n = num_coas
	
	def __str__(self):
		s = "Waitting time = " + repr(self.length) + "\n"
		s = s + str(self.spe) + "\n"
		s = s + str(self.coas) + "\n"
		s = s + "---------------------------------------------------------\n"
		return s 
	
	def update(self, sp_rate, sp_p, coa_rate, coa_p):
		self.spe_rate = sp_rate 
		self.spe_p = sp_p
		self.coa_rate = coa_rate 
		self.coa_p = coa_p
		self.spe.update(sp_rate, sp_p)
		self.coas.update(coa_rate, coa_p)
		self.b = self.spe.getBirthRate() + self.coas.getSumCoaRate() 
	
	def update_p(self, sp_p, coa_p):
		self.spe_p = sp_p
		self.coa_p = coa_p
		self.spe.update_p(sp_p)
		self.coas.update_p(coa_p)
	
	def update_rate(self, sp_rate, coa_rate):
		self.spe_rate = sp_rate 
		self.coa_rate = coa_rate 
		self.spe.update_rate(sp_rate)
		self.coas.update_rate(coa_rate)
	
	def logl(self):
		self.b = self.spe.getBirthRate() + self.coas.getSumCoaRate()
		prob = self.b * math.exp (-1.0 * self.b * self.length)
		if prob <=0:
			return None
		else:
			return math.log(prob)
	
	def scaleSpeBranchL(self):
		return math.pow(self.spe_n, self.spe_p) * self.length
	
	def scaleCoaBranchL(self):
		br = 0
		for coa in self.coas.coa_list:
			#print(coa.num_individual)
			br = br + math.pow (coa.num_individual * (coa.num_individual-1), self.coa_p) * self.length 
		#print("branch length in wt: " + repr(br))
		return br


class tree_time:
	def __init__(self, wtimes, num_spe = 1):
		self.w_time_list = wtimes
		self.llh = 0
		self.spe_rate = 0.001
		self.spe_p = 1
		self.coa_rate = 0.001
		self.coa_p = 1 
		self.num_species = num_spe
		self.numSpeEvent = self.num_species - 1
		last_wc = self.w_time_list[-1]
		self.numCoaEvent = 0
		for coa in last_wc.coas.coa_list:
			self.numCoaEvent = self.numCoaEvent + coa.getNumDivEvent()
		
	def show(self):
		print("This is tree_time with spe event: " + repr(self.numSpeEvent) + ", coa event: " + repr(self.numCoaEvent)) 
	
	def sum_llh(self):
		self.llh = 0.0
		for w_time in self.w_time_list:
			if w_time.logl() == None:
				self.llh = -sys.float_info.max
				break
			else:
				self.llh = self.llh + w_time.logl()
		return self.llh
	
	def update(self, spe_p, coa_p):
		for w_time in self.w_time_list:
			w_time.update_p(spe_p, coa_p)
		spe_rate_dn = 0
		coa_rate_dn = 0 
		for w_time in self.w_time_list:
			spe_rate_dn = spe_rate_dn + w_time.scaleSpeBranchL()
			coa_rate_dn = coa_rate_dn + w_time.scaleCoaBranchL()
		
		if spe_rate_dn == 0:
			self.spe_rate = 0
		else:
			self.spe_rate = self.numSpeEvent/spe_rate_dn
			
		if coa_rate_dn == 0:
			self.coa_rate =0
		else:
			self.coa_rate = self.numCoaEvent/coa_rate_dn
		
		#print("spe rate dn:" + repr(spe_rate_dn))
		#print("cao_rate dn:" + repr(coa_rate_dn))
		
		for w_time in self.w_time_list:
			w_time.update_rate(self.spe_rate, self.coa_rate)


def tar_fun(x, *args):
	"""args[0] is tree_time"""
	spe_p = x[0]
	coa_p = x[1]
	args[0].update(spe_p, coa_p)
	return args[0].sum_llh() * (-1.0)


def optimize_all(tree):
	min_change = 1
	max_iters = 100
	best_llh = float("-inf")
	best_num_spe = -1
	utree = Ultrametric_tree.um_tree(tree)
	for tnode in utree.nodes:
		wt_list, num_spe = utree.get_waiting_times(threshold_node = tnode)
		tt = tree_time(wt_list)
		
		last_llh = float("-inf")
		change = float("inf")
		cnt = 0
		while change > min_change and cnt < max_iters:
			cnt = cnt + 1
			para = fmin(tar_fun, [1, 1], [tt])
			tt.update(para[0], para[1])
			change = abs(tt.sum_llh() - last_llh)
			
		if tt.sum_llh() > best_llh:
			best_llh = tt.sum_llh()
			best_num_spe = num_spe
		
	
	print("Highest llh:" + repr(best_llh))
	print("Num spe:" + repr(best_num_spe))

if __name__ == "__main__":
        #if len(sys.argv) != 6: 
        #    print("usage: ./ncbi_taxonomy.py <tree_of_life.tre> <id_name.txt> <id_rank.txt> <name_tax.txt> <outputfile>")
        #    sys.exit() 2mtree.tre
        #sf = species_finder("2mtree.tre")
        #sf = species_finder("test.tree.tre")
        
        #num_spe = sf.search()
        #print("Final No. Spe" )
        #print(num_spe)
        #utree = Ultrametric_tree.um_tree("test.tree.tre")
        #wt_list, num_spe = utree.get_waiting_times(threshold_node_idx = 0)
        #his = tree_time(wt_list)
        #his.optimize_all()
        #val = fmin(tar_fun, [1,1], [his,1])
        #print(val)
        optimize_all("test.tree.tre")
        #optimize_all("2mtree.tre")


