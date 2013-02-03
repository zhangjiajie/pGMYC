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
#import scipy
from scipy.optimize import fmin

class exp_distribution:
    def __init__(self):
		pass 
		
    def __str__(self):
        return "Nothing at all"
    
	def coalescent_prob(self, r, n , x):
		"""
		r: coalesce rate = 1/2N
		n: number of lineages 
		x: x
		"""
		prob = r * n * (n - 1.0) * math.exp( -1.0 * r * n * (n - 1.0) * x)
		return math.log(prob)
		
	def yule_prob(self,r, n , x):
		"""
		r: speciation rate
		n: number of lineages 
		x: x
		"""
		prob = r * n * math.exp(-1.0 * r * n * x)
		return math.log(prob)
		
	def gmyc_prob(self, b, x):
		prob = b * math.exp (-1.0 * b * x)
		return math.log(prob)
		
	
	def sum_log_l(self, bl, xl):
		s = 0.0
		for i in range(len(bl)):
			s = s + gmyc_prob(bl[i], xl[i])
		return s 


class speciation:
	def __init__(self, num_lineage, rate = 0 , p = 1.0):
		#self.length = 0 #x
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
			br = br + math.pow (coa.num_individual * (coa.num_individual-1), coa.p) * self.length 
		print("br:" + repr(br))
		return br
			
	def spe_rate_deriv1(self):
		deriv1 = self.b * math.exp (-1.0 * self.b * self.length) * (-1.0) * math.pow(self.spe_n, self.spe_p)
		return deriv1
	
	def spe_rate_deriv2(self):
		deriv2 = self.b * math.exp (-1.0 * self.b * self.length) * math.pow(self.spe_n, self.spe_p * 2.0)
		return deriv2 


class optimization:
	def __init__(self, range_a, range_b, step = 1):
		self.range_a = range_a
		self.range_b = range_b
		self.max_val = float("-inf") 
		self.max_x = range_a
		self.step = step
		self.curr_x = range_a
		self.last_val = float("-inf") 
		self.cont_incease = True
		
	
	def next_search(self, curr_val):
		#self.curr_x = self.curr_x + step
		if curr_val > self.max_val:
			self.max_val = curr_val
			self.max_x = self.curr_x
		
		self.curr_x = self.curr_x + self.step
		if self.curr_x > self.range_b:
			return False, self.curr_x
		else:
			return True, self.curr_x
			
	def next_newton(self, curr_val, d1_d2):
		pass


def tar_fun(x, *args):
	spe_p = x[0]
	coa_p = x[1]
	args[0].update(spe_p, coa_p)
	return -args[0].sum_llh()


class tree_time:
	def __init__(self, wtimes, num_spe = 1, step = 1, maxiters = 100):
		self.w_time_list = wtimes
		self.llh = 0
		#self.spe_rate = random.random()
		self.spe_rate = 0.001
		self.spe_p = 1
		#self.coa_rate = random.random()
		self.coa_rate = 0.001
		self.coa_p = 1 
		self.step = step
		self.maxiters = maxiters
		self.num_species = num_spe
		self.numSpeEvent = self.num_species - 1
		last_wc = self.w_time_list[-1]
		self.numCoaEvent = 0
		for coa in last_wc.coas.coa_list:
			self.numCoaEvent = self.numCoaEvent + coa.getNumDivEvent()
		
	def show(self):
		print("This is tree_time") 
	
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
		
		print(spe_rate_dn)
		print(coa_rate_dn)
		
		for w_time in self.w_time_list:
			w_time.update_rate(self.spe_rate, self.coa_rate)


class species_finder:
	def __init__(self, tree_file):
		self.utree = Ultrametric_tree.um_tree(tree_file)
		self.bestll = float("-inf")
		self.threshold = None
		self.num_spe = 2
	
	def search(self):
		curr_spe = 2
		for tnode in self.utree.nodes:
			print("Next node ----------------------------------------------")
			wtimes, numspe= self.utree.get_waiting_times(threshold_node = tnode)
			tt = tree_time(wtimes, step = 0.001)
			llh = tt.optimize_all(min_change = 0.01)

			#print(tt.coa_p)
			#print(tt.coa_p)
			#print(tt.spe_rate)
			if llh > self.bestll:
				self.bestll = llh
				self.threshold = tnode
				self.num_spe = numspe
				#curr_spe = curr_spe + 1
			print(llh)
			#print(numspe)
			print(self.num_spe)
			print(tt.spe_rate)
			print(tt.coa_rate)
		return self.num_spe


def optimize_all(tree):
	best_llh = float("-inf")
	best_num_spe = -1
	utree = Ultrametric_tree.um_tree(tree)
	for tnode in utree.nodes:
		wt_list, num_spe = utree.get_waiting_times(threshold_node = tnode)
		tt = tree_time(wt_list)
		para = fmin(tar_fun, [1, 1], [tt])
		tt.update(para[0], para[1])
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


