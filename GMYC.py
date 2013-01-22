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
	def __init__(self, num_lineage, rate = 0 , p = 1):
		self.length = 0 #x
		self.num_lineages = num_lineage #n
		self.rate = rate #speciation rate 
		self.p = p
		
	def update(self,spe_rate, spe_p):
		self.rate = spe_rate
		self.p = spe_p
	
	def getBirthRate(self): # this is b
		return self.rate * self.num_lineages
		

class coalescent:
	def __init__(self, num_individual, rate, p):
		 self.num_individual = num_individual #n
		 self.rate = rate # 1/2N
		 self.p = p
	
	def getCoalesecntRate(self):
		return self.rate * self.num_individual * (self.num_individual - 1)


class coalescents:
	def __init__(self, num_coalescent, rate = -1, p = -1, const_rate = True, const_p = True):
		self.coa_list = []
		self.num_coa =  num_coalescent
		self.rate = rate 
		self.p = p
		self.const_rate = const_rate
		self.const_p = const_p
		
	def add_coalescent(self,coa):
		if self.const_rate:
			coa.rate = self.rate 
		if self.const_p:
			coa.p = self.p	
		self.coa_list.append(coa)
		
	def getSumCoaRate(self): #this is b
		if self.num_coa == len(self.coa_list):
			sr = 0
			for coa in self.coa_list:
				sr = sr + coa.getCoalesecntRate()
			return sr
		else:
			print("Coalescents events are not proper init")
			return None
	
	def update(self, rate = 0, p = 1):
		for coa in self.coa_list:
			if self.const_rate:
				coa.rate = rate
			if self.const_p:
				coa.p = p
	
	
class waiting_time:
	def __init__(self, length, num_coas = 0, num_lines = 0):
		self.length = length #x
		self.spe = speciation(num_lineage = num_lines)
		self.coas  =  coalescents(num_coalescent = num_coas)
		self.b = 0.0
	
	def update(self, sp_rate, sp_p, coa_rate, coa_p):
		self.spe.update(sp_rate, sp_p)
		self.coas.update(coa_rate, coa_p)
		self.b = self.spe.getBirthRate() + self.coas.getSumCoaRate() 
	
	def logl(self):
		self.b = self.spe.getBirthRate() + self.coas.getSumCoaRate() 
		prob = self.b * math.exp (-1.0 * self.b * self.length)
		return math.log(prob)
			
			
class optimization:
	def __init__(self, range_a, range_b, step = 0.001):
		self.range_a = range_a
		self.range_b = range_b
		self.max_val = 0 
		self.max_x = range_a
		self.step = step
		self.curr_x = range_a
	
	def next_search(self, curr_val):
		#self.curr_x = self.curr_x + step
		if curr_val > self.max_val:
			self.max_val = curr_val
			self.max_x = self.curr_x
		
		self.curr_x = self.curr_x + step
		if self.curr_x > self.range_b:
			return False, self.curr_x
		else:
			return True, self.curr_x
			
		
class tree_time:
	
	def __init__(self, step = 0.01, maxiters = 100):
		self.w_time_list = []
		self.llh = 0
		self.spe_rate = 0.1
		self.spe_p = 1
		self.coa_rate = 0.1
		self.coa_p = 1 
		self.step = step
		self.maxiters = maxiters
	
	def init_w_times(self, tree, threshold):
		"""This should init the waiting times from a tree"""
		pass
	
	def sum_llh(self):
		for w_time in w_time_list:
			self.llh = self.llh + w_time.logl()
		return self.llh
	
	def update(self, spe_rate, spe_p, coa_rate, coa_p):
		for w_time in self.w_time_list:
			w_time.update(spe_rate, spe_p, coa_rate, coa_p)

	
	def optimize_spe_rate(self):
		optim = optimization(range_a = 0, range_b = 2, step = self.step)
		opt_flag = True
		next_sp_rate = 0
		while opt_flag:
			self.update(next_sp_rate, self.spe_p, self.coa_rate, self.coa_p)
			val = self.sum_llh()
			opt_flag, next_sp_rate = optim.next_search(val)
		self.spe_rate = optim.max_x
		return optim.max_val
	
	def optimize_spe_p(self):
		optim = optimization(range_a = 0, range_b = 2, step = self.step)
		opt_flag = True
		next_sp_p = 0
		while opt_flag:
			self.update(self.spe_rate, next_sp_p, self.coa_rate, self.coa_p)
			val = self.sum_llh()
			opt_flag, next_sp_p = optim.next_search(val)
		self.spe_p = optim.max_x
		return optim.max_val
		
	def optimize_coa_rate(self):
		optim = optimization(range_a = 0, range_b = 1, step = self.step)
		opt_flag = True
		next_coa_rate = 0
		while opt_flag:
			self.update(self.spe_rate, self.spe_p, next_coa_rate, self.coa_p)
			val = self.sum_llh()
			opt_flag, next_coa_rate = optim.next_search(val)
		self.coa_rate = optim.max_x
		return optim.max_val
	
	def optimize_coa_p(self):
		optim = optimization(range_a = 0, range_b = 2, step = self.step)
		opt_flag = True
		next_coa_p = 0
		while opt_flag:
			self.update(self.spe_rate, self.spe_p, self.coa_rate, next_coa_p)
			val = self.sum_llh()
			opt_flag, next_coa_p = optim.next_search(val)
		self.coa_p = optim.max_x
		return optim.max_val
	
	def optimize_all(self, min_change = 0.1):
		change = 1.0
		last_val = 0 
		cnt = 0
		while change > min_change and cnt < self.maxiters:
			self.optimize_spe_rate()
			self.optimize_coa_rate()
			self.optimize_spe_p()
			curr_val = self.optimize_coa_p()
			change = curr_val - last_val
			last_val = curr_val
			cnt = cnt + 1
		
		self.llh = last_val
		return last_val
			
		
					
		
		
		
	
		
		

    




if __name__ == "__main__":
        if len(sys.argv) != 6: 
            print("usage: ./ncbi_taxonomy.py <tree_of_life.tre> <id_name.txt> <id_rank.txt> <name_tax.txt> <outputfile>")
            sys.exit()
        t = ncbi_taxa()
        t.init_tax_tree(sys.argv[1], sys.argv[2], sys.argv[3])
        t.extract_sub_tax_tree(sys.argv[4], sys.argv[5])