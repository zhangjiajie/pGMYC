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
	
	def getBirthRate(self): # this is b
		return self.rate * math.pow (self.num_lineages, self.p)


class coalescent:
	def __init__(self, num_individual, rate = 0, p = 1.0):
		 self.num_individual = num_individual #n
		 self.rate = rate # 1/2N
		 self.p = p
	
	def getCoalesecntRate(self):
		return self.rate * math.pow(self.num_individual * (self.num_individual - 1), self.p)


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


class waiting_time:
	def __init__(self, length, num_coas = 0, num_lines = 0):
		self.length = length #x
		self.spe = speciation(num_lineage = num_lines)
		self.coas  =  coalescents(num_coalescent = num_coas)
		self.b = self.spe.getBirthRate() + self.coas.getSumCoaRate()
		self.spe_rate = 0
		self.spe_p = 0
		self.spe_n = num_lines
		self.coa_rate = 0
		self.coa_p = 0
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
	
	def logl(self):
		self.b = self.spe.getBirthRate() + self.coas.getSumCoaRate()
		prob = self.b * math.exp (-1.0 * self.b * self.length)
		#print("b")
		#print(self.b)
		#print("L")
		#print(self.length)
		#print("prob")
		#print(prob)
		if prob <=0:
			return None
		else:
			return math.log(prob)
			
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
	spe_rate = x[0]
	spe_p = x[1]
	coa_rate = x[2]
	coa_p = x[3]
	args[0].update(spe_rate, spe_p, coa_rate, coa_p)
	return args[0].sum_llh()

class tree_time:
	def __init__(self, wtimes, step = 1, maxiters = 100):
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
	
	def update(self, spe_rate, spe_p, coa_rate, coa_p):
		for w_time in self.w_time_list:
			w_time.update(spe_rate, spe_p, coa_rate, coa_p)
	
	def optimize_spe_rate(self):
		optim = optimization(range_a = 0, range_b = 1, step = self.step)
		opt_flag = True
		next_sp_rate = 0
		while opt_flag:
			self.update(next_sp_rate, self.spe_p, self.coa_rate, self.coa_p)
			val = self.sum_llh()
			opt_flag, next_sp_rate = optim.next_search(val)
		self.spe_rate = optim.max_x
		self.update(self.spe_rate, self.spe_p, self.coa_rate, self.coa_p)
		return optim.max_val
	
	def optimize_spe_p(self):
		optim = optimization(range_a = 0.01, range_b = 10, step = self.step)
		opt_flag = True
		next_sp_p = 0
		while opt_flag:
			self.update(self.spe_rate, next_sp_p, self.coa_rate, self.coa_p)
			val = self.sum_llh()
			opt_flag, next_sp_p = optim.next_search(val)
		self.spe_p = optim.max_x
		self.update(self.spe_rate, self.spe_p, self.coa_rate, self.coa_p)
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
		self.update(self.spe_rate, self.spe_p, self.coa_rate, self.coa_p)
		return optim.max_val
	
	def optimize_coa_p(self):
		optim = optimization(range_a = 0.01, range_b = 10, step = self.step)
		opt_flag = True
		next_coa_p = 0
		while opt_flag:
			self.update(self.spe_rate, self.spe_p, self.coa_rate, next_coa_p)
			val = self.sum_llh()
			opt_flag, next_coa_p = optim.next_search(val)
		self.coa_p = optim.max_x
		self.update(self.spe_rate, self.spe_p, self.coa_rate, self.coa_p)
		return optim.max_val
	
	def optimize_all(self, min_change = 1):
		change = 100000000.0
		last_val = -999999999999999
		cnt = 0
		while cnt < self.maxiters and change > min_change:
			self.optimize_spe_rate()
			self.optimize_coa_rate()
			self.optimize_spe_p()
			self.optimize_coa_p()
			#print("highest llh:" + repr(self.sum_llh()) + ", with spe_rate = " + repr(self.spe_rate) + " coa_rate = " + repr(self.coa_rate))
			change = abs(self.sum_llh() - last_val)
			last_val = self.sum_llh()
			cnt = cnt + 1
		self.llh = last_val
		return self.llh


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


if __name__ == "__main__":
        #if len(sys.argv) != 6: 
        #    print("usage: ./ncbi_taxonomy.py <tree_of_life.tre> <id_name.txt> <id_rank.txt> <name_tax.txt> <outputfile>")
        #    sys.exit() 2mtree.tre
        #sf = species_finder("2mtree.tre")
        #sf = species_finder("test.tree.tre")
        
        #num_spe = sf.search()
        #print("Final No. Spe" )
        #print(num_spe)
        utree = Ultrametric_tree.um_tree("test.tree.tre")
        wt_list, num_spe = utree.get_waiting_times(threshold_node_idx = 0)
        his = tree_time(wt_list)
        #his.optimize_all()
        val = fmin(tar_fun, [0.01, 0.01, 0.01, 0.01], [his,1])
        print(val)


