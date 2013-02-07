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
from scipy.optimize import fmin_powell
from scipy.optimize import fmin_l_bfgs_b
from scipy.optimize import fmin_tnc
from scipy import stats

class lh_ratio_test:
	def __init__(self, null_llh, llh, df):
		self.lr = 2.0 * (llh - null_llh)
		self.p = 1 - stats.chi2.cdf(self.lr, df)
	
	def get_p_value(self):
		return self.p


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
		
	def bprime(self):
		if self.num_lineages == 0:
			return 0.0
		else:
			bp = self.p * self.rate * math.pow(self.num_lineages, (self.p - 1.0))
			return bp


class coalescent:
	def __init__(self, num_individual, rate = 0, p = 1.0):
		 self.num_individual = num_individual #n
		 self.rate = rate # 1/2N
		 self.p = p
	
	def getCoalesecntRate(self):
		return self.rate * math.pow(self.num_individual * (self.num_individual - 1.0), self.p)
		
	def getNumDivEvent(self):
		return self.num_individual -1
	
	def bprime(self):
		bp = self.p * self.rate * math.pow(self.num_individual * (self.num_individual - 1.0), (self.p - 1.0))
		return bp


class coalescents:
	def __init__(self, num_coalescent, rate = 0, p = 1, const_rate = True, const_p = True):
		self.coa_list = []
		self.num_coa =  num_coalescent
		self.rate = rate 
		self.p = p
		self.const_rate = const_rate
		self.const_p = const_p
		self.coas_idx = -1
	
	def __str__(self):
		s = ""
		cnt = 1
		for coa in self.coa_list:
			s = s + "coa_event" + repr(cnt) + ", with n= "+repr(coa.num_individual) + ", rate="+ repr(coa.rate)+ "\n" 
			cnt = cnt + 1
		return s
	
	def add_coalescent(self,coa):
		coa.rate = self.rate 
		coa.p = self.p
		self.coa_list.append(coa)
		
	def getSumCoaRate(self): #this is b
		sr = 0
		for coa in self.coa_list:
			sr = sr + coa.getCoalesecntRate()
		return sr
	
	def update(self, rate = 0, p = 1):
		self.p = p
		self.rate = rate
		for coa in self.coa_list:
			coa.rate = rate
			coa.p = p
	
	def update_p(self, p):
		self.p = p
		for coa in self.coa_list:
				coa.p = p
				
	
	def update_rate(self, rate):
		self.rate = rate
		for coa in self.coa_list:
				coa.rate = rate
				
				
	def getNumDivEvent(self):
		ndiv = 0
		for coa in self.coa_list:
			ndiv = ndiv + coa.getNumDivEvent()
		return ndiv
		
	def bprime(self):
		bp = 0.0
		for coa in self.coa_list:
			bp = bp + coa.bprime()
		return bp


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
		self.num_lines = None
		self.num_curr_coa = 0
	
	def count_num_lines(self):
		"""This is used for the null model only!!!"""
		self.num_lines = self.spe_n
		for coa in self.coas.coa_list:
			self.num_lines = self.num_lines + coa.num_individual
		self.num_lines = self.num_lines * (self.num_lines - 1)
		
	
	def __str__(self):
		s = "Waitting time = " + repr(self.length) + "\n"
		s = s + str(self.spe) + "\n"
		s = s + str(self.coas)# + "\n"
		s = s + "Numer curr coa lines:"+ str(self.num_curr_coa) + "\n"
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
			print("branch lenghth:" + repr(self.length))
			print("b:" + repr(self.b))
			return None
		else:
			return math.log(prob)
	
	def scaleSpeBranchL(self):
		return math.pow(self.spe_n, self.spe_p) * self.length
	
	def scaleCoaBranchL(self):
		br = 0
		for coa in self.coas.coa_list:
			br = br + math.pow(coa.num_individual * (coa.num_individual-1), self.coa_p) * self.length 
		return br
		
	def bprime_spe(self):
		bp = self.spe.bprime()
		bps = bp * math.exp(-1.0 * self.b * self.length) + self.b * math.exp (-1.0 * self.b * self.length) * -1.0 * self.length * bp
		return bps
		
	def bprime_coa(self):
		bp = self.coas.bprime()
		bps = bp * math.exp(-1.0 * self.b * self.length) + self.b * math.exp (-1.0 * self.b * self.length) * -1.0 * self.length * bp
		return bps


class tree_time:
	def __init__(self, wtimes, num_spe):
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
			self.coa_rate = 0
		else:
			self.coa_rate = self.numCoaEvent/coa_rate_dn
		
	def show(self):
		print("This is tree_time with spe event: " + repr(self.numSpeEvent) + ", coa event: " + repr(self.numCoaEvent)) 
	
	def sum_llh(self):
		self.llh = 0.0
		for w_time in self.w_time_list:
			if w_time.logl() == None:
				self.llh = -sys.float_info.max
				print("wtime logl infinity!!")
				break
			else:
				self.llh = self.llh + w_time.logl()
		return self.llh
	
	def update(self, spe_p, coa_p):
		self.spe_p = spe_p
		self.coa_p = coa_p
		
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
		
		
		for w_time in self.w_time_list:
			w_time.update_rate(self.spe_rate, self.coa_rate)
	
	def bprime_spe(self):
		bp = 0.0
		for wt in self.w_time_list:
			bp = bp + wt.bprime_spe()
		return bp
		
	def bprime_coa(self):
		bp = 0.0
		for wt in self.w_time_list:
			bp = bp + wt.bprime_coa()
		return bp


#The null model using Coalescent
class null_model:
	def __init__(self, wt_list, tree):
		nodes = tree.get_leaves()
		self.num_speEvent = len(nodes) - 1
		self.wt_list = wt_list
		self.p = 1.0
		self.rate = 0.0
	
	def logl(self, p = 1.0):
		self.p = p
		br_de = 0.0
		for wt in self.wt_list:
			br_de = br_de + wt.length * math.pow(wt.num_lines, self.p)
		self.rate = self.num_speEvent / br_de
		logl = 0
		for wt in self.wt_list:
			logl = logl + math.log(self.rate * math.pow(wt.num_lines, self.p) * math.exp(self.rate * math.pow(wt.num_lines, self.p) * wt.length * -1.0))
		return logl


def tar_fun(x, *args):
	"""args[0] is tree_time"""
	spe_p = x[0]
	coa_p = x[1]
	args[0].update(spe_p, coa_p)
	return args[0].sum_llh() * (-1.0)


def tar_fun_null(x, *args):
	return args[0].logl(p = x[0]) * (-1.0)


def prime_fun(x, *args):
	spe_p = x[0]
	coa_p = x[1]
	args[0].update(spe_p, coa_p)
	return [args[0].bprime_spe() * (-1.0) , args[0].bprime_coa() * (-1.0)]


def optimize_null_model(umtree):
	min_change = 0.1
	max_iters = 100
	wt_list, num_spe = umtree.get_waiting_times(threshold_node_idx = 0)
	nm = null_model(wt_list, umtree.tree)
	last_llh = float("-inf")
	change = float("inf")
	cnt = 0
	while change > min_change and cnt < max_iters:
		cnt = cnt + 1
		para, nn, cc = fmin_l_bfgs_b(tar_fun_null, [1], args = [nm], disp = False, bounds = [[0, 100]], approx_grad = True)
		curr_logl = nm.logl(p = para[0])
		change = abs(curr_logl - last_llh)
		last_llh = curr_logl
	return last_llh


def optimize_all(tree, print_detail = False):
	min_change = 0.1
	max_iters = 100
	best_llh = float("-inf")
	best_num_spe = -1
	utree = Ultrametric_tree.um_tree(tree)
	for tnode in utree.nodes:
		wt_list, num_spe = utree.get_waiting_times(threshold_node = tnode)
		tt = tree_time(wt_list, num_spe)
		last_llh = float("-inf")
		change = float("inf")
		cnt = 0
		
		while change > min_change and cnt < max_iters:
			cnt = cnt + 1
			para, nn, cc = fmin_l_bfgs_b(tar_fun, [1, 1], args = [tt], disp = False, bounds = [[0, 100], [0, 100]], approx_grad = True)
			#para, nn, cc = fmin_tnc(tar_fun, [0, 0], args = [tt], disp = False, bounds = [[0, 10], [0, 10]], approx_grad = True)
			#para, nn, cc = fmin_tnc(tar_fun, [1.0, 1.0], args = [tt], disp = False, bounds = [[0, 10], [0, 10]], fprime = prime_fun ) #, approx_grad = True)
			tt.update(para[0], para[1])
			
			change = abs(tt.sum_llh() - last_llh)
			last_llh = tt.sum_llh()
		
		if print_detail:
			print("Num spe:" + repr(num_spe) + ": " + repr(tt.sum_llh()))
			print("spe_lambda:" + repr(tt.spe_rate))
			print("coa_lambda:" + repr(tt.coa_rate))
			print("spe_p:" + repr(tt.spe_p))
			print("coa_p:" + repr(tt.coa_p))
			print("-----------------------------------------------------")
		if tt.sum_llh() > best_llh:
			best_llh = tt.sum_llh()
			best_num_spe = num_spe
			
	print("Highest llh:" + repr(best_llh))
	print("Num spe:" + repr(best_num_spe))
	
	null_logl = optimize_null_model (utree)
	print("Null llh:" + repr(null_logl))
	
	lrt = lh_ratio_test(null_llh = null_logl, llh = best_llh, df = 2)
	print("P-value:" + repr(lrt.get_p_value()))
	if lrt.get_p_value() >= 0.5:
		return 0
	else:
		return best_num_spe


if __name__ == "__main__":
	if len(sys.argv) != 3: 
		print("usage: ./GMYC.py <um_tree.tre>  <print p/n>")
		sys.exit()
	if sys.argv[1] == "p":
		optimize_all(sys.argv[1], print_detail = True)
	else:
		optimize_all(sys.argv[1])


