#! /usr/bin/env python
import sys
import math
import collections
import re
import json
import os 
import glob
import EXP
import subprocess
import random
from ete2 import Tree, TreeStyle, TextFace, SeqGroup, NodeStyle
from scipy.optimize import fmin
from collections import deque
from scipy import stats
from numpy import array
from subprocess import call


def trim_alignment(sfin, sfout, size = 100, seql = 1000):
	fin = open(sfin,"r")
	fout = open(sfout,"w")
	head=fin.readline()
	fout.write(head)
	l = fin.readline().strip()
	while l!="":
		ls=l.split()
		tname=ls[0]
		if tname.endswith("r"):
			fout.write(ls[0] + "	" + ls[1] + "\n")
		else:
			gaps = "-" * (seql - size)
			newseq = ls[1][:size] + gaps
			s = ls[0] + "	" + newseq + "\n"
			fout.write(s)
			
		l = fin.readline().strip()
	fin.close()
	fout.close()
	
def batch_trim(folder, size = 100):
	phyl = glob.glob(folder + "*.phy")
	os.makedirs(folder + "l"+repr(size))
	newfolder = folder + "l"+repr(size) + "/"
	cnt = 1
	for phy in phyl:
		trim_alignment(sfin = phy, sfout = newfolder + "t" + repr(cnt) + ".phy", size = size, seql = 1000)
		cnt = cnt + 1


if __name__ == "__main__":
	#trim_alignment(sfin = "simulated_set_1_1.phy", sfout = "test2/t.phy", size = 100, seql = 1000)
	
	batch_trim("/home/zhangje/Research/SpeciesCount/simulation/set_1/", size = 400)
	batch_trim("/home/zhangje/Research/SpeciesCount/simulation/set_1/", size = 500)
	batch_trim("/home/zhangje/Research/SpeciesCount/simulation/set_1/", size = 600)
	batch_trim("/home/zhangje/Research/SpeciesCount/simulation/set_1/", size = 700)
	batch_trim("/home/zhangje/Research/SpeciesCount/simulation/set_1/", size = 800)
	
	
	batch_trim("/home/zhangje/Research/SpeciesCount/simulation/set_3/", size = 400)
	batch_trim("/home/zhangje/Research/SpeciesCount/simulation/set_3/", size = 500)
	batch_trim("/home/zhangje/Research/SpeciesCount/simulation/set_3/", size = 600)
	batch_trim("/home/zhangje/Research/SpeciesCount/simulation/set_3/", size = 700)
	batch_trim("/home/zhangje/Research/SpeciesCount/simulation/set_3/", size = 800)
	
	
	batch_trim("/home/zhangje/Research/SpeciesCount/simulation/set_5/", size = 400)
	batch_trim("/home/zhangje/Research/SpeciesCount/simulation/set_5/", size = 500)
	batch_trim("/home/zhangje/Research/SpeciesCount/simulation/set_5/", size = 600)
	batch_trim("/home/zhangje/Research/SpeciesCount/simulation/set_5/", size = 700)
	batch_trim("/home/zhangje/Research/SpeciesCount/simulation/set_5/", size = 800)
	
	
	batch_trim("/home/zhangje/Research/SpeciesCount/simulation/set_05/", size = 400)
	batch_trim("/home/zhangje/Research/SpeciesCount/simulation/set_05/", size = 500)
	batch_trim("/home/zhangje/Research/SpeciesCount/simulation/set_05/", size = 600)
	batch_trim("/home/zhangje/Research/SpeciesCount/simulation/set_05/", size = 700)
	batch_trim("/home/zhangje/Research/SpeciesCount/simulation/set_05/", size = 800)
	
	batch_trim("/home/zhangje/Research/SpeciesCount/simulation/set_10/", size = 400)
	batch_trim("/home/zhangje/Research/SpeciesCount/simulation/set_10/", size = 500)
	batch_trim("/home/zhangje/Research/SpeciesCount/simulation/set_10/", size = 600)
	batch_trim("/home/zhangje/Research/SpeciesCount/simulation/set_10/", size = 700)
	batch_trim("/home/zhangje/Research/SpeciesCount/simulation/set_10/", size = 800)
	
	batch_trim("/home/zhangje/Research/SpeciesCount/simulation/set_100/", size = 400)
	batch_trim("/home/zhangje/Research/SpeciesCount/simulation/set_100/", size = 500)
	batch_trim("/home/zhangje/Research/SpeciesCount/simulation/set_100/", size = 600)
	batch_trim("/home/zhangje/Research/SpeciesCount/simulation/set_100/", size = 700)
	batch_trim("/home/zhangje/Research/SpeciesCount/simulation/set_100/", size = 800)
	#batch_trim("/home/zhangje/Research/SpeciesCount/simulation/set_1/", size = 350)
