#! /usr/bin/env python
import math
import random
from ete2 import Tree, TreeStyle, TextFace, SeqGroup
from subprocess import call
from scipy.optimize import fmin
import collections
from collections import deque
from scipy import stats


t = Tree("test.tree.tre")
print(t)
t.set_outgroup( t&"A" )
print(t)
t.set_outgroup( t&"C" )
print(t)
#lvs = t.get_leaves()
#for lv in lvs:
#	lv.dist = lv.dist + 0.004



#t.write(outfile = "pk.tre", format=5)

