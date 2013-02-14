#! /usr/bin/env python
import Ultrametric_tree
import math
import random
from ete2 import Tree, TreeStyle, TextFace, SeqGroup
from subprocess import call
from scipy.optimize import fmin
import collections
from collections import deque
from scipy import stats


t = Tree("Randomly_resolved.tre")
lvs = t.get_leaves()
for lv in lvs:
	lv.dist = lv.dist + 0.004



t.write(outfile = "pk.tre", format=5)

