#! /usr/bin/env python
import sys
import os
import json
import operator
from ete2 import Tree, TreeStyle, TextFace, SeqGroup

class um_tree:
	def __init__(self, tree):
		self.tree = Tree(tree)
	
	
	
if __name__ == "__main__":
	print("main function")
	t = Tree("2mtree.tre")
	t.show()
