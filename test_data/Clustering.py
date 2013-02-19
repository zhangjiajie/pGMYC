#! /usr/bin/env python
from scipy.cluster.vq import vq, kmeans2, whiten
from numpy import array

class cluster:
	def __init__(self, feature, k , method = "kmean"):
		self.centroid = None 
		self.label = None 
		self.k = k
		self.feature = feature
	
	def find_cluster(self):
		self.centroid, self.label = kmeans2(self.feature, self.k)
		return self.label
