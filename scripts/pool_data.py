'''
pool_data.py

Pools the data collected in the preliminary and validation screens. Saves data as a pkl file

Written by Nicolas Quach
'''
#Import packages
import matplotlib
matplotlib.use('Agg')
import numpy as np 
import os
import tifffile as tiff
from skimage.io import imread
from skimage.measure import label, regionprops
import scipy
import matplotlib.pyplot as plt
import cPickle as pickle
from SLIP_functions import analyze_well, analyze_plate, segment_SLIP, plot_slip_well, plot_slip_wells, plot_slip_wells_gmm, plot_slip_wells_lysis_posterior, plot_slip_wells_MOI_posterior, save_slip_wells
from SLIP_functions import plot_slip_joint_plot, fit_kde, compute_p_values, save_classified_wells
from SLIP_functions import classify_infections_gmm, classify_infections_gmm2
from keio_names import get_keio_names, pos_to_strain
import seaborn as sns
import pandas as pd
import pymc3 as pm
import json

#this function is for quality control, removes lists of lists
def flatten_list(l):
	flattened = []
	for item in l:
		if type(item) is list:
			for i in item:
				flattened.append(i)
		else:
			flattened.append(item)
	return flattened

#Set root directory
direc = "/scratch/users/nquach/datatxt/"
#open pickled files for pass 1 and pass 2 of the keio screen
p1_classified = open(os.path.join(direc, 'pass1_all_data.pkl'), 'r')
p2_classified = open(os.path.join(direc, 'pass2_all_data.pkl'), 'r')

p1_dict = pickle.load(p1_classified)
p2_dict = pickle.load(p2_classified)

pooled_dict = {}
#load the first pass hit gene names
titles = p2_dict.keys()

#find the genes that were not first pass hits
p1_titles = np.setdiff1d(p1_dict.keys(),titles)

for title in titles:
	print title, (title in p1_dict.keys()), (title in p2_dict.keys())
	#Uncomment if pooling classified dataset
	#pooled_dict[title] = [flatten_list(p1_dict[title][0] + p2_dict[title][0]), flatten_list(p1_dict[title][1] + p2_dict[title][1]), flatten_list(p1_dict[title][2] + p2_dict[title][2])]
	#Uncomment if pooling unclassified dataset
	pooled_dict[title] = [flatten_list(list(p1_dict[title][0]) + list(p2_dict[title][0])), flatten_list(list(p1_dict[title][1]) + list(p2_dict[title][1]))]
	if title == 'pdhR':
		print list(p2_dict[title][0])

for title in p1_titles:
	print title
	pooled_dict[title] = p1_dict[title]

pickle.dump(pooled_dict, open(os.path.join(direc, 'all_pooled_data.pkl'),'w+'), pickle.HIGHEST_PROTOCOL)

