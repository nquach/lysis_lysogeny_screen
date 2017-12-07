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

#Set root directory
direc = "/scratch/users/nquach/datatxt/"
#Set directory where you want to save the .txt JSON file
save_direc = "/scratch/users/nquach/plots/"
#Load dataset
pooled_file = open(os.path.join(direc, 'pooled_classified_data.pkl'), 'r')

pooled_dict = pickle.load(pooled_file)

titles = list(pooled_dict.keys())

index = 0

stats_dict = {}

for i in range(len(titles)):
	title = titles[i]
	print title
	#Calculate cell fate statistics
	lytic_list = pooled_dict[title][0]
	lyso_list = pooled_dict[title][1]
	uninfected_list = pooled_dict[title][2]
	num_lytic = len(lytic_list)
	num_lyso = len(lyso_list)
	num_uninfected = len(uninfected_list)
	total_cells = num_lytic + num_lyso + num_uninfected
	num_infected = num_lytic + num_lyso	
	a = float(1+num_lytic)
	b = float(1+num_lyso)
	mode = 0.5
	if a+b > 2:
		mode = (a-1)/(a+b-2)
		var = a*b/(np.square(a+b)*(a+b+1))
	stats_dict[title] = [total_cells, num_infected, num_uninfected, num_lytic, num_lyso, mode, var]

#Save dict containing cell fate statistics, indexed by gene name
stats_dict_name = 'pooled_hits_stats.txt'
stats_dict_path = os.path.join(direc, stats_dict_name)
stats_dict_file = open(stats_dict_path, 'w+')
json.dump(stats_dict, stats_dict_file)

stats_dict_file.close()
pooled_file.close()
