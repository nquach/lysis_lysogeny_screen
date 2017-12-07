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

#set root directory
direc = "/scratch/users/nquach/datatxt/"
#list numbers of the plates you wish to incorporate into the compiled .pkl file
plate_numbers =['1_1','1_2', 3, 5, 7, '9_1','9_2', 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47, 49, 51, 53, 55, 57, 59, 61, 63, 65, 67, 69, 71, 73, 75, 77, 79, 81, 83, 85, 87, 89, 91, 93, 95]

all_data_dict = {}
classified_data_dict = {}

#iterate through all plates
for plate_number in plate_numbers:
	print plate_number
	#create a file for unclassified data
	all_file_name = 'all_data_K' + str(plate_number) + '.pkl'
	classified_file_name = 'classified_data_K' + str(plate_number) + '.pkl'
	all_file_path = os.path.join(direc, all_file_name)
	classified_file_path = os.path.join(direc, classified_file_name)
	all_file = open(all_file_path, 'r')
	#create a file for classified data
	classified_file = open(classified_file_path, 'r')
	all_data = pickle.load(all_file)
	classified_data = pickle.load(classified_file)
	all_titles = all_data.keys()
	classified_titles = classified_data.keys()
	for title in all_titles:
		#if a gene name appears twice, pool the unclassified data
		if title in all_data_dict.keys():
			all_data_dict[title][0] = all_data_dict[title][0] + list(all_data[title][0])  
			all_data_dict[title][1] = all_data_dict[title][1] + list(all_data[title][1])
		else:
			all_data_dict[title] = all_data[title]

	for title in classified_titles:
		#if a gene name appears twice, pool the classified data
		if title in classified_data_dict.keys():
			classified_data_dict[title][0] = classified_data_dict[title][0] + classified_data[title][0]
			classified_data_dict[title][1] = classified_data_dict[title][1] + classified_data[title][1]
			classified_data_dict[title][2] = classified_data_dict[title][2] + classified_data[title][2]
		else:
			classified_data_dict[title] = classified_data[title]

all_data_name = 'pass1_all_data.pkl'
classified_name = 'pass1_classified_data.pkl'

pickle.dump(all_data_dict, open(os.path.join(direc, all_data_name), 'w+'), pickle.HIGHEST_PROTOCOL)
pickle.dump(classified_data_dict, open(os.path.join(direc, classified_name), 'w+'), pickle.HIGHEST_PROTOCOL)
