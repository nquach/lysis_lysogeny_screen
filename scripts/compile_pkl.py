'''
compile_pkl.py
Compiles pkl files containing the raw and classified fluorescence values for cells in a plate.
Expect the script to hit errors on certain plates that require specific classification wells. Those that require K29H12 as the classification well, 
change classify_infections_gmm to classify_infections_gmm2. See the Appendix of the manuscript for the specific classification wells.
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
import utils

ROOT_DIREC = utils.ROOT_DIREC

#Define root directory path
direc = "/scratch/users/nquach/"
#list the numbers of the plates you wish to turn into .pkl files
plate_numbers = ['1_1','1_2', 3, 5, 7, '9_1','9_2', 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 39, 41, 45, 47, 49, 51, 53, 55, 57, 59, 61, 63, 65, 67, 69, 71, 73, 75, 77, 79, 81, 83, 85, 89, 91, 95]
#list classification wells to use
classification_wells = None

#Import gene name to position map for first pass hits
datatxt_direc = os.path.join(ROOT_DIREC, 'datatxt')
titles_dict_path = os.path.join(datatxt_direc, 'hits_pos_to_name.txt')
titles_file = open(titles_dict_path,'r')
titles_dict = json.load(titles_file)
hits_dict = titles_dict['4']

data = []
for plate_number in plate_numbers:
	plate_name = 'keio' + str(plate_number)
	data.append(os.path.join(direc, plate_name))

for root_direc, plate_number in zip(data, plate_numbers):
	print root_direc
	#store all data in this dictionary. Keys = position
	all_data_dict = {}
	#Define directory path to infection data (all positions)
	data_direc = os.path.join(root_direc, 'data')

	#Define directory path to control data (all positions)
	control_direc = os.path.join(root_direc, 'control')

	#Define directory path to where you want to store neural net outputs. 
	#mask directories must exist at run time!
	mask_direc = os.path.join(root_direc, 'masks')
	control_mask_direc = os.path.join(root_direc,'control_masks')

	#Load saved data
	mean_FITC_name = os.path.join(root_direc, 'mean_FITC.pkl')
	mean_cherry_name = os.path.join(root_direc, 'mean_cherry.pkl')
	fitc_dict = pickle.load(open(mean_FITC_name, 'rb'))
	cherry_dict = pickle.load(open(mean_cherry_name, 'rb'))

	wells = sorted(fitc_dict.keys())

	#classify data
	lytic_dict, lysogenic_dict, uninfected_dict = classify_infections_gmm(fitc_dict, cherry_dict, wells = wells, classification_wells = classification_wells)

	#turn unclassified and classified data into a compiled dictionary
	all_data_dict = save_slip_wells(fitc_dict, cherry_dict, plate_number) 
	classified_data_dict = save_classified_wells(lytic_dict, lysogenic_dict, uninfected_dict, plate_number = plate_number)
	file_name = 'classified_data_K' + str(plate_number) + '.pkl'
	all_file_name = 'all_data_K' + str(plate_number) + '.pkl'
	save_direc = os.path.join(direc, 'datatxt')
	save_path = os.path.join(save_direc, file_name)
	all_save_path = os.path.join(save_direc, all_file_name)
	save_file = open(save_path, 'w+')
	all_save_file = open(all_save_path, 'w+')
	pickle.dump(all_data_dict, all_save_file, pickle.HIGHEST_PROTOCOL)
	pickle.dump(classified_data_dict,save_file, pickle.HIGHEST_PROTOCOL)
	save_file.close()
	all_save_file.close()


	