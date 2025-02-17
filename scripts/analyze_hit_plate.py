'''
analyze_hit_plate.py 
Extracts fluorescence values for a 96 well plate (or plates) of a preliminary hits plate and classifies them.
The classified fluorescence values are saved as a pkl file.
Written by Nicolas Quach
'''
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
from SLIP_functions import analyze_well, analyze_plate, segment_SLIP, plot_slip_well, plot_slip_wells, plot_slip_wells_gmm, plot_slip_wells_lysis_posterior, plot_slip_wells_MOI_posterior
from SLIP_functions import plot_slip_joint_plot, fit_kde, compute_p_values
from SLIP_functions import classify_infections, compute_p_lysis_posterior, compute_MOI_posterior
from keio_names import get_keio_names, pos_to_strain
import seaborn as sns
import pandas as pd
import pymc3 as pm
import json
import utils

ROOT_DIREC = utils.ROOT_DIREC

#Define root directory for the data folders
direc = "/scratch/users/nquach/"
#List of strings of the paths to the data for the hits plate you want to analyze
data = [os.path.join(direc, 'hits31')]
#List of strings of the hits plates you wish to analyze
plate_numbers = ['hits31'] 

#Load the map of gene names for hits plates.
datatxt_direc = os.path.join(ROOT_DIREC, 'datatxt')
titles_dict_path = os.path.join(datatxt_direc, 'hits_pos_to_name.txt')
titles_file = open(titles_dict_path,'r')
titles_dict = json.load(titles_file)
#Define which plate is being analyzed. So hits plate 3 is keyed as '3'
hits_dict = titles_dict['3']

for root_direc, plate_number in zip(data, plate_numbers):
	print root_direc
	#Define directory path to infection data (all positions)
	data_direc = os.path.join(root_direc, 'data')

	#Define directory path to control data (all positions)
	control_direc = os.path.join(root_direc, 'control')

	#Define directory path to where you want to store neural net outputs. 
	#mask directories must exist at run time!
	mask_direc = os.path.join(root_direc, 'masks')
	control_mask_direc = os.path.join(root_direc,'control_masks')
	row_data = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
	col_data = range(1,13)
	wells = []
	for row in row_data:
		for col in col_data:
			well = row + str(col)
			wells += [well]
	# Quantify the data from the infection wells
	mean_FITC, mean_cherry = analyze_plate(data_direc, mask_direc, pos_list = range(16), wells = wells, panorama=False)
	mean_FITC_name = os.path.join(root_direc, 'mean_FITC.pkl')
	mean_cherry_name = os.path.join(root_direc, 'mean_cherry.pkl')
	pickle.dump(mean_FITC, open(mean_FITC_name, 'wb'))
	pickle.dump(mean_cherry, open(mean_cherry_name, 'wb'))

