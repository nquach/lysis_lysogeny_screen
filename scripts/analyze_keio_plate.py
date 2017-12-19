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

#Root directory
direc = "/scratch/users/nquach/"

#List of paths to the plates you wish to analyze
data = [os.path.join(direc, 'keio79')]

#List of plate numbers for the plates you wish to analyze
plate_numbers = [79] 

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
	mean_FITC, mean_cherry = analyze_plate(data_direc, mask_direc, pos_list = range(25), wells = wells, panorama = True)
	mean_FITC_name = os.path.join(root_direc, 'mean_FITC.pkl')
	mean_cherry_name = os.path.join(root_direc, 'mean_cherry.pkl')
	pickle.dump(mean_FITC, open(mean_FITC_name, 'wb'))
	pickle.dump(mean_cherry, open(mean_cherry_name, 'wb'))



