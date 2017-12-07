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
from SLIP_functions import analyze_well, analyze_plate, segment_SLIP, plot_slip_well, plot_slip_wells, plot_slip_wells_gmm, plot_slip_wells_lysis_posterior, plot_slip_wells_MOI_posterior
from SLIP_functions import plot_slip_joint_plot, fit_kde, compute_p_values
from SLIP_functions import classify_infections, compute_p_lysis_posterior, compute_MOI_posterior, compute_stats
from keio_names import get_keio_names, pos_to_strain
import seaborn as sns
import pandas as pd
import pymc3 as pm
import json

#Define root directory path
direc = "/scratch/users/nquach/"
plate_numbers = ['hits41'] #, 39]
classification_wells = ['H12']#['G11','G12']   #

#Uncomment if plotting pass1 hit plate.
titles_dict_path = "/scratch/users/nquach/datatxt/hits_pos_to_name.txt"
titles_file = open(titles_dict_path,'r')
titles_dict = json.load(titles_file)
hits_dict = titles_dict['4']

data = []
for plate_number in plate_numbers:
	plate_name = str(plate_number) 
	#Uncomment if plotting a Keio library plate
	#plate_name ='keio' + str(plate_number)
	data.append(os.path.join(direc, plate_name))

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

	#Define which wells were used

	row_data = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
	col_data = range(1,13)#[10, 11, 12]
#

	#Load saved data
	mean_FITC_name = os.path.join(root_direc, 'mean_FITC.pkl')
	mean_cherry_name = os.path.join(root_direc, 'mean_cherry.pkl')
	mean_FITC = pickle.load(open(mean_FITC_name, 'rb'))
	mean_cherry = pickle.load(open(mean_cherry_name, 'rb'))

	#Plot the scatter plot of intensities
	#Uncomment if plotting a plate of the Keio collection
	'''
	wells = []
	titles = []
	keio_names_array = get_keio_names()

	for row in row_data:
		for col in col_data:
			well = row + str(col)
			wells += [well]
			titles += [pos_to_strain(keio_names_array, plate_number, well)]
	'''
	#Uncomment if plotting a plate of the pass1 hits collection
	wells = []
	titles = []
	keio_names_array = get_keio_names()
	i = 0
	for row in row_data:
		for col in col_data:
			if i < 89:
				well = row + str(col)
				wells += [well]
				titles += [hits_dict[well]]
				i = i+1

	print "Plotting unclassified data..."
	plot_slip_wells(mean_FITC, mean_cherry, wells = wells, titles = titles, save_direc = os.path.join(direc, 'scatterplot'), plate_number = plate_number, save_fig = True)
	
	print "Plotting GMM classification..."
	plot_slip_wells_gmm(mean_FITC, mean_cherry, wells = wells, titles = titles, save_direc = os.path.join(direc, 'gmm_scatterplot'), plate_number = plate_number, classification_wells = classification_wells)
	
	print "Plotting Lysis Posteriors..."
	plot_slip_wells_lysis_posterior(mean_FITC, mean_cherry, wells = wells, titles = titles, save_direc = os.path.join(direc, 'lysis_posterior'), plate_number = plate_number, classification_wells = classification_wells)
	print "Plotting MOI Posteriors..."
	plot_slip_wells_MOI_posterior(mean_FITC, mean_cherry, wells = wells, titles = titles, save_direc = os.path.join(direc, 'MOI_posterior'), plate_number = plate_number, classification_wells = classification_wells)
	print "Converting to JSON..."
	compute_stats(mean_FITC, mean_cherry, wells = wells, titles = titles, save_direc = '/scratch/users/nquach/datatxt/', plate_number = plate_number, classification_wells = classification_wells)
	