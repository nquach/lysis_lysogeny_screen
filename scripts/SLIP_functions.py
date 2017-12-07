import matplotlib
matplotlib.use('Agg')
import h5py
import tifffile as tiff
from keras.backend.common import _UID_PREFIXES

from cnn_functions import nikon_getfiles, get_image, run_models_on_directory, get_image_sizes, segment_nuclei, segment_cytoplasm, dice_jaccard_indices, run_model
from model_zoo import sparse_bn_feature_net_61x61 as fn 
from keio_names import get_keio_names, pos_to_strain

import os
import fnmatch
import numpy as np
from skimage.io import imread
from scipy.stats import mode
from skimage.measure import label, regionprops
import scipy
from scipy import stats, ndimage
from matplotlib import pyplot as plt
import seaborn as sns
import pymc3 as pm
import math
from sklearn import mixture
import json
import cPickle as pickle

#Load gene name to position maps
keio_names_array = get_keio_names()
titles_dict_path = "/scratch/users/nquach/datatxt/hits_pos_to_name.txt"
titles_file = open(titles_dict_path,'r')
titles_dict = json.load(titles_file)
hits_dict = titles_dict['4']


def segment_SLIP(data_direc, mask_direc, alphabet = ['A','B','C','D','E','F','G','H'], columns = range(1,13), replace = True, pos_list = range(25)):
	"""
    Function: segment_SLIP()

    Segments all image files within a directory tree corresponding to an entire plate.

    Parameters
    ----------
    data_direc : str
        directory containing image data. Directory tree must be of the format plate->position->images
    mask_direc : str
        directory to which masks will be saved. Directory tree doesn't have to be made at run time. A directory tree similar to that containing the data images will be created.
    alphabet: list
    	rows you wish to segment
    columns: list
    	columns you wish to segment
    replace: bool
    	bool specifying whether or not to resegment images that have already been segmented
    pos_list: list
    	list of position indexes per well position. For an n by m grid of position per well, input range(m*n)


    Returns
    -------
    None

	"""
	print alphabet, columns

	for row in alphabet:
		for column in columns:
			pos = row + str(column)
			print "Segmenting Position " + pos
			current_data_direc = os.path.join(data_direc, pos)
			current_mask_direc = os.path.join(mask_direc, pos)

			#make position directory for masks if necessary
			try:
				os.stat(current_mask_direc)
			except:
				os.mkdir(current_mask_direc)

			file_list = os.listdir(current_mask_direc)
			file_list = fnmatch.filter(file_list, 'feature*')
			number_of_files_if_finished = 3 * len(pos_list)

			if replace is True:
				segment_SLIP_plate(current_data_direc, current_mask_direc)

			else:
				if len(file_list) == number_of_files_if_finished:
					print "Well " + pos + " already segmented ... skipping ..."

				else:
					segment_SLIP_plate(current_data_direc, current_mask_direc)



def segment_SLIP_plate(data_direc, mask_direc):
	"""
    Function: segment_SLIP_plate()

    Helper function for segment_SLIP(). Actually runs images through neural network.

    Parameters
    ----------
    data_direc : str
        directory containing image data. Directory tree must be of the format plate->position->images
    mask_direc : str
        directory to which masks will be saved. Directory tree doesn't have to be made at run time. A directory tree similar to that containing the data images will be created.

    Returns
    -------
    None

	"""
	data_location = data_direc
	phase_location = mask_direc

	phase_channel_names = ['Phase']#['channel000']

	trained_network_phase_directory = "/home/nquach/DeepCell/trained_networks/slip/"   

	phase_prefix = "2017-06-06_slip_61x61_bn_feature_net_61x61_"
	#"2017-02-12_ecoli_90x_31x31_ecoli_90x_feature_net_31x31_"

	win_phase = 30

	image_size_x, image_size_y = get_image_sizes(data_location, phase_channel_names)
	image_size_x /= 2
	image_size_y /= 2

	list_of_phase_weights = []
	for j in xrange(1):
		phase_weights = os.path.join(trained_network_phase_directory,  phase_prefix + str(j) + ".h5")
		list_of_phase_weights += [phase_weights]

	phase_predictions = run_models_on_directory(data_location, phase_channel_names, phase_location, model_fn = fn, 
		list_of_weights = list_of_phase_weights, image_size_x = image_size_x, image_size_y = image_size_y, 
		win_x = win_phase, win_y = win_phase, std = False, split = False)

	
def screen_masks(list_of_masks, confidence = 0.75, debris_size = 20, area_filter = True, eccen_filter = True, minor_axis_filter = False, major_axis_filter = False, solidity_filter = True):
	"""
    Function: screen_masks()

    Helper function for analyze_plate(). Performs post-processing of neural network outputs to convert them into binary masks.

    Parameters
    ----------
   	list_of_masks : list
   		list of mask images as Numpy arrays.
   	confidence : float
   		Confidence cutoff for converting the neural network prediction map into a binary mask
   	debris_size : int
   		Pixel area cutoff for removing ROI's that are too small (and are likely debris)
   	area_filter : bool
   		Bool that determine whether the ROI's are screen for outliers in area (>2 std dev from the mean)
   	eccen_filter : bool
   		Bool that determine whether the ROI's are screen for outliers in eccentricity (>2 std dev from the mean)
   	minor_axis_filter : bool
   		Bool that determine whether the ROI's are screen for outliers in minor axis (>2 std dev from the mean)
   	major_axis_filter : bool
   		Bool that determine whether the ROI's are screen for outliers in major axis (>2 std dev from the mean)
   	solidity_filter : bool
   		Bool that determine whether the ROI's are screen for outliers in solidity (>2 std dev from the mean)
	
    Returns
    -------
    list
    	list of screened masks as binary Numpy arrays.
	"""
	mask_area = []
	mask_ecc = []
	mask_minor_axis = []
	mask_major_axis = []
	mask_solidity = []
	list_of_screened_masks = []

	print len(list_of_masks)
	for mask in list_of_masks:
		print mask
		mask = mask > confidence
		label_mask = label(mask)
		mask_props = regionprops(label_mask, mask)

		for prop in mask_props:
			if prop.area > debris_size:
				mask_area.append(prop.area)
			mask_ecc.append(prop.eccentricity)
			mask_minor_axis.append(prop.minor_axis_length)
			mask_major_axis.append(prop.major_axis_length)
			mask_solidity.append(prop.solidity)

	area_limit = [np.mean(mask_area) - np.std(mask_area), np.mean(mask_area) + 2*np.std(mask_area)]
	ecc_limit = [np.mean(mask_ecc) - np.std(mask_ecc), np.mean(mask_ecc) + np.std(mask_ecc)]
	minor_limit = [np.mean(mask_minor_axis) - np.std(mask_minor_axis), np.mean(mask_minor_axis) + np.std(mask_minor_axis)]
	major_limit = [np.mean(mask_major_axis) - np.std(mask_major_axis), np.mean(mask_major_axis) + np.std(mask_major_axis)]
	solidity_limit = [np.mean(mask_solidity) - np.std(mask_solidity), np.mean(mask_solidity) + np.std(mask_solidity)]

	for mask in list_of_masks:
		mask = mask > confidence
		label_mask = label(mask)
		mask_props = regionprops(label_mask)
		for prop in mask_props:
			if area_filter:
				if prop.area < area_limit[0] or prop.area > area_limit[1]:
					mask[label_mask == prop.label] = 0
			if eccen_filter:
				if prop.eccentricity < ecc_limit[0] or prop.eccentricity > ecc_limit[1]:
					mask[label_mask == prop.label] = 0
			if minor_axis_filter:
				if prop.minor_axis_length < minor_limit[0] or prop.minor_axis_length > minor_limit[1]:
					mask[label_mask == prop.label] = 0
			if major_axis_filter:
				if prop.major_axis_length < major_limit[0] or prop.major_axis_length > major_limit[1]:
					mask[label_mask == prop.label] = 0
			if solidity_filter:
				if prop.solidity < solidity_limit[0] or prop.solidity > solidity_limit[1]:
					mask[label_mask == prop.label] = 0
		list_of_screened_masks.append(mask)

	return list_of_screened_masks


def background_subtraction(image, median = True, mask = None, confidence = 0.75):
	"""
    Function: background_subtraction()

    Helper function for analyze_plate(). Normalizes fluorescence pixel values against the background.

    Parameters
    ----------
   	image : Numpy array
   		Fluorescence image to normalize against background
   	median : bool
   		Bool that determines whether the median background pixel value is used as the normalizing value
   	confidence : float
   		Confidence cutoff for converting the neural network prediction map into a binary mask

    Returns
    -------
    Numpy array
    	Returns background corrected image.
	"""
	if mask is None:
		if median:
			background = np.median(image.flatten())
		else:
			avg_kernel = np.ones((61,61))
			background = ndimage.convolve(image, avg_kernel)/avg_kernel.size
	else:
		if median:
			mask = mask > confidence
			background = np.median(image[mask == 0].flatten())

	return image - background

def analyze_well(data_direc, mask_direc, pos_list, panorama = True):
	"""
    Function: analyze_well()

    Helper function for analyze_plate(). Extracts the mean fluorescence values for each cells given the images and masks for a single well.

    Parameters
    ----------
   	data_direc : str
        directory containing image data. 
    mask_direc : str
        directory to which masks will be saved. 
	pos_list : list of int
		list containing the FOV positions to be analyzed (e.g. range(16))
	panorama : bool
		Bool that turns on panoramic stitching of images. Increases runtime substantially.

    Returns
    -------
    mean_FITC
    	Returns a list of mean GFP channel fluorescence values for cells in the well
    mean_cherry
    	Returns a list of mean Cherry channel fluorescence values for cells in the well
	"""
	if panorama is True:
		list_of_masks = []
		FITC_list = []
		cherry_list = []
		phase_list = []

		for pos in pos_list:
			print 'Position ' + str(pos)
			mask_name = os.path.join(mask_direc, 'feature_1_frame_' + str(pos) + '.tif')
			FITC_name = os.path.join(data_direc, 'img_000000000_EGFP_' +  str(pos).zfill(3) + '.tif')
			cherry_name = os.path.join(data_direc, 'img_000000000_mCherry_' +  str(pos).zfill(3) + '.tif')
			phase_name = os.path.join(data_direc, 'img_000000000_Phase_' +  str(pos).zfill(3) + '.tif')

			print 'Reading image...'
			mask = np.float32(imread(mask_name))[40:-40, 140:-140]
			FITC = np.float32(imread(FITC_name))[40:-40, 140:-140]
			cherry = np.float32(imread(cherry_name))[40:-40, 140:-140]
			phase = np.float32(imread(phase_name))[40:-40, 140:-140]

			FITC = background_subtraction(FITC, mask = mask)
			cherry = background_subtraction(cherry, mask = mask)

			list_of_masks.append(mask)
			FITC_list.append(FITC)
			cherry_list.append(cherry)
			phase_list.append(phase)
		print 'Stitching...'
		# Check the stitching parameters - if off, use pre computed stitching parameters
		mask_pan, h, v = merge_images_v2(list_of_masks)
		
		print h, v

		if len(pos_list) == 9:
			h_pre = [490, 498, 491, 499, 490, 499, 10, 11]
			v_pre = [9, 9, 9, 9, 9, 8, 595, 596]

		if len(pos_list) == 25:
			h_pre = [491, 497, 493, 493, 492, 497, 494, 492, 491, 497, 493, 493, 492, 497, 493, 492, 491, 497, 493, 492, 10, 11, 11, 11]
			v_pre = [10, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 591, 610, 585, 596]

		if len(pos_list) == 16:
			h_pre = [894, 209, 884, 644, 860, 943, 145, 442, 907, 731, 913, 851, 1, 1, 1]
			v_pre = [206, 796, 133, 874, 90, 971, 757, 894, 289, 410, 559, 23, 1, 1, 1]

		replace = False  #CHANGE ME TO FALSE FOR NON PRECOMPUTED PARAMS
		for h_c, h_p, v_c, v_p in zip(h, h_pre, v, v_pre):
			if np.abs(h_c-h_p) > 10 or np.abs(h_c-h_p) > 10:
				replace = True

		if replace == True:
			h = h_pre
			v = v_pre
		
		print 'Making mask panos...'
		mask_panorama = merge_images_v2(list_of_masks, h = h, v = v)[0]
		print mask_panorama
		mask_panorama = screen_masks([mask_panorama])[0]
		print 'Making phase panos...'
		phase_panorama = merge_images_v2(phase_list, h = h, v = v)[0]
		print 'Making FITC panos...'
		fitc_panorama = merge_images_v2(FITC_list, h = h, v = v)[0]
		print 'Making cherry panos...'
		cherry_panorama = merge_images_v2(cherry_list, h = h, v = v)[0]
		print 'Saving panos...'
		# Save panoramic images
		mask_name = os.path.join(mask_direc, 'panorama_mask.tif')
		phase_name = os.path.join(mask_direc, 'panorama_phase.tif')
		fitc_name = os.path.join(mask_direc, 'panorama_fitc.tif')
		cherry_name = os.path.join(mask_direc, 'panorama_cherry.tif')

		tiff.imsave(mask_name, np.float32(mask_panorama))
		tiff.imsave(phase_name, np.float32(phase_panorama))
		tiff.imsave(fitc_name, np.float32(fitc_panorama))
		tiff.imsave(cherry_name, np.float32(cherry_panorama))


		# Collect data points
		mean_FITC = []
		mean_cherry = []

		label_mask = label(mask_panorama)
		FITC_props = regionprops(label_mask, fitc_panorama)
		cherry_props = regionprops(label_mask, cherry_panorama)

		for props in FITC_props:
			mean_FITC.append(props.mean_intensity)

		for props in cherry_props:
			mean_cherry.append(props.mean_intensity)


	if panorama is False:
		list_of_masks = []
		FITC_list = []
		cherry_list = []
		for pos in pos_list:
			mask_name = os.path.join(mask_direc, 'feature_1_frame_' + str(pos) + '.tif')
			FITC_name = os.path.join(data_direc, 'img_000000000_EGFP_' +  str(pos).zfill(3) + '.tif')
			cherry_name = os.path.join(data_direc, 'img_000000000_mCherry_' +  str(pos).zfill(3) + '.tif')

			mask = np.float32(imread(mask_name))
			FITC = np.float32(imread(FITC_name))
			cherry = np.float32(imread(cherry_name))

			FITC_norm = background_subtraction(FITC)
			cherry_norm = background_subtraction(cherry)

			list_of_masks.append(mask)
			FITC_list.append(FITC_norm)
			cherry_list.append(cherry_norm)

		list_of_screened_masks = screen_masks(list_of_masks)

		# Save screened masks
		for pos, mask in zip(pos_list, list_of_screened_masks):
			mask_name = os.path.join(mask_direc, 'mask_' + str(pos) + '.tif')
			tiff.imsave(mask_name, np.float32(mask))

		# Collect data points
		mean_FITC = []
		mean_cherry = []

		for mask, FITC, cherry in zip(list_of_screened_masks, FITC_list, cherry_list):
			label_mask = label(mask)

			FITC_props = regionprops(label_mask, FITC)
			cherry_props = regionprops(label_mask, cherry)

			for props in FITC_props:
				mean_FITC.append(props.mean_intensity)

			for props in cherry_props:
				mean_cherry.append(props.mean_intensity)

	return mean_FITC, mean_cherry

def analyze_plate(data_direc, mask_direc, pos_list, wells, panorama = True):
	"""
    Function: analyze_plate()

    Extracts the mean fluorescence values for each cells given the images and masks for all wells in a 96 plate.

    Parameters
    ----------
   	data_direc : str
        directory containing image data. 
    mask_direc : str
        directory to which masks will be saved. 
	pos_list : list of str
		list containing the positions to be analyzed (e.g. ['A1','A2'])
	panorama : bool
		Bool that turns on panoramic stitching of images. Increases runtime substantially.

    Returns
    -------
    mean_FITC
    	Returns a list of mean GFP channel fluorescence values for cells in the well
    mean_cherry
    	Returns a list of mean Cherry channel fluorescence values for cells in the well
	"""
	mean_FITC = {}
	mean_cherry = {}
	
	for well in wells:
		print 'Processing well ' + well
		data_directory = os.path.join(data_direc, well)
		mask_directory = os.path.join(mask_direc, well)
		fitc, cherry = analyze_well(data_directory, mask_directory, pos_list, panorama=panorama)
		mean_FITC[well] = fitc
		mean_cherry[well] = cherry
	return mean_FITC, mean_cherry

def plot_slip_wells(fitc_dict, cherry_dict, wells, titles, infected_cells = None, plate_number = None, save_direc = '/home/vanvalen/keio_screen/scatter_plots/', save_fig = False):
	"""
    Function: plot_slip_wells()

    Given mean fluorescence data for all cells in each well of a 96 well plate, plots scatterplots of fluorescence data for each well.

    Parameters
    ----------
   	fitc_dict : dict
        dict of lists containing fluorescence data for each well. key = position (e.g. 'A1'), element = list of FITC fluorescence data
    cherry_dict : dict
        dict of lists containing fluorescence data for each well. key = position (e.g. 'A1'), element = list of Cherry fluorescence data
	titles : list
		list containing the gene names corresponding to each well
	infected_cells : bool
		Bool that turns on panoramic stitching of images. Increases runtime substantially.
	plate_number : int
		Keio plate number
	save_direc : str
		String defining file path to where you want to save the plot
	save_fig : bool
		Bool that determines whether the plot is saved as a .pdf file in save_direc.

    Returns
    -------
    fig
    	scatterplot of the fluorescence data for all 96 wells.
    axes
    	axes of the plot
	"""

	fig, axes = plt.subplots(8,12, figsize = (4*12, 4*8))

	xmax = 0
	xmin = 0
	ymax = 0
	ymin = 0

	for well in wells:
		if infected_cells is None:
			fitc_list = fitc_dict[well]
			cherry_list = cherry_dict[well]
		else:
			fitc_list = np.array(fitc_dict[well])[infected_cells[well]]
			cherry_list = np.array(cherry_dict[well])[infected_cells[well]]

	for well, title in zip(wells, titles):
		if infected_cells is None:
			fitc_list = fitc_dict[well]
			cherry_list = cherry_dict[well]
		else:
			fitc_list = np.array(fitc_dict[well])[infected_cells[well]]
			cherry_list = np.array(cherry_dict[well])[infected_cells[well]]

		if len(fitc_list) > 0:
			alphabet = ['A','B','C','D','E','F','G','H']
			chars = list(well)
			row = alphabet.index(chars[0])
			if len(chars) == 2:
				column = int(chars[1])-1
			if len(chars) == 3:
				column = int(chars[1] + chars[2])-1

			axes[row,column].plot(fitc_list, cherry_list,'o')
			axes[row,column].set_xlim([min(fitc_list), max(fitc_list)])
			axes[row,column].set_ylim([min(cherry_list), max(cherry_list)]) 
			axes[row,column].set_xlabel('FITC Pixel Intensity (au)')
			axes[row,column].set_ylabel('Cherry Pixel Intensity (au)')
			axes[row,column].set_title(title)
	plt.tight_layout()

	if save_fig:
		figure_name = os.path.join(save_direc, 'raw_data_plate_' + str(plate_number) + '.pdf')
		fig.savefig(figure_name, bbox_inches = 'tight')
	return fig, axes

def save_slip_wells(fitc_dict, cherry_dict, plate_number = None, is_hits = False):
	"""
    Function: save_slip_wells()

    Given the fluorescence data for a particular plate, returns a dict of fluorescence data organized by position.  

    Parameters
    ----------
   	fitc_dict : dict
        dict of lists containing fluorescence data for each well. key = position (e.g. 'A1'), element = list of FITC fluorescence data
    cherry_dict : dict
        dict of lists containing fluorescence data for each well. key = position (e.g. 'A1'), element = list of Cherry fluorescence data
	plate_number : int
		Keio plate number
	is_hits : bool
		Bool. Is true if the plate being analyzed is from the verification round.

    Returns
    -------
    data_dict : dict
    	Dict containing fluorescence data. keys = gene name (e.g. 'hflK'), item = list of lists of the form [fitc_list, cherry_list]
	"""
	if type(plate_number) is str:
		plate_num = int(str(plate_number)[0])
	else:
		plate_num = plate_number
	
	wells = sorted(fitc_dict.keys())
	titles = []
	print wells
	for well in wells:
		if is_hits:
			titles += [hits_dict[well]]
		else:
			titles += [pos_to_strain(keio_names_array, plate_num, well)]
	
	data_dict = {}

	for well, title in zip(wells, titles):
		fitc_list = fitc_dict[well]
		cherry_list = cherry_dict[well]
		print well, title
		if title in data_dict.keys():
			data_dict[title][0] = data_dict[title][0] + [fitc_list]   
			data_dict[title][1] = data_dict[title][1] + [cherry_list]
		else:
			data_dict[title] = [fitc_list, cherry_list]
	return data_dict

def save_classified_wells(lytic_dict, lyso_dict, uninfected_dict, plate_number = None, is_hits = False):
	"""
    Function: save_classified_wells()

    Given the classified fluorescence data for a particular plate, returns a dict of fluorescence data organized by position.  

    Parameters
    ----------
   	lytic_dict : dict
        dict of lists containing lytic fluorescence data for each well. key = position (e.g. 'A1'), element = list of lists in the form [fitc_data, cherry_data]
    lyso_dict : dict
        dict of lists containing lysogenic fluorescence data for each well. key = position (e.g. 'A1'), element = list of lists in the form [fitc_data, cherry_data]
	uninfected_dict : dict
        dict of lists containing uninfected fluorescence data for each well. key = position (e.g. 'A1'), element = list of lists in the form [fitc_data, cherry_data]
	plate_number : int
		Keio plate number
	is_hits : bool
		Bool. Is true if the plate being analyzed is from the verification round.

    Returns
    -------
    data_dict : dict
    	Dict containing fluorescence data. keys = gene name (e.g. 'hflK'), item = list of list of lists of the form [lytic_data([fitc_data, cherry_data]), lyso_data([fitc_data, cherry_data]), uninfected_data ([fitc_data, cherry_data])] 
	"""
	if type(plate_number) is str:
		plate_num = int(str(plate_number)[0])
	else:
		plate_num = plate_number
	
	wells = sorted(lytic_dict.keys())
	print wells
	titles = []

	for well in wells:
		if is_hits:
			titles += [hits_dict[well]]
		else:
			titles += [pos_to_strain(keio_names_array, plate_num, well)]

	data_dict = {}

	for well, title in zip(wells, titles):
		print well, title
		if title in data_dict.keys():
			data_dict[title][0] = data_dict[title][0] + list(lytic_dict[well])
			data_dict[title][1] = data_dict[title][1] + list(lyso_dict[well])
			data_dict[title][2] = data_dict[title][2] + list(uninfected_dict[well])
		else:
			data_dict[title] = [list(lytic_dict[well]), list(lyso_dict[well]), list(uninfected_dict[well])]
	
	return data_dict

def classify_infections_gmm(fitc_dict, cherry_dict, wells, classification_wells = None):
	"""
    Function: classify_infections_gmm()

    Helper function. Given the dicts of cell fluorescence values, returns dicts of the fluorescences value classified into lytic, lysogenic and uninfected dicts. 

    Parameters
    ----------
   	fitc_dict : dict
        dict of lists containing GFP fluorescence data for each well. key = position (e.g. 'A1'), element = list of GFP fluorescence values
    cherry_dict : dict
        dict of lists containing Cherry fluorescence data for each well. key = position (e.g. 'A1'), element = list of Cherry fluorescence values
	wells : list
        list of well positions to classify (e.g. ['A1','A2'])
	classification_wells : list
		list of well positions containing data to train GMM model if you don't want to train the model on the aggregated plate data.
	
    Returns
    -------
    lytic_dict : dict
    	dict of lists containing lytic fluorescence data for each well. key = position (e.g. 'A1'), element = list of lists in the form [fitc_data, cherry_data]
    lyso_dict : dict
    	dict of lists containing lysogenic fluorescence data for each well. key = position (e.g. 'A1'), element = list of lists in the form [fitc_data, cherry_data]
    uninfected_dict : dict
    	dict of lists containing uninfected fluorescence data for each well. key = position (e.g. 'A1'), element = list of lists in the form [fitc_data, cherry_data]
  	"""
	fitc_list = []
	cherry_list = []

	if classification_wells is None:
		classification_wells = wells

	for well in classification_wells:
		if len(fitc_dict[well]) != 0:
			fitc_list += fitc_dict[well]
		if len(cherry_dict[well]) != 0:
			cherry_list += cherry_dict[well]
	data = np.stack((fitc_list,cherry_list), axis = 1)

	gmm = mixture.GMM(n_components = 3).fit(data)

	#identify populations
	means = gmm.means_

	fitc_means = np.array([mean[0] for mean in means])
	cherry_means = np.array([mean[1] for mean in means])
	index_list = [0,1,2]
	
	lysogenic_id = np.argmax(cherry_means)
	lytic_id = np.argmax(fitc_means)
	print lytic_id, lysogenic_id
	print fitc_means, cherry_means

	index_list.remove(lytic_id)
	index_list.remove(lysogenic_id)
	uninfected_id = index_list[0]
	
	lytic_dict = {}
	lysogenic_dict = {}
	uninfected_dict = {}

	for well in wells:
		fitc_list = np.array(fitc_dict[well])
		cherry_list = np.array(cherry_dict[well])

		if fitc_list.size == 0:
			continue

		data_well = np.stack((fitc_list, cherry_list), axis = 1)
		probs = gmm.predict_proba(data_well)

		lytic_dict[well] = []
		lysogenic_dict[well] = []
		uninfected_dict[well] = []

		if fitc_list.size == 0 or cherry_list.size == 0:
			lytic_dict[well] = np.array(lytic_dict[well])
			lysogenic_dict[well] = np.array(lysogenic_dict[well])
			uninfected_dict[well] = np.array(uninfected_dict[well])
			continue

		for j in xrange(data_well.shape[0]):

			if probs[j][lytic_id] > 0.95:
				if data_well[j,0] > means[uninfected_id][0]:
					lytic_dict[well] += [data_well[j,:]]
				else:
					uninfected_dict[well] += [data_well[j,:]]

			elif probs[j][lysogenic_id] > 0.95:
				if data_well[j,1] > means[uninfected_id][1]:
					lysogenic_dict[well] += [data_well[j,:]]
				else:
					uninfected_dict[well] += [data_well[j,:]]

			elif probs[j][uninfected_id] > 0.95:
				uninfected_dict[well] += [data_well[j,:]]


		lytic_dict[well] = np.array(lytic_dict[well])
		lysogenic_dict[well] = np.array(lysogenic_dict[well])
		uninfected_dict[well] = np.array(uninfected_dict[well])

	return lytic_dict, lysogenic_dict, uninfected_dict

def classify_infections_gmm2(fitc_dict, cherry_dict, wells, classification_wells = None):
	"""
    Function: classify_infections_gmm2()

    Helper function. Given the dicts of cell fluorescence values, returns dicts of the fluorescences value classified into lytic, lysogenic and uninfected dicts. 
    Only uses plate 29 with H12 as classification well as training data for GMM model.

    Parameters
    ----------
   	fitc_dict : dict
        dict of lists containing GFP fluorescence data for each well. key = position (e.g. 'A1'), element = list of GFP fluorescence values
    cherry_dict : dict
        dict of lists containing Cherry fluorescence data for each well. key = position (e.g. 'A1'), element = list of Cherry fluorescence values
	wells : list
        list of well positions to classify (e.g. ['A1','A2'])
	classification_wells : list
		list of well positions containing data to train GMM model if you don't want to train the model on the aggregated plate data. None=means aggregate plate data.
	
    Returns
    -------
    lytic_dict : dict
    	dict of lists containing lytic fluorescence data for each well. key = position (e.g. 'A1'), element = list of lists in the form [fitc_data, cherry_data]
    lyso_dict : dict
    	dict of lists containing lysogenic fluorescence data for each well. key = position (e.g. 'A1'), element = list of lists in the form [fitc_data, cherry_data]
    uninfected_dict : dict
    	dict of lists containing uninfected fluorescence data for each well. key = position (e.g. 'A1'), element = list of lists in the form [fitc_data, cherry_data]
  	"""
	root_direc = "/scratch/users/nquach/keio29"

	mean_FITC_control_name = os.path.join(root_direc, 'mean_FITC.pkl')
	mean_cherry_control_name = os.path.join(root_direc, 'mean_cherry.pkl')
	mean_FITC_control = pickle.load(open(mean_FITC_control_name, 'rb'))
	mean_cherry_control = pickle.load(open(mean_cherry_control_name, 'rb'))

	control_fitc_list = []
	control_cherry_list = []

	if classification_wells is None:
		classification_wells = wells

	for well in classification_wells:
		if len(mean_FITC_control[well]) != 0:
			control_fitc_list += mean_FITC_control[well]
		if len(mean_cherry_control[well]) != 0:
			control_cherry_list += mean_cherry_control[well]
	control_data = np.stack((control_fitc_list,control_cherry_list), axis = 1)

	gmm = mixture.GMM(n_components = 3).fit(control_data)


	fitc_list = []
	cherry_list = []

	#identify populations
	means = gmm.means_

	fitc_means = np.array([mean[0] for mean in means])
	cherry_means = np.array([mean[1] for mean in means])
	index_list = [0,1,2]
	
	lysogenic_id = np.argmax(cherry_means)
	lytic_id = np.argmax(fitc_means)
	print lytic_id, lysogenic_id
	print fitc_means, cherry_means

	index_list.remove(lytic_id)
	index_list.remove(lysogenic_id)
	uninfected_id = index_list[0]
	
	lytic_dict = {}
	lysogenic_dict = {}
	uninfected_dict = {}

	for well in wells:
		fitc_list = np.array(fitc_dict[well])
		cherry_list = np.array(cherry_dict[well])

		if fitc_list.size == 0:
			continue

		data_well = np.stack((fitc_list, cherry_list), axis = 1)
		probs = gmm.predict_proba(data_well)

		lytic_dict[well] = []
		lysogenic_dict[well] = []
		uninfected_dict[well] = []

		if fitc_list.size == 0 or cherry_list.size == 0:
			lytic_dict[well] = np.array(lytic_dict[well])
			lysogenic_dict[well] = np.array(lysogenic_dict[well])
			uninfected_dict[well] = np.array(uninfected_dict[well])
			continue

		for j in xrange(data_well.shape[0]):

			if probs[j][lytic_id] > 0.95:
				if data_well[j,0] > means[uninfected_id][0]:
					lytic_dict[well] += [data_well[j,:]]
				else:
					uninfected_dict[well] += [data_well[j,:]]

			elif probs[j][lysogenic_id] > 0.95:
				if data_well[j,1] > means[uninfected_id][1]:
					lysogenic_dict[well] += [data_well[j,:]]
				else:
					uninfected_dict[well] += [data_well[j,:]]

			elif probs[j][uninfected_id] > 0.95:
				uninfected_dict[well] += [data_well[j,:]]


		lytic_dict[well] = np.array(lytic_dict[well])
		lysogenic_dict[well] = np.array(lysogenic_dict[well])
		uninfected_dict[well] = np.array(uninfected_dict[well])

	return lytic_dict, lysogenic_dict, uninfected_dict

def compute_stats(fitc_dict, cherry_dict, wells, titles, save_direc = '/scratch/users/nquach/datatxt/', plate_number = 9, classification_wells = None):
	"""
    Function: compute_stats()

    Given the dicts of cell fluorescence values, computes the cell-fate statistics of every strain in a plate and saves the data in a .txt file configured in JSON format 

    Parameters
    ----------
   	fitc_dict : dict
        dict of lists containing GFP fluorescence data for each well. key = position (e.g. 'A1'), element = list of GFP fluorescence values
    cherry_dict : dict
        dict of lists containing Cherry fluorescence data for each well. key = position (e.g. 'A1'), element = list of Cherry fluorescence values
	wells : list
        list of well positions to classify (e.g. ['A1','A2'])
	titles : list
		list of gene names corresponding to the positions listed in wells.
	save_direc : str
		String containing the full path to the directory where you want to save the .txt JSON file.
	plate_number : int or str
		Int or str of which plate of the Keio collection is being analyzed (e.g. 9 or 'hits11')
	classification_wells : list
		list of well positions containing data to train GMM model if you don't want to train the model on the aggregated plate data. None=means aggregate plate data.
	
    Returns
    -------
   	None : None
   		Writes cell fate statistics for all strains in a plate into a .txt file in JSON format. Data is stored as a dict of lists. 
   		key=position, item = list([gene name, number of total cells, number of infected cels, number of uninfected cells, number of lytic cells, number of lysogenic cells, lysis ratio, variance of lysis ratio estimate])
  	"""
	lytic_dict, lysogenic_dict, uninfected_dict = classify_infections_gmm2(fitc_dict, cherry_dict, wells = wells, classification_wells = classification_wells)
	stats_file_name = "keio" + str(plate_number) + ".txt"
	stats_file = open(os.path.join(save_direc, stats_file_name), "w")
	hits_file_name = "hits_keio" + str(plate_number) + ".txt"
	hits_file = open(os.path.join(save_direc, hits_file_name), "w")
	print "Processing Keio " + str(plate_number)
	stats_dict = {}
	hits_dict = {}
	for well, title in zip(wells, titles):
		print well
		if well not in lytic_dict.keys():
			continue
		num_lytic = lytic_dict[well].shape[0]
		num_lyso = lysogenic_dict[well].shape[0]
		num_uninfected = uninfected_dict[well].shape[0]
		total_cells = num_lytic + num_lyso + num_uninfected
		num_infected = num_lytic + num_lyso	
		a = float(1+num_lytic)
		b = float(1+num_lyso)
		mode = 0.5
		if a+b > 2:
			mode = (a-1)/(a+b-2)
		var = a*b/(np.square(a+b)*(a+b+1))
		stats_dict[well] = [title, total_cells, num_infected, num_uninfected, num_lytic, num_lyso, mode, var]
		#print stats_dict
		if mode > 0.8 or mode < 0.2:
			if title != None:
				hits_dict[well] = [title, total_cells, num_infected, num_uninfected, num_lytic, num_lyso, mode, var]

	json.dump(stats_dict, stats_file)
	json.dump(hits_dict, hits_file)
	stats_file.close()
	hits_file.close()
	return None

def plot_slip_wells_gmm(fitc_dict, cherry_dict, wells, titles, save_direc = '/scratch/users/nquach/scatterplot_gmm', plate_number = 9, save_fig = True, classification_wells = None):
	"""
    Function: plot_slip_wells_gmm()

    Given the dicts of cell fluorescence values, plots the classified fluorescence data as a scatterplot for all strains in a plate. 

    Parameters
    ----------
   	fitc_dict : dict
        dict of lists containing GFP fluorescence data for each well. key = position (e.g. 'A1'), element = list of GFP fluorescence values
    cherry_dict : dict
        dict of lists containing Cherry fluorescence data for each well. key = position (e.g. 'A1'), element = list of Cherry fluorescence values
	wells : list
        list of well positions to classify (e.g. ['A1','A2'])
	titles : list
		list of gene names corresponding to the positions listed in wells.
	save_direc : str
		String containing the full path to the directory where you want to save the .txt JSON file.
	plate_number : int or str
		Int or str of which plate of the Keio collection is being analyzed (e.g. 9 or 'hits11')
	save_fig : bool
		bool determining whether the plotted figure is saved. Defaults to TRUE=save figure
	classification_wells : list
		list of well positions containing data to train GMM model if you don't want to train the model on the aggregated plate data. None=means aggregate plate data.
	
    Returns
    -------
    None : None
    	Saves a scatterplot of the classified fluorescence data for all strains in a plate.
  	"""
	sns.set_context('notebook', font_scale = 1.1)
	sns.set_style('white')
	sns.set_style('ticks')

	sky_blue = (86./255., 180./255., 233./255.)
	bluish_green = (0, 158./255., 115./255.)
	reddish_purple = (204./255.,121./255.,167./255.)
	black = (0.,0.,0.)

	lytic_dict, lysogenic_dict, uninfected_dict = classify_infections_gmm2(fitc_dict, cherry_dict, wells = wells, classification_wells = classification_wells)
	fig, axes = plt.subplots(8,12, figsize = (4*12, 4*8))

	for well, title in zip(wells, titles):
 		print well
 		if well not in lytic_dict.keys():
			continue

		alphabet = ['A','B','C','D','E','F','G','H']
		chars = list(well)
		row = alphabet.index(chars[0])
		if len(chars) == 2:
			column = int(chars[1])-1
		if len(chars) == 3:
			column = int(chars[1] + chars[2])-1

		if lytic_dict[well].shape[0] > 0:
			axes[row,column].scatter(lytic_dict[well][:,0], lytic_dict[well][:,1], c = bluish_green, alpha = 0.75, edgecolors = 'none')

		if lysogenic_dict[well].shape[0] > 0:
			axes[row,column].scatter(lysogenic_dict[well][:,0], lysogenic_dict[well][:,1], c = reddish_purple, alpha = 0.75, edgecolors = 'none')

		if uninfected_dict[well].shape[0] > 0:
			axes[row,column].scatter(uninfected_dict[well][:,0], uninfected_dict[well][:,1], c = black, alpha = 0.75, edgecolors = 'none')

		axes[row,column].set_xlabel('Lytic Reporter Intensity (au)')
		axes[row,column].set_ylabel('Lysogenic Reporter Intensity (au)')
		axes[row,column].set_xmargin(0.1)
		axes[row,column].set_ymargin(0.1)
		axes[row,column].autoscale(enable = True, axis = 'both', tight = True)
		axes[row,column].set_title(title)
		sns.despine()

	plt.tight_layout()

	if save_fig:
		figure_name = os.path.join(save_direc, 'classified_infections_plate_' + str(plate_number) + '.pdf')
		fig.savefig(figure_name, bbox_inches = 'tight')

	return None

def plot_slip_wells_lysis_posterior(fitc_dict, cherry_dict, wells, titles, save_direc = '/scratch/users/nquach/lysis_posteriors/', plate_number = 9, save_fig = True, classification_wells = None):
	"""
    Function: plot_slip_lysis_posterior()

    Given the dicts of cell fluorescence values, plots the lysis posterior distributions for all strains in a plate. 

    Parameters
    ----------
   	fitc_dict : dict
        dict of lists containing GFP fluorescence data for each well. key = position (e.g. 'A1'), element = list of GFP fluorescence values
    cherry_dict : dict
        dict of lists containing Cherry fluorescence data for each well. key = position (e.g. 'A1'), element = list of Cherry fluorescence values
	wells : list
        list of well positions to classify (e.g. ['A1','A2'])
	titles : list
		list of gene names corresponding to the positions listed in wells.
	save_direc : str
		String containing the full path to the directory where you want to save the .txt JSON file.
	plate_number : int or str
		Int or str of which plate of the Keio collection is being analyzed (e.g. 9 or 'hits11')
	save_fig : bool
		Bool determining whether the plotted figure is saved. Defaults to TRUE=save figure
	classification_wells : list
		list of well positions containing data to train GMM model if you don't want to train the model on the aggregated plate data. None=means aggregate plate data.
	
    Returns
    -------
    None : None
    	Saves a plot of the lysis posterior distributions for all strains in a plate.
  	"""
	sns.set_context('notebook', font_scale = 1.1)
	sns.set_style('white')
	sns.set_style('ticks')

	lytic_dict, lysogenic_dict, uninfected_dict = classify_infections_gmm2(fitc_dict, cherry_dict, wells = wells, classification_wells = classification_wells)
	fig, axes = plt.subplots(8,12, figsize = (4*12, 4*8))

	for well, title in zip(wells, titles):
		print well
		if well not in lytic_dict.keys():
			continue

		alphabet = ['A','B','C','D','E','F','G','H']
		chars = list(well)
		row = alphabet.index(chars[0])
		if len(chars) == 2:
			column = int(chars[1])-1
		if len(chars) == 3:
			column = int(chars[1] + chars[2])-1

		if lytic_dict[well].shape[0] > 0:
			N_lytic = lytic_dict[well].shape[0]
		else:
			N_lytic = 0	

		if lysogenic_dict[well].shape[0] > 0:
			N_lysogenic = lysogenic_dict[well].shape[0]
		else:
			N_lysogenic = 0	

		x, posterior = compute_p_lysis_posterior(N_lytic, N_lysogenic)

		axes[row, column].plot(x, posterior, color = 'k', linewidth = 2)
		axes[row, column].set_xlim([0, 1])
		axes[row, column].set_xlabel('Probability of lysis')
		axes[row, column].set_ylabel('Probability density')
		axes[row, column].set_title(title)

		sns.despine()

	plt.tight_layout()

	if save_fig:
		figure_name = os.path.join(save_direc, 'p_lysis_posterior_plate_' + str(plate_number) + '.pdf')
		fig.savefig(figure_name, bbox_inches = 'tight')

	return None

def compute_inverse_MOI_posterior(N_infected, N_cells):
	"""
    Function: compute_inverse_MOI_posterior()

    Given the number of infected cells and the total number of cells for a particular strain, returns points belonging to the MOI posterior distribution. 

    Parameters
    ----------
   	N_infected : int
   		number of infected cells for a particular strain
   	N_cells : int
   		number of total cells for a particular strain

    Returns
    -------
    x : list of floats
    	list of MOI values for x axis.
    posterior : list of floats
    	list of points on the MOI posterior distribution
  	"""
    x = np.linspace(0,5,200)
    gamma = np.float(N_cells)*np.log(1-1/np.float(N_cells))
    posterior = np.abs(gamma*np.exp(gamma*x))*scipy.stats.beta.pdf(np.exp(gamma*x), 1+N_cells-N_infected, 1+N_infected)

    return x, posterior

def plot_slip_wells_MOI_posterior(fitc_dict, cherry_dict, wells, titles, save_direc = '/scratch/users/nquach/MOI_posteriors/', plate_number = 9, save_fig = True, classification_wells = None):
	"""
    Function: plot_slip_wells_MOI_posterior()

    Given the dicts of cell fluorescence values, plots the MOI posterior distributions for all strains in a plate. 

    Parameters
    ----------
   	fitc_dict : dict
        dict of lists containing GFP fluorescence data for each well. key = position (e.g. 'A1'), element = list of GFP fluorescence values
    cherry_dict : dict
        dict of lists containing Cherry fluorescence data for each well. key = position (e.g. 'A1'), element = list of Cherry fluorescence values
	wells : list
        list of well positions to classify (e.g. ['A1','A2'])
	titles : list
		list of gene names corresponding to the positions listed in wells.
	save_direc : str
		String containing the full path to the directory where you want to save the .txt JSON file.
	plate_number : int or str
		Int or str of which plate of the Keio collection is being analyzed (e.g. 9 or 'hits11')
	save_fig : bool
		Bool determining whether the plotted figure is saved. Defaults to TRUE=save figure
	classification_wells : list
		list of well positions containing data to train GMM model if you don't want to train the model on the aggregated plate data. None=means aggregate plate data.
	
    Returns
    -------
    None : None
    	Saves a plot of the MOI posterior distributions for all strains in a plate.
  	"""
	sns.set_context('notebook', font_scale = 1.1)
	sns.set_style('white')
	sns.set_style('ticks')

	lytic_dict, lysogenic_dict, uninfected_dict = classify_infections_gmm2(fitc_dict, cherry_dict, wells = wells, classification_wells = classification_wells)
	fig, axes = plt.subplots(8,12, figsize = (4*12, 4*8))

	for well, title in zip(wells, titles):
		print well
		if well not in lytic_dict.keys():
			continue
		alphabet = ['A','B','C','D','E','F','G','H']
		chars = list(well)
		row = alphabet.index(chars[0])
		if len(chars) == 2:
			column = int(chars[1])-1
		if len(chars) == 3:
			column = int(chars[1] + chars[2])-1

		if uninfected_dict[well].shape[0] > 0:
			N_uninfected = uninfected_dict[well].shape[0]
		else:
			N_uninfected = 0

		if lytic_dict[well].shape[0] > 0:
			N_lytic = lytic_dict[well].shape[0]
		else:
			N_lytic = 0	

		if lysogenic_dict[well].shape[0] > 0:
			N_lysogenic = lysogenic_dict[well].shape[0]
		else:
			N_lysogenic = 0	

		if N_uninfected == 0:
			continue

		x, posterior = compute_inverse_MOI_posterior(len(fitc_dict[well]) - N_uninfected, len(fitc_dict[well]))

		axes[row, column].plot(x, posterior, color = 'k', linewidth = 2)
		axes[row, column].set_xlim([0, 5])
		axes[row, column].set_xlabel('MOI')
		axes[row, column].set_ylabel('Probability density')
		axes[row, column].set_title(title)
		sns.despine()

	plt.tight_layout()

	if save_fig:
		figure_name = os.path.join(save_direc, 'MOI_posterior_plate_' + str(plate_number) + '.pdf')
		fig.savefig(figure_name, bbox_inches = 'tight')

	return None

def compute_p_lysis_posterior(N_lysis, N_lysogeny):
	"""
    Function: compute_p_lysis_posterior()

    Given the number of lytic and lysogenic events, computes the lysis posterior for a particular strain.

    Parameters
    ----------
   	N_lysis : int
   		number of lytic events for a particular strain
   	N_lysogeny: int
   		number of lysogenic events for a particular strain
	
    Returns
    -------
    x : list of floats
    	list of probability of lysis values
    posterior : list of floats
    	list of points discribing the probability of lysis posterior distribution
  	"""
    x = np.linspace(0,1,100)
    posterior= scipy.stats.beta.pdf(x, 1+N_lysis, 1+N_lysogeny)
    return x, posterior

def cross_corr(im0, im1):
	"""
    Function: cross_corr()

    Helper function for merge_images(). Given two matrices, compute the stitching parameters t0, and t1

    Parameters
    ----------
   	im0 : Numpy array
   		first image to stitch
   	im1 : Numpy array
   		second image to stitch
	
    Returns
    -------
    t0 : int
    	stitching parameter
    t1 : int
    	stitching parameter
  	"""
	from numpy.fft import fft2, ifft2
	f0 = fft2(im0)
	f1 = fft2(im1)
	shape = im0.shape
	ir = abs(ifft2((f0 * f1.conjugate()) / (abs(f0) * abs(f1))))
	t0, t1 = np.unravel_index(np.argmax(ir), shape)
	return t0, t1

def is_square(integer):
	"""
    Function: is_square()

    Helper function for merge_images(). Given an integer, determines if it is a square number

    Parameters
    ----------
   	integer : int
   		integer to test
	
    Returns
    -------
    bool : bool
		Bool value. True if the number is square. 
  	"""
	root = math.sqrt(integer)
	if int(root + 0.5) ** 2 == integer:
		return True
	else:
		return False

def merge_images_v2(img_list, h = None, v = None):
	"""
    Function: merge_images_v2()

    Helper function for analyze_plate(). Given a list of image, stitches the images into one large image.

    Parameters
    ----------
   	img_list : list of Numpy arrays
   		list of images to be stitched together
   	h : list of ints 
   		Stitching parameters for the horizontal axis
   	v : list of ints
		Stitching parameters for the vertical axis
		
    Returns
    -------
    merged : Numpy array
    	merged imaged
    h : list of ints
    	Stitching parameters for the horizontal axis
    v : list of ints
    	Stitching parameters for the vertical axis
  	"""
	if is_square(len(img_list)) is True:
		num = np.int(math.sqrt(len(img_list))+0.5)

		cols = []
		for col in xrange(num):
			imgs_to_merge = []
			for row in xrange(num):
				imgs_to_merge += [img_list[row + col*num]]

			if col % 2 == 1:
				imgs_to_merge.reverse()

			cols += [imgs_to_merge]

		is_h_none = h is None

		if is_h_none is True:
			h = []
			v = []
			for col in cols:
				for row in xrange(num-1):
					h_temp, v_temp = cross_corr(col[row], col[row+1])
					if h_temp == 0:
						h_temp += 1
					if v_temp == 0:
						v_temp += 1
					h += [h_temp]
					v += [v_temp]

		# Merge rows using the offsets
		merged_cols = []
		for j in xrange(num):
			merged_col = cols[j][0]
			h_temp = h[j*(num-1):j*(num-1) + num-1]
			v_temp = v[j*(num-1):j*(num-1) + num-1]

			for row in xrange(num-1):
				merged_col = np.concatenate((cols[j][row+1][:-h_temp[row],:-np.sum(v_temp[0:row+1])], merged_col[:,v_temp[row]:]), axis = 0)
			merged_cols += [merged_col]
		xmins = [col.shape[0] for col in merged_cols]
		ymins = [col.shape[1] for col in merged_cols]
		xmin = min(xmins)
		ymin = min(ymins)

		merged_cols_v2 = [merged_col[0:xmin,0:ymin] for merged_col in merged_cols]
		merged_cols = merged_cols_v2

		if is_h_none is True:
			for j in xrange(num-1):
				if merged_cols[j].shape[0] == 0 or merged_cols[j].shape[1] == 0:
					h_temp = 1
					v_temp = 1
				if merged_cols[j+1].shape[0] == 0 or merged_cols[j+1].shape[1] == 0:
					h_temp = 1
					v_temp = 1
				else:
	 				h_temp, v_temp = cross_corr(merged_cols[j+1], merged_cols[j])
				if h_temp == 0:
					h_temp += 1
				if v_temp == 0:
					v_temp += 1
				h +=[h_temp]
				v +=[v_temp]

		# Merge the merged rows by column together using the offsets
		merged = merged_cols[0]
		h_temp = h[num*(num-1):]
		v_temp = v[num*(num-1):]

		for j in xrange(num-1):
			merged = np.concatenate((merged_cols[j+1][np.sum(h_temp[0:j+1]):,:v_temp[j]],merged[:-h_temp[j],:]), axis = 1)

		return merged, h, v

	else:
		print "Not a square grid!"
