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
from SLIP_functions import plot_slip_joint_plot, fit_kde, compute_p_values, save_classified_wells, compute_p_lysis_posterior
from SLIP_functions import classify_infections_gmm, classify_infections_gmm2
from keio_names import get_keio_names, pos_to_strain
import seaborn as sns
import pandas as pd
import pymc3 as pm
import json

direc = "/scratch/users/nquach/datatxt/"
save_direc = "/scratch/users/nquach/plots/"
pooled_file = open(os.path.join(direc, 'pooled_classified_data.pkl'), 'r')

pooled_dict = pickle.load(pooled_file)

titles = list(pooled_dict.keys())
	
sns.set_context('notebook', font_scale = 1.1)
sns.set_style('white')
sns.set_style('ticks')

sky_blue = (86./255., 180./255., 233./255.)
bluish_green = (0, 158./255., 115./255.)
reddish_purple = (204./255.,121./255.,167./255.)
black = (0.,0.,0.)

fig, axes = plt.subplots(16,12, figsize = (4*12, 4*16))

index = 0

'''
for row in range(16):
	for col in range(12):
		if index < len(titles):
			title = titles[index]
			print title
			lytic_list = zip(*pooled_dict[title][0])
			lyso_list = zip(*pooled_dict[title][1])
			uninfected_list = zip(*pooled_dict[title][2])
			if len(lytic_list) > 0:
				axes[row,col].scatter(lytic_list[0], lytic_list[1], c = bluish_green, alpha = 0.75, edgecolors = 'none')
			if len(lyso_list) > 0:
				axes[row,col].scatter(lyso_list[0], lyso_list[1], c = reddish_purple, alpha = 0.75, edgecolors = 'none')
			if len(uninfected_list) > 0:
				axes[row,col].scatter(uninfected_list[0], uninfected_list[1], c = black, alpha = 0.75, edgecolors = 'none')
			axes[row,col].set_xlabel('Lytic Reporter Intensity (au)')
			axes[row,col].set_ylabel('Lysogenic Reporter Intensity (au)')
			axes[row,col].set_xmargin(0.1)
			axes[row,col].set_ymargin(0.1)
			axes[row,col].autoscale(enable = True, axis = 'both', tight = True)
			axes[row,col].set_title(title)
			sns.despine()
			index = index + 1
'''
for row in range(16):
	for col in range(12):
		if index < len(titles):
			title = titles[index]
			print title
			lytic_list = zip(*pooled_dict[title][0])
			lyso_list = zip(*pooled_dict[title][1])
			uninfected_list = zip(*pooled_dict[title][2])
			N_lytic = len(lytic_list)
			N_lysogenic = len(lyso_list)
			x, posterior = compute_p_lysis_posterior(N_lytic, N_lysogenic)
			axes[row, col].plot(x, posterior, color = 'k', linewidth = 2)
			axes[row, col].set_xlim([0, 1])
			axes[row, col].set_xlabel('Probability of lysis')
			axes[row, col].set_ylabel('Probability density')
			axes[row, col].set_title(title)
			sns.despine()
			index = index + 1


plt.tight_layout()
figure_name = os.path.join(save_direc, 'pooled_classified_hits_lysis_posteriors.pdf')
fig.savefig(figure_name, bbox_inches = 'tight')

pooled_file.close()
