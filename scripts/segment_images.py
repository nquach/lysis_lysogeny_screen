'''
segment_images.py

Script used to analyze raw SLIP images. Segments images using DeepCell algorithm.

Written by Nicolas Quach and David Van Valen
'''
from SLIP_functions import segment_SLIP
import os

#Define root directory path to data folders
root_direc = '/scratch/users/nquach/hits12/'

#define directory path to infection data (all positions)
data_direc = os.path.join(root_direc, 'data')

#define directory path to where you want to save masks
mask_direc = os.path.join(root_direc, 'masks')


segment_SLIP(data_direc, mask_direc, alphabet = ['A','B','C','D','E','F','G','H'], columns = range(1,13), pos_list = range(16))