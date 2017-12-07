import matplotlib
matplotlib.use('Agg')
import json
import os
import numpy as np
import xlrd as xls
import matplotlib.pyplot as plt
from keio_names import get_keio_names, pos_to_strain

plate_numbers = [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47, 49, 51, 53, 55, 57, 59, 61, 63, 65, 67, 69, 71, 73, 75, 77, 79, 81, 83, 85, 89, 91, 93, 95]
direc = '/scratch/users/nquach/datatxt/'

keio_names_array = get_keio_names()

name_dict = {}
rows = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
columns = range(1,13)

for plate in plate_numbers:
	for row in rows:
		for col in columns:
			well = row + str(col)
			pos_str = 'K' + str(plate) + well
			print pos_str, pos_to_strain(keio_names_array, plate, well)
			name_dict[pos_str] = pos_to_strain(keio_names_array, plate, well)

file_path = os.path.join(direc, 'pos_to_name.txt')
file = open(file_path, 'w+')
json.dump(name_dict, file)
file.close()

print name_dict

