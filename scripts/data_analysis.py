'''
data_analysis.py

Contains definitions for the functions used to analyze the classified data and identify hits
'''
import matplotlib
matplotlib.use('Agg')
import json
import os
import numpy as np
import xlrd as xls
import matplotlib.pyplot as plt
import utils

ROOT_DIREC = utils.ROOT_DIREC

#List all Keio plate numbers
plate_numbers = ['1_1','1_2', 3, 5, 7, '9_1','9_2', 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 39, 41, 43, 45, 47, 49, 51, 53, 55, 57, 59, 61, 63, 65, 67, 69, 71, 73, 75, 77, 79, 81, 83, 85, 89, 95]
#set root directory
direc = os.path.join(ROOT_DIREC, 'datatxt')

def get_OD_data(plate_numbers):
	"""
    Function: get_OD_data()

    Given plate numbers, produces a dictionary of OD600 plate reader values for each position of each plates, 
    for all plates in the given list of plate numbers, and saves the dict of dicts as a .txt file in JSON format. 
    Outer dict: key = plate number as str, item = inner dict
    Inner dict: key = position, item = OD600 value

    Parameters
    ----------
   	plate_numbers : list 
        list of plate numbers to compile
    
    Returns
    -------
    None : None
		saves the compiled dictionary (dict of dicts) as a .txt file in JSON format
	"""
	#set root direc where platereader data is
	direc = '/scratch/users/nquach/'
	data_direc = os.path.join(direc, 'platereader')
	save_direc = os.path.join(direc, 'datatxt')
	compiled_file_path = os.path.join(save_direc, 'OD_data.txt')
	compiled_file = open(compiled_file_path, 'w+')
	compiled_dict = {}
	letters = ['A','B','C','D','E','F','G','H']
	numbers = range(1,13)
	positions = []
	#compile list of positions
	for letter in letters:
		for number in numbers:
			pos = letter + str(number)
			positions.append(pos)

	for plate_number in plate_numbers:
		print plate_number
		plate_dict = {}
		file_name = 'keio' + str(plate_number) + '.xls'
		file_path = os.path.join(data_direc, file_name)
		wb = xls.open_workbook(file_path)
		sheet = wb.sheet_by_index(0)
		for i in range(0,96):
			pos = positions[i]
			plate_dict[pos] = sheet.cell_value(rowx=i+1, colx=5)

		compiled_dict[plate_number] = plate_dict
	
	json.dump(compiled_dict, compiled_file)
	compiled_file.close()


def compile_dataset(direc, plate_numbers):
	"""
	***DEFUNCT***
    Function: compile_dataset()

    Given a root directory path and a list of plate numbers, produces a dict of dicts of cell fate statistic summaries, indexed by plate and position.
    Outer dict: key = plate number as str, item = inner dict
    Inner dict: key = position, item = list of cell fate statistics

    Parameters
    ----------
   	direc : str
   		str defining the full path to the root directory
   	plate_numbers : list 
        list of plate numbers to compile
    
    Returns
    -------
    None : None
		saves the compiled dictionary as a .txt file in JSON format
	"""
	dataset = {}
	compiled_file_name = 'all_data.txt'
	compiled_file = open(os.path.join(direc, compiled_file_name), "w+")
	for plate_number in plate_numbers:
		print 'Adding Keio ' + str(plate_number)
		file_name = 'keio' + str(plate_number) + '.txt'
		plate_file = open(os.path.join(direc, file_name), 'r')
		plate_dict = json.load(plate_file)
		plate_dict = invert_plate_index(plate_dict)
		dataset[plate_number] = [plate_dict]
		plate_file.close()

	json.dump(dataset, compiled_file)
	compiled_file.close()

def compile_dataset_by_name(direc, plate_numbers):
	"""
	***DEFUNCT***
    Function: compile_dataset()

    Given a root directory path and a list of plate numbers, produces a dict of lists of cell fate statistic summaries, indexed by gene name.
    key = gene name, item = list(list([gene name, number of total cells, number of infected cels, number of uninfected cells, number of lytic cells, number of lysogenic cells, lysis ratio, variance of lysis ratio estimate]))

    Parameters
    ----------
   	direc : str
   		str defining the full path to the root directory
   	plate_numbers : list 
        list of plate numbers to compile
    
    Returns
    -------
    None : None
		saves the compiled dictionary (dict of list) as a .txt file in JSON format
	"""
	dataset = {}
	compiled_file_name = 'data_by_name.txt'
	compiled_file = open(os.path.join(direc, compiled_file_name), "w+")
	for plate_number in plate_numbers:
		print 'Adding Keio ' + str(plate_number)
		file_name = 'keio' + str(plate_number) + '.txt'
		plate_file = open(os.path.join(direc, file_name), 'r')
		plate_dict = json.load(plate_file)
		plate_dict = invert_plate_index(plate_dict)
		for pos in plate_dict.keys():
			stats_list = plate_dict[pos]
			title = stats_list[0]
			total_cells = stats_list[1]
			num_infected = stats_list[2]
			num_uninfected = stats_list[3]
			num_lytic = stats_list[4]
			num_lyso = stats_list[5]
			mode = stats_list[6]
			var = stats_list[7]
			pos_str = 'K' + str(plate_number) + pos
			dataset[title] = [pos_str, total_cells, num_infected, num_uninfected, num_lytic, num_lyso, mode, var]
		plate_file.close()

	json.dump(dataset, compiled_file)
	compiled_file.close()

def load_name_dataset(direc):
	"""
    Function: load_name_dataset()

    Given the path to the directory containing data_by_name.txt, returns the dict of lists of cell fate summaries indexed by gene name. 
    key = gene name, item = list([gene name, number of total cells, number of infected cels, number of uninfected cells, number of lytic cells, number of lysogenic cells, lysis ratio, variance of lysis ratio estimate])

    Parameters
    ----------
   	direc : str
   		str defining the full path to the directory containing data_by_name.txt
    
    Returns
    -------
    compiled_data : dict of lists
		data structure containing the cell fate summaries indexed by gene name.
	"""
	compiled_file_name = 'data_by_name.txt'
	compiled_file = open(os.path.join(direc, compiled_file_name), "r")
	compiled_data = json.load(compiled_file)
	compiled_file.close()
	return compiled_data

def load_maynard_hits(direc):
	"""
    Function: load_maynard_hits()

    Given the path to the directory containing maynard_hits.txt, returns a list of gene names found to be hits in the Maynard screen

    Parameters
    ----------
   	direc : str
   		str defining the full path to the directory containing maynard_hits.txt
    
    Returns
    -------
    data : list of str
		list of gene names found to be essential for phage infection.
	"""
	file_name = 'maynard_hits.txt'
	file = open(os.path.join(direc, file_name), "r")
	data = json.load(file)
	file.close()
	return data


def load_PPIs(direc):
	"""
    Function: load_PPIs()

    Given the path to the directory containing PPI.txt, returns a list of E. coli gene names found to have protein-protein interactions
    with lambda phage proteins in the Blasche screen.

    Parameters
    ----------
   	direc : str
   		str defining the full path to the directory containing PPI.txt
    
    Returns
    -------
    data : list of str
		list of E. coli gene names found to have protein-protein interactions
    	with lambda phage proteins in the Blasche screen.
	"""
	file_name = 'PPI.txt'
	file = open(os.path.join(direc, file_name), "r")
	data = json.load(file)
	file.close()
	return data


def load_transcription_f(direc):
	"""
    Function: load_transcription_f()

    Given the path to the directory containing transcription.txt, returns genes encoding E. coli transcription factors.

    Parameters
    ----------
   	direc : str
   		str defining the full path to the directory containing PPI.txt
    
    Returns
    -------
    data : list of str
		list of E. coli transcription factor genes
    """
	file_name = 'transcription.txt'
	file = open(os.path.join(direc, file_name), "r")
	data = json.load(file)
	file.close()
	return data

def load_dataset(direc):
	"""
    Function: load_dataset()

    Given the path to the directory containing all_data.txt, returns the dict of dicts containing the cell fate summaries indexed by plate and position
    Outer dict: key = plate number as str, item = inner dict
    Inner dict: key = position, item = list([gene name, number of total cells, number of infected cels, number of uninfected cells, number of lytic cells, number of lysogenic cells, lysis ratio, variance of lysis ratio estimate])

    Parameters
    ----------
   	direc : str
   		str defining the full path to the directory containing all_data.txt
    
    Returns
    -------
    compiled_data : dict of dicts
		Dict of dicts containing the cell fate summaries indexed by plate and position
    """
	compiled_file_name = 'all_data.txt'
	compiled_file = open(os.path.join(direc, compiled_file_name), "r")
	compiled_data = json.load(compiled_file)
	compiled_file.close()
	return compiled_data

def load_pooled_hits(direc):
	"""
    Function: load_pooled_hits()

    Given the path to the directory containing pooled_hits_stats.txt, returns the dict of lists containing the cell fate summaries of 
    pooled data of pass1 hits, index by gene name. 
    key = gene name, item = list([gene name, number of total cells, number of infected cels, number of uninfected cells, number of lytic cells, number of lysogenic cells, lysis ratio, variance of lysis ratio estimate])

    
    Parameters
    ----------
   	direc : str
   		str defining the full path to the directory containing pooled_hits_stats.txt
    
    Returns
    -------
    compiled_data : dict of dicts
		Dict of lists containing the cell fate summaries of pooled data of pass1 hits, indexed by gene name
    """
	compiled_file_name = 'pooled_hits_stats.txt'
	compiled_file = open(os.path.join(direc, compiled_file_name), "r")
	compiled_data = json.load(compiled_file)
	compiled_file.close()
	return compiled_data

def load_OD_data(direc):
	"""
    Function: load_OD_data()

    Given the path to the directory containing OD_data.txt, returns dict of dicts for the OD600 platereader data, indexed by plate and position
    Outer dict: key = plate, item = inner dict
    Inner dict: key = position, item = OD600 reading

    Parameters
    ----------
   	direc : str
   		str defining the full path to the directory containing OD_data.txt
    
    Returns
    -------
    data : dict of dicts
		dict of dicts containing OD600 platereader data.
    """
	data_file_name = 'OD_data.txt'
	data_file = open(os.path.join(direc, data_file_name), 'r')
	data = json.load(data_file)
	data_file.close()
	return data

def invert_index(old_index):
	"""
    Function: invert_index()

    Given a position (e.g. 'A1'), returns the inverted column position. (e.g. A1 -> A12, B2 -> B11)

    Parameters
    ----------
   	old_index : str
   		index to be inverted.
    
    Returns
    -------
    new_index : str
		inverted index
    """
	letter = old_index[0]
	old_num = int(old_index[1:])
	new_num = str(13 - old_num)
	new_index = letter + new_num
	return new_index

def invert_plate_index(old_dict):
	"""
    Function: invert_plate_index()

    Given a dictionary indexed by position, returns a dictionary with keys replaced with the inverted index (e.g. A1 -> A12, B2 -> B11)

    Parameters
    ----------
   	old_dict : dict
   		Dict whose keys you wish to invert
    
    Returns
    -------
    new_dict : dict
		dict with inverted position keys.
    """
	new_dict = {}
	for key in old_dict.keys():
		new_key = invert_index(key)
		new_dict[new_key] = old_dict[key]
	return new_dict

def print_data(direc, plate = None):
	"""
    Function: print_data()

    Given the directory containing all_data.txt and OD_data.txt, and a plate number, prints out the cell fate summaries of all strains on the plate

    Parameters
    ----------
   	direc : str
   		path to directory containing all_data.txt and OD_data.txt
   	plate : int or str
   		Plate number you which to print
    
    Returns
    -------
    None : None
		Prints cell fate summaries of all strains on the selected plate.
    """
	dataset = load_dataset(direc)
	OD_data = load_OD_data(direc)
	for plate_num in sorted(dataset.keys()):
		if plate != None and str(plate) != str(plate_num):
			continue
		else:
			plate_dict = dataset[plate_num][0]
			OD_dict = OD_data[plate_num]
			print 'Data for Keio ' + str(plate_num)
			for pos in sorted(plate_dict.keys()):
				stats_list = plate_dict[pos]
				#print stats_list
				OD = OD_dict[pos]
				title = stats_list[0]
				total_cells = stats_list[1]
				num_infected = stats_list[2]
				num_uninfected = stats_list[3]
				num_lytic = stats_list[4]
				mode = stats_list[6]
				var = stats_list[7]
				if num_infected+num_uninfected > 0:
					frac_infected = float(num_infected)/float(num_infected+num_uninfected)
				else:
					frac_infected = 0
				if title != None:
					pos_str = 'K' + str(plate_num) + pos + ': ' + title
					print pos_str + ' mode=' + str(mode) + ' n_lytic=' + str(num_lytic) + ' infected=' + str(num_infected) + ' total cells=' + str(num_infected+num_uninfected) + ' frac_infected=' + str(frac_infected) + ' OD600=' + str(OD)

def filter_hits(name = None, mode_lb = 0.0, mode_ub = 1.0, n_infected_lb = 0, n_infected_ub = 10000, n_cell_lb = 0, n_cell_ub = 10000, frac_infected_lb = 0, frac_infected_ub = 1.0, var_lb = 0, var_ub = 10000, n_lytic_lb = 0, n_lytic_ub = 10000, n_lyso_lb = 0, n_lyso_ub = 10000, OD_lb = 0, OD_ub = 1, in_maynard = False, is_TF = False, is_PPI = False, save_hits = False):
	"""
    Function: filter_hits()

	Searches the pass1 dataset for strains meeting certain criteria specified by parameters. If you search by name, no other criteria is considered.
    
    Parameters
    ----------
   	name : str
   		gene name. leave as None if you aren't searching by name.
   	mode_lb : float
   		lower bound for lysis ratio
   	mode_ub : float
   		upper bound for lysis ratio
   	n_infected_lb : int
   		lower bound for number of infected cells
   	n_infected_ub : int 
   		upper bound for number of infected cells
   	n_cell_lb : int
   		lower bound for total cells 
   	n_cell_ub : int
   		upper bound for total cells
   	frac_infected_lb : float
   		lower bound for fraction infected.
   	frac_infected_ub : float
   		upper bound for fraction infected.
   	var_lb : float
   		lower bound on the variance of the lysis ratio estimate 
   	var_ub : float
   		upper bound on the variance of the lysis ratio estimate
   	n_lytic_lb : int
   		lower bound on the number of lytic cells
   	n_lytic_ub : int
   		upper bound on the number of lytic cells
   	n_lyso_lb : int
   		lower bound on the number of lysogenic cells
   	n_lyso_ub : int
   		upper bound on the number of lysogenic cells
   	OD_lb : float
   		lower bound on the plate reader OD600
   	OD_ub : float
   		upper bound on the plate reader OD600
   	in_maynard : bool
   		True if you want to search only for hits that are also hits in the Maynard screen
   	is_TF : bool
   		True if you want to search only for hits that are also transcription factors
   	is_PPI : bool
   		True if you want to search only for hits that also have PPIs with lambda proteins (identified by Blasche screen)
   	save_hits : bool
   		True if you want to save the results of the search as a .txt file in JSON format

    Returns
    -------
    None : None
		Prints cell fate summaries of all strains fulfilling all the search criteria.
    """
	print ' '
	print 'Searching through Pass 1 dataset...'
	#keeps track of how many strains meet search criteria
	counter = 0
	#stores the strains that meet search criteria
	hits_list = []
	#search by name
	if name != None:
		#load dataset
		dataset = load_name_dataset(direc)
		stats_list = dataset[name]
		counter += 1
		pos_str = stats_list[0]
		total_cells = stats_list[1]
		num_infected = stats_list[2]
		num_uninfected = stats_list[3]
		num_lytic = stats_list[4]
		mode = stats_list[6]
		var = stats_list[7]
		#calculate fraction infected
		if num_infected+num_uninfected > 0:
			frac_infected = float(num_infected)/float(total_cells)
		else:
			frac_infected = 0
		#print cell fate statistic summary
		print pos_str + ' ' + name + ' mode=' + str(mode) + ' n_lytic=' + str(num_lytic) + ' infected=' + str(num_infected) + ' total cells=' + str(num_infected+num_uninfected) + ' frac_infected=' + str(frac_infected) 
	#search by cell fate statistic
	else:
		#load datasets
		dataset = load_dataset(direc)
		OD_data = load_OD_data(direc)
		#search through all elements in the datasets
		for plate_num in sorted(dataset.keys()):
			plate_dict = dataset[plate_num][0]
			OD_dict = OD_data[plate_num]
			for pos in sorted(plate_dict.keys()):
				stats_list = plate_dict[pos]
				OD = OD_dict[pos]
				title = stats_list[0]
				total_cells = stats_list[1]
				num_infected = stats_list[2]
				num_uninfected = stats_list[3]
				num_lytic = stats_list[4]
				num_lyso = stats_list[5]
				mode = stats_list[6]
				var = stats_list[7]
				if num_infected+num_uninfected > 0:
					frac_infected = float(num_infected)/float(total_cells)
				else:
					frac_infected = 0

				if title != None:
					#check if search criteria are met
					maynard_hits = load_maynard_hits(direc)
					transcription_f = load_transcription_f(direc)
					PPI = load_PPIs(direc)
					pos_str = 'K' + str(plate_num) + pos + ': ' + title
					mode_criteria = (mode >= mode_lb) and (mode <= mode_ub)
					n_infected_criteria = (num_infected >= n_infected_lb) and (num_infected <= n_infected_ub)
					n_cell_criteria = (total_cells >= n_cell_lb) and (total_cells <= n_cell_ub)
					frac_infected_criteria = (frac_infected >= frac_infected_lb) and (frac_infected <= frac_infected_ub)
					var_criteria = (var >= var_lb) and (var <= var_ub)
					lytic_criteria = (num_lytic >= n_lytic_lb) and (num_lytic <= n_lytic_ub)
					lyso_criteria = (num_lyso >= n_lyso_lb) and (num_lyso <= n_lyso_ub)
					OD_criteria = (OD >= OD_lb) and (OD <= OD_ub)
					maynard_criteria = True
					TF_criteria = True
					PPI_criteria = True
					if in_maynard:
						maynard_criteria = title in maynard_hits
					if is_TF:
						TF_criteria = title.lower() in (name.lower() for name in transcription_f)
					if is_PPI:
						PPI_criteria = title.lower() in (name.lower() for name in PPI)
					criteria_met = mode_criteria and n_infected_criteria and n_cell_criteria and frac_infected_criteria and var_criteria and lytic_criteria and lyso_criteria and maynard_criteria and TF_criteria and PPI_criteria
					if criteria_met:
						counter += 1
						hits_list.append(title)
						#uncomment if you want to store the positions of the strains instead of the gene names
						#hits_list.append('K' + str(plate_num) + pos)
						print pos_str + ' mode=' + str(mode) + ' n_lytic=' + str(num_lytic) + ' infected=' + str(num_infected) + ' total cells=' + str(num_infected+num_uninfected) + ' frac_infected=' + str(frac_infected) + ' OD600=' + str(OD)	
	print ' '
	print 'Criteria:' 
	print 'mode=[', mode_lb, mode_ub, ']' 
	print 'n_infected=[', n_infected_lb, n_infected_ub, ']' 
	print 'total_cell=[', n_cell_lb, n_cell_ub, ']'
	print 'frac_infected=[', frac_infected_lb, frac_infected_ub, ']'
	print 'var=[', var_lb, var_ub, ']' 
	print 'n_lytic=[', n_lytic_lb, n_lytic_ub, ']'
	print 'n_lyso=[', n_lyso_lb, n_lyso_ub, ']'
	print 'OD_criteria=[', OD_lb, OD_ub, ']'
	print 'Maynard criteria', in_maynard
	print 'TF criteria', is_TF
	print 'PPI criteria', is_PPI
	print 'Number of hits meeting criteria: ' + str(counter)
	#save the name or position of strains meeting search criteria
	if save_hits:
		save_path = os.path.join(direc,'hit_list.txt')
		save_file = open(save_path, 'w+')
		json.dump(hits_list, save_file)
		save_file.close()

def filter_pooled_hits(name = None, mode_lb = 0.0, mode_ub = 1.0, n_infected_lb = 0, n_infected_ub = 10000, n_cell_lb = 0, n_cell_ub = 10000, frac_infected_lb = 0, frac_infected_ub = 1.0, var_lb = 0, var_ub = 10000, n_lytic_lb = 0, n_lytic_ub = 10000, n_lyso_lb = 0, n_lyso_ub = 10000, in_maynard = False, is_TF = False, is_PPI = False, save_hits = False):
	"""
    Function: filter_pooled_hits()

	Searches the pooled dataset for pass1 hits for strains meeting certain criteria specified by parameters. If you search by name, no other criteria is considered.
    
    Parameters
    ----------
   	name : str
   		gene name. leave as None if you aren't searching by name.
   	mode_lb : float
   		lower bound for lysis ratio
   	mode_ub : float
   		upper bound for lysis ratio
   	n_infected_lb : int
   		lower bound for number of infected cells
   	n_infected_ub : int 
   		upper bound for number of infected cells
   	n_cell_lb : int
   		lower bound for total cells 
   	n_cell_ub : int
   		upper bound for total cells
   	frac_infected_lb : float
   		lower bound for fraction infected.
   	frac_infected_ub : float
   		upper bound for fraction infected.
   	var_lb : float
   		lower bound on the variance of the lysis ratio estimate 
   	var_ub : float
   		upper bound on the variance of the lysis ratio estimate
   	n_lytic_lb : int
   		lower bound on the number of lytic cells
   	n_lytic_ub : int
   		upper bound on the number of lytic cells
   	n_lyso_lb : int
   		lower bound on the number of lysogenic cells
   	n_lyso_ub : int
   		upper bound on the number of lysogenic cells
   	in_maynard : bool
   		True if you want to search only for hits that are also hits in the Maynard screen
   	is_TF : bool
   		True if you want to search only for hits that are also transcription factors
   	is_PPI : bool
   		True if you want to search only for hits that also have PPIs with lambda proteins (identified by Blasche screen)
   	save_hits : bool
   		True if you want to save the results of the search as a .txt file in JSON format

    Returns
    -------
    None : None
		Prints cell fate summaries of all strains fulfilling all the search criteria.
    """
	print ' '
	print 'Searching through pooled hits...'
	counter = 0
	hits_list = []
	#search by name
	if name != None:
		#load dataset
		dataset = load_pooled_hits(direc)
		stats_list = dataset[name]
		counter += 1
		total_cells = stats_list[0]
		num_infected = stats_list[1]
		num_uninfected = stats_list[2]
		num_lytic = stats_list[3]
		mode = stats_list[5]
		var = stats_list[6]
		if num_infected+num_uninfected > 0:
			frac_infected = float(num_infected)/float(total_cells)
		else:
			frac_infected = 0
		print name + ' mode=' + str(mode) + ' n_lytic=' + str(num_lytic) + ' infected=' + str(num_infected) + ' total cells=' + str(num_infected+num_uninfected) + ' frac_infected=' + str(frac_infected) 
	#search by cell fate statistic
	else:
		#load dataset
		dataset = load_pooled_hits(direc)
		for name in sorted(dataset.keys()):
			stats_list = dataset[name]
			title = name
			total_cells = stats_list[0]
			num_infected = stats_list[1]
			num_uninfected = stats_list[2]
			num_lytic = stats_list[3]
			num_lyso = stats_list[4]
			mode = stats_list[5]
			var = stats_list[6]
			if num_infected+num_uninfected > 0:
				frac_infected = float(num_infected)/float(total_cells)
			else:
				frac_infected = 0
			if title != None:
				#Check if search criteria is met
				maynard_hits = load_maynard_hits(direc)
				transcription_f = load_transcription_f(direc)
				PPI = load_PPIs(direc)
				mode_criteria = (mode >= mode_lb) and (mode <= mode_ub)
				n_infected_criteria = (num_infected >= n_infected_lb) and (num_infected <= n_infected_ub)
				n_cell_criteria = (total_cells >= n_cell_lb) and (total_cells <= n_cell_ub)
				frac_infected_criteria = (frac_infected >= frac_infected_lb) and (frac_infected <= frac_infected_ub)
				var_criteria = (var >= var_lb) and (var <= var_ub)
				lytic_criteria = (num_lytic >= n_lytic_lb) and (num_lytic <= n_lytic_ub)
				lyso_criteria = (num_lyso >= n_lyso_lb) and (num_lyso <= n_lyso_ub)
				maynard_criteria = True
				TF_criteria = True
				PPI_criteria = True
				if in_maynard:
					maynard_criteria = title in maynard_hits
				if is_TF:
					TF_criteria = title.lower() in (name.lower() for name in transcription_f)
				if is_PPI:
					PPI_criteria = title.lower() in (name.lower() for name in PPI)
				criteria_met = mode_criteria and n_infected_criteria and n_cell_criteria and frac_infected_criteria and var_criteria and lytic_criteria and lyso_criteria and maynard_criteria and TF_criteria and PPI_criteria
				if criteria_met:
					counter += 1
					hits_list.append(title)
					#uncomment if you want to store the positions of the strains instead of the gene names
					#hits_list.append('K' + str(plate_num) + pos)
					print title + ' mode=' + str(mode) + ' n_lytic=' + str(num_lytic) + ' infected=' + str(num_infected) + ' total cells=' + str(num_infected+num_uninfected) + ' frac_infected=' + str(frac_infected)	
	print ' '
	print 'Criteria:' 
	print 'mode=[', mode_lb, mode_ub, ']' 
	print 'n_infected=[', n_infected_lb, n_infected_ub, ']' 
	print 'total_cell=[', n_cell_lb, n_cell_ub, ']'
	print 'frac_infected=[', frac_infected_lb, frac_infected_ub, ']'
	print 'var=[', var_lb, var_ub, ']' 
	print 'n_lytic=[', n_lytic_lb, n_lytic_ub, ']'
	print 'n_lyso=[', n_lyso_lb, n_lyso_ub, ']'
	print 'Maynard criteria', in_maynard
	print 'TF criteria', is_TF
	print 'PPI criteria', is_PPI
	print 'Number of hits meeting criteria: ' + str(counter)
	#Save hits 
	if save_hits:
		save_path = os.path.join(direc,'hit_list.txt')
		save_file = open(save_path, 'w+')
		json.dump(hits_list, save_file)
		save_file.close()


