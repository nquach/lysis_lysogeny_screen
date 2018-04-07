# Keio Collection Lysis Lysogeny Screen
This repository contains data and scripts for the lysis-lysogeny screen of the Keio collection.

## **Organization of the repository**
* Python scripts for analyzing SLIP microscope data is contained in the /scripts folder.
* all scatterplots and posterior distribution plots are contained in the /plots folder.
* JSON .txt files containing cell fate statistics summaries are contained in the /datatxt folder.
* hits_map_name.xlsx contains the map of the final hits arrayed in 96 well format.
* /all_data_pkl folder contains all data from pass1 and pass2 of the Keio screen.

## How to run images through the analysis pipeline:
0. Set the `ROOT_DIREC` variable to the local absolute path of the cloned repo.
1. Run `segment_images.py` to segment images for a particular plate. Data must be organized in the directory tree of format root->position->images. Directory into which the masks will go into must be already created at runtime.

2. Run `analyze_keio_plate.py` or `analyze_hits_plate.py` to extract fluorescence values from images using the masks.

3. Run `plot.py` to produce scatterplots of classified and unclassified fluorescence data, and the lysis and MOI posterior distributions for each strain on the plate. This should also save cell fate statistics summary for each strain on the plate in JSON formatted .txt files.

## **Discription of .pkl files**
* **all_pooled_classified_data.pkl** : contains all fluorescence data collected (i.e. pass 1 and pass 2 data pooled), classified into lytic, lysogenic, and uninfected events. Use the python cPickle package to unzip file. Unzips into a dict of list of lists. 
Key = gene name. Item = list([lytic data, lysogenic data, uninfected data]), where each data list is of the form list([GFP data, Cherry data])
* **all_pooled_data.pkl** : contains all fluorescence data collected (i.e. pass 1 and pass 2 data pooled), unclassified. Use python cPickle package to unzip file. Unzips into a dict of lists. Key = gene name. Item = list([GFP data, Cherry data])
* **pass1_all_data.pkl** : contains fluorescence data from pass 1 of the library screen, unclassified. Use python cPickle to unzip file. Unzips into a dict of lists. Key = gene name. Item = list([GFP data, Cherry data])
* **pass1_classified_data.pkl** : contains fluorescence data from pass 1 of the library screen, classified into lytic, lysogenic and uninfected events. Use the python cPickle package to unzip file. Unzips into a dict of list of lists. 
Key = gene name. Item = list([lytic data, lysogenic data, uninfected data]), where each data list is of the form list([GFP data, Cherry data])
* **pass2_all_data.pkl** : contains fluorescence data from the verification pass of the library screen, unclassified. Use python cPickle to unzip file. Unzips into a dict of lists. Key = gene name. Item = list([GFP data, Cherry data])
* **pass2_classified_data.pkl** : contains fluorescence data from pass 1 of the library screen, classified into lytic, lysogenic and uninfected events. Use the python cPickle package to unzip file. Unzips into a dict of list of lists. 
Key = gene name. Item = list([lytic data, lysogenic data, uninfected data]), where each data list is of the form list([GFP data, Cherry data])
* **pass2_pooled_classified_data** : contains pooled fluorescence data for the strains included in the verification pass (i.e. hits found in pass 1 of the screen), classified nto lytic, lysogenic and uninfected events. Use the python cPickle package to unzip file. Unzips into a dict of list of lists. Key = gene name. Item = list([lytic data, lysogenic data, uninfected data]), where each data list is of the form list([GFP data, Cherry data])
* **K29_mean_cherry** : contains Cherry fluorescence data for the strains in Keio Plate 29. Use the python cPickle package to unzip file. Unzips into a dict of lists. Key = position string (e.g. 'A1'). Item = list of fluorescence values)
* **K29_mean_FITC** : contains FITC fluorescence data for the strains in Keio Plate 29. Use the python cPickle package to unzip file. Unzips into a dict of lists. Key = position string (e.g. 'A1'). Item = list of fluorescence values)
* **/pkl_files/CS_duplicate1** : contains FITC and Cherry fluorescence data for the strains for the carbon source experiments Trial 1. Use the python cPickle package to unzip the files. Unzips into a dict of lists. Key = position string (e.g. 'A1'). Item = list of fluorescence values for either FITC or Cherry channel)
* **/pkl_files/CS_duplicate1** : contains FITC and Cherry fluorescence data for the strains for the carbon source experiments Trial 2. Use the python cPickle package to unzip the files. Unzips into a dict of lists. Key = position string (e.g. 'A1'). Item = list of fluorescence values for either FITC or Cherry channel)

## **Description of scripts**

### iPython Notebooks:
`SearchHits.ipynb` : Interactive notebook to search through the pooled data for genes displaying cell fate statistics meeting user-defined criteria

`CarbonSourcePlotter.ipynb` : Notebook that plots the carbon source experiment data using different x-axis modes.

`HitsHistogramMOI-LysisRatio` : Notebook that generates plots for Figure 3.9 in the manuscript

`HitsScatterplotGenerator` : Notebook that generates plots for Figures 4.1 and 4.2 in the manuscript. Also generates MOI and lysis ratio plots for hits.

`BulkAssayPlotter.ipynb` : Notebook that generates histograms for the bulk assay (Figure 3.5).

`DataVisualizer.ipynb` : iPython notebook that allows visualization of all data collected in the lysis-lysogeny screen.

### Analysis pipeline scripts:

`analyze_hit_plate.py` : perform post-segmentation data analysis on pass1 hits

`analyze_keio_plate.py` : perform post-segmentation data analysis on a Keio library plate

`compile_pass_pkl.py`: saves the classified and unclassified data for all plates within a pass of the Keio collection (either pass1 or 	pass2=verification pass)

`compile_pkl.py`: saves the classified and unclassified data for each plate into a .pkl file.

`compile_pooled_hits_dict.py` : converts the compiled classified data for the pooled data for pass1 hits into a dict of lists.
Saves the data as a .txt file in JSON format

`data_analysis.py` : contains all functions for data organization, compilation and visualization. Functions are run in search_hits.py. 

`keio_names.py` : contains all the functions for mapping positions to a gene name based off the map of the Keio collection

`plot.py` : produces plots for individual plates of the Keio collection or of pass1 hits.

`pool_data.py` : merges all classified and unclassified data into large data structures indexed by gene name.

`segment_images.py` : script that segments all the images for a plate

`SLIP_function.py` : contains all functions for the analysis pipeline

### DeepCell & Analysis Pipeline dependencies:

`cnn_functions.py` : contains all the functions necessary for running conv-nets

`model_zoo.py` : contains the conv-net architecture definitions

`tifffile.py` : functions to manipulate .tif files



