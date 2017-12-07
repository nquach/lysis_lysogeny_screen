# lysis_lysogeny_screen
All data and scripts for the lysis-lysogeny screen of the Keio collection

ORGANIZATION:
************
- Python scripts for analyzing SLIP microscope data is contained in the /scripts folder.
- all scatterplots and posterior distribution plots are contained in the /plots folder.
- JSON .txt files containing cell fate statistics summaries are contained in the /datatxt folder.
- hits_map_name.xlsx contains the map of the final hits arrayed in 96 well format.
- /all_data_pkl folder contains all data from pass1 and pass2 of the Keio screen.

SCRIPTS:
*********
__________________________
Analysis pipeline scripts:
__________________________

analyze_hit_plate.py : perform post-segmentation data analysis on pass1 hits

analyze_keio_plate.py : perform post-segmentation data analysis on a Keio library plate

compile_pass_pkl.py: saves the classified and unclassified data for all plates within a pass of the Keio collection (either pass1 or 	pass2=verification pass)

compile_pkl.py: saves the classified and unclassified data for each plate into a .pkl file.

compile_pooled_hits_dict.py : converts the compiled classified data for the pooled data for pass1 hits into a dict of lists.
Saves the data as a .txt file in JSON format

data_analysis.py : contains all functions for data organization, compilation and visualization. Functions are run in search_hits.py. 

DataVisualizer.ipynb : iPython notebook that allows visualization of all data collected in the lysis-lysogeny screen.

keio_names.py : contains all the functions for mapping positions to a gene name based off the map of the Keio collection

plot_pooled_hits.py : produces the plots for all the final hits, using pooled data between pass1 and pass2

plot.py : produces plots for individual plates of the Keio collection or of pass1 hits.

pool_data.py : merges all classified and unclassified data into large data structures indexed by gene name.

pos_to_name.py : produces a dictionary mapping positions and plates to a gene name

search_hits.py : script that allows you to search the dataset generated in pass1 and the pooled dataset for pass1 hits.

segment_images.py : script that segments all the images for a plate

SLIP_function.py: contains all functions for the analysis pipeline
____________________________
DeepCell dependency scripts:
____________________________

cnn_functions.py : contains all the functions necessary for running conv-nets

model_zoo.py : contains the conv-net architecture definitions

tifffile.py : functions to manipulate .tif files


RUNNING IMAGES THROUGH THE ANALYSIS PIPELINE:
*********************************************
1) Run segment_images.py to segment images for a particular plate. Data must be organized in the directory tree of format root->position->images. Directory into which the masks will go into must be already created at runtime.

2) Run analyze_keio_plate.py or analyze_hits_plate.py to extract fluorescence values from images using the masks.

3) Run plot.py to produce scatterplots of classified and unclassified fluorescence data, and the lysis and MOI posterior distributions for each strain on the plate. This should also save cell fate statistics summary for each strain on the plate in JSON formatted .txt files.


