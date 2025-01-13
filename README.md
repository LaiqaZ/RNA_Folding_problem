# RNA_Folding_problem
----------------------------------------------------------------------Data Preparation-----------------------------------------------------------------------------
***Dataset for Training***
This project uses RNA 3D structures obtained from the RCSB Protein Data Bank (PDB). The data includes 202 RNA structures downloaded based on the following criteria:

    Macromolecule: RNA
    Experimental Method: X-ray diffraction

I copied the accession IDs for the RNA structures in a text file 'pdb_ids.txt'. I used a shell script 'batch_download.sh' provided on RCSB to download the pdb structures. The structures were downloaded into a directory named pdb_files

following command was run to batch download the pdb.gz files.
./batch_download.sh -f pdb_ids.txt -o pdb_files -p

decompressed the .gz files using :
gunzip *.gz

***Dataset for Testing (predicted RNA structures)***
The following dataset was downloaded (random choice) from RNA puzzles. 
https://github.com/RNA-Puzzles/raw_dataset_and_for_assessment/blob/master/raw/PZ13.tar.gz

This data set has predicted structures and those structures were moved to a directory 'rna_structures'.

-----------------------------------------------------------------------------SCRIPTS----------------------------------------------------------------------------


							                                                  	***************************************
						                                                  		* Training script: training_script.py *
						                                                  		***************************************
								
This script calculates interaction scores for RNA base pairs based on the distances between C3' atoms extracted from RNA 3D structures in PDB files. These scores are saved as text files for later use in analysis and scoring.

Input : pdb_files directory containing .pdb files.

Methods used:
*compute_distance: Calculates Euclidean distance between two 3D points.
*parse_pdb: Extracts C3' atom coordinates from a PDB file.
*count_distances: Counts distances between atom pairs, grouped into bins.
*compute_scores: Computes log-ratio scores for base pairs based on distance counts.
*save_scores: Saves computed scores to text files for each base pair.

Output:
A directory 'scores' containing text files for the score of each base pair.

command to run the script: python training_script.py
					                                                			************************************
				                                                				* Plotting: Interaction_profile.py *
					                                                			************************************

Input: a directory 'scores' which was created in the last step.

Methods used:
*read_scores: Reads and returns scores from a specified text file.
*plot_profile: Generates and saves a plot for the interaction profile of a given base pair.
*generate_plots: Reads scores for all base pairs and generates plots, saving them in the plots folder.

Output: the script saves plots as .png files in the plots directory. Each file is named after the base pair (e.g., AA_profile.png).

command to run the script: python interaction_profile.py

							                                                        	***************************
					                                                        			* Scoring: gibbs_score.py *
						                                                        		***************************
								

Encapsulates the logic for loading precomputed scores, parsing RNA structure files, calculating distances, interpolating scores, and computing Gibbs free energy.
Input: a directory 'score', that was generated in previous steps, containing pre-computed score for each base pair.
	a directory 'rna_structures' that contains the pdb files for predicted structures.

Methods:
*load_scores: Loads precomputed scores for base pairs from text files.
*interpolate_score: Estimates a score for a given base pair at a specific distance using linear interpolation.
*parse_pdb: Extracts C3' atom coordinates from a PDB file.
*compute_gibbs_free_energy: Calculates normalized Gibbs free energy for an RNA structure based on valid atom pairs.

Output:
terminal displays all the results and the structure with lowest gibbs free energy. The results are also saved in a file 'gibbs_free_energy_results.txt'.

command to run the script: python gibbs_score.py

*** NOTE: The total Gibbs free energy is divided by the number of valid atom pairs. This ensures that the energy is comparable across structures of varying sizes, providing a per-interaction average rather than an absolute total. Structures with fewer valid pairs will not inherently score lower, making the comparison fair and size-independent. I did this normalization step because the sizes of the predicted structures were very big comparatively.***
