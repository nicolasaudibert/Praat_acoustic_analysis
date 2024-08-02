#####################################################################################################################
# extract_acoustic_parameters_from_selected_labels_v20_fileslist.praat
#
# This script extracts from a set of matched sound files and textgrid objects every label
# matching a set defined in a text file with their start and end time and duration. If no text file is specified,
# all non-empty intervals are processed.
# Aligned labels on other tiers can be optionally extracted (either one selected secondary tier or all other),
# as well as previous and next labels on each target tier.
# The following acoustic parameters can be extracted on request (on a variable number of user-defined points,
# plus mean value on the whole interval, min, max with time and standard deviation on request):
# - F0
# - Intensity
# - Formants 1 to 4
# - Spectral center of gravity (CoG), optionally limited to frequency bands, with energy in each band (may be used for the extraction of the level of energy in frequency bands)
# - Spectral dispersion
# - Spectral tilt (difference of energy between frequency bands)
# - Harmonics-to-noise ratio (HNR)
# - Zero-crossing rate (ZCR)
# - Cepstral Peak Prominence (CPP) + smoothed version CPPS
# - MFCC coefficients
# - Signal value (may be useful when applied to physiological signals)
#
# Selected features and target extraction points are defined in an external text file.
#
# The script assumes that all textgrids have the same structure.
#
# This version of the script selects pairs of .TextGrid and .wav files listed in the input file
# (tab-separated file with columns named wav and TextGrid). It may be useful as an alternative to
# other versions when names of wav and TextGrid files are not paired following a regular enough pattern.
#
# Author: Nicolas Audibert, LPP UMR7018, January 2011 - last modified November 2023
# https://lpp.in2p3.fr/nicolas-audibert/
#####################################################################################################################

form Extract_acoustic_parameters
	comment Folder with textgrids (all textgrids must have the same structure)
	text textgrids_folder extraits_NCCFr
	comment Text file with the list of paired .wav and .TextGrid files to be processed
	sentence paired_files_list liste_paires_wav_textgrids_extraits_NCCFr.txt
	comment Folder with sounds (leave empty if same as textgrids folder or to extract only duration and context)
	text wavefiles_folder 
	comment Output file
	text results_file resultats_analyse_acoustique.txt
	comment Index of the tier with labels to be processed
	positive reference_tier 1
	comment Path to the parameters file
	text parameters_file parametres_extraction_mesures_acoustiques_voix_femme.txt
	comment File with relative positions of the target points for parameters extraction
	text extraction_points_definition_file positions_5points.txt
	comment Text file that contains the labels to be processed (leave empty to process all non-empty labels)
	text dictionary_file voyelles_francais_LIMSI.txt
endform

# Clear info window
clearinfo

# Read the list of textgrids from the specified file
flistTable = Read Table from tab-separated file: paired_files_list$
wavColIndex = Get column index: "wav"
tgColIndex = Get column index: "TextGrid"

# Check that columns "wav" and "TextGrid" both exist
if wavColIndex>0 and tgColIndex>0
	# Call the main script
	runScript: "extract_acoustic_parameters_from_selected_labels_v19.praat", textgrids_folder$, wavefiles_folder$, results_file$, reference_tier, parameters_file$, extraction_points_definition_file$, dictionary_file$
else
	appendInfoLine: "File ", paired_files_list$, " does not have columns named wav and TextGrid. Aborting."
endif
