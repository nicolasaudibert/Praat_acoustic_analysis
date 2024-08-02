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
# The script assumes that matched textgrid and sound files have the same name, and that all textgrids
# have the same structure.
#
# This version of the script selects .TextGrid files listed in the input file.
#
# Author: Nicolas Audibert, LPP UMR7018, January 2011 - last modified November 2023
# https://lpp.cnrs.fr/nicolas-audibert/
#####################################################################################################################

form Extract_acoustic_parameters
	comment Folder with textgrids (all textgrids must have the same structure)
	text textgrids_folder .
	comment Text file with the list of .TextGrid files to be processed
	sentence textgrid_files_list myFileList.txt
	comment Optional suffix and extension of TextGrid files
	# When the script is called from command line, set suffix to * to use an empty string (no suffix)
	# This is because empty strings given as command line arguments will be ignored, shifting the arguments list and generating unexpected behaviour
	text textgrids_suffix 
	text textgrids_extension .TextGrid
	comment Folder with sounds (leave empty if same as textgrids folder or to extract only duration and context)
	text wavefiles_folder
	comment Optional suffix and extension of sound files
	# When the script is called from command line, set suffix to * to use an empty string (no suffix)
	# This is because empty strings given as command line arguments will be ignored, shifting the arguments list and generating unexpected behaviour
	text wavefiles_suffix 
	text wavefiles_extension .wav
	comment Output file
	text results_file acoustic_results.txt
	comment Index of the tier with labels to be processed
	positive reference_tier 1
	comment Path to the parameters file
	text parameters_file extract_durations_F0.txt
	comment File with relative positions of the target points for parameters extraction
	text extraction_points_definition_file positions_3points_10_50_90.txt
	comment Text file that contains the labels to be processed (leave empty to process all non-empty labels)
	text dictionary_file French_vowels_SAMPA.txt
endform

# Clear info window
clearinfo

# Adjust textgrid and/or sound files suffix if needed
if textgrids_extension$ = "*"
	textgrids_extension$ = ""
endif
if wavefiles_suffix$ = "*"
	wavefiles_suffix$ = ""
endif

# Read the list of textgrids from the specified file
flist = Read Strings from raw text file: textgrid_files_list$

# Build a table with paired .TextGrid and .wav file names
nFiles = Get number of strings
flistWav = Replace all: textgrids_suffix$ + textgrids_extension$ + "$", wavefiles_suffix$ + wavefiles_extension$, 1, "regular expressions"
flistTable = Create Table with column names: "filelist", nFiles, "wav TextGrid"
for iFile from 1 to nFiles
	selectObject: flistWav
	wavFile$ = Get string: iFile
	selectObject: flist
	tgFile$ = Get string: iFile
	selectObject: flistTable
	Set string value: iFile, "wav", wavFile$
	Set string value: iFile, "TextGrid", tgFile$
endfor
removeObject: flist, flistWav
selectObject: flistTable

# Call the main script
runScript: "extract_acoustic_parameters_from_selected_labels_v20.praat", textgrids_folder$, wavefiles_folder$, results_file$, reference_tier, parameters_file$, extraction_points_definition_file$, dictionary_file$
