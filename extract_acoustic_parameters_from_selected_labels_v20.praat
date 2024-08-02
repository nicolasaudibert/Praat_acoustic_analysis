#####################################################################################################################
# extract_acoustic_parameters_from_selected_labels_v20.praat
#
# This script extracts from a set of matched sound files and textgrid objects every label
# matching a set defined in a text file with their start and end time and duration. If no text file is specified,
# all non-empty intervals are processed.
# Aligned labels on other tiers can be optionally extracted (either one selected secondary tier or all other),
# as well as previous and next labels on each target tier.
# The following acoustic features can be extracted on request (on a variable number of user-defined points,
# plus mean value on the whole interval, min, max with time and standard deviation on request):
# - F0
# - Intensity
# - Formants 1 to 4
# - Spectral center of gravity (CoG), optionally limited to frequency bands, with energy in each band (may be used for the extraction of the level of energy in frequency bands)
# - Spectral dispersion in the same frequency band as CoG
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
# This version of the script assumes that a Table object with the list of WAV and TextGrid files to be processed is selected.
# It is typically called by either extract_acoustic_parameters_from_selected_labels_v20_regex.praat, extract_acoustic_parameters_from_selected_labels_v20_fileslist.praat or extract_acoustic_parameters_from_selected_labels_v20_paired_files_list.praat.
#
# Author: Nicolas Audibert, LPP UMR7018, January 2011 - last modified November 2023
# https://lpp.cnrs.fr/nicolas-audibert/
#####################################################################################################################

form Extract_acoustic_parameters
	comment Folder with textgrids (all textgrids must have the same structure)
	text textgrids_folder 
	comment Folder with sounds (leave empty if same as textgrids folder or to extract only duration and context)
	text wavefiles_folder
	comment Output file
	text results_file 
	comment Index of the tier with labels to be processed
	positive reference_tier 
	comment Path to the parameters file
	text parameters_file extract_all_parameters_default_settings.txt
	comment File with relative positions of the target points for parameters extraction
	text extraction_points_definition_file 
	comment Text file that contains the labels to be processed (leave empty to process all non-empty labels)
	text dictionary_file 
endform

# Clear info window
clearinfo

# Get the list of textgrids (selected Strings object)
flist = selected("Table")
ntextgrids = Get number of rows

# Check values of parameters defined in form
if wavefiles_folder$ = ""
	wavefiles_folder$ = textgrids_folder$
endif
if dictionary_file$ = ""
	filter_labels = 0
else
	filter_labels = 1
	# Get list of segments to be extracted from reference tier
	stringsDictionary = Read Strings from raw text file: dictionary_file$
	dictionarySize = Get number of strings
endif

# Read the external parameters file
parametersTable = Read Table from tab-separated file: parameters_file$
nParameters = Get number of rows
for iParam from 1 to nParameters
	currentParamName$ = Get value: iParam, "Parameter"
	currentParamType$ = Get value: iParam, "Type"
	if currentParamType$ = "Txt"
		currentParamName$ = currentParamName$ + "$"
		'currentParamName$' = Get value: iParam, "Value"
	elsif currentParamType$ = "Num"
		'currentParamName$' = Get value: iParam, "Value"
	else
		appendInfoLine: "File ",parameters_file$ ,", line ", iParam+1," - unknown type for parameter ", currentParamName$
	endif
endfor
removeObject: parametersTable

# Get the names and relative position of the points in a Matrix object
extracted_points_relative_times = Read Matrix from raw text file: extraction_points_definition_file$
n_extracted_points = Get number of rows

# Write the results file header
#fileappend 'results_file$' textgrid_file'tab$'label'tab$'previousLabel'tab$'followingLabel'tab$'start'tab$'end'tab$'duration(s)'tab$'mean_f0(Hz)'tab$'f0_point1(Hz)'tab$'f0_point2(Hz)'tab$'f0_point3(Hz)'newline$'
writeFile: results_file$, "textgrid_file"
if export_sound_files_name
	appendFile: results_file$, tab$, "sound_file"
endif
appendFile: results_file$, tab$, "label", tab$, "start_time", tab$,"end_time", tab$, "duration(s)"
if extract_left_and_right_context
	appendFile: results_file$, tab$, "previousLabel", tab$, "followingLabel"
endif
if extract_F0
	for ipoint from 1 to n_extracted_points
		selectObject: extracted_points_relative_times
		current_point_relative_postition = Get value in cell: ipoint, 1
		appendFile: results_file$, tab$, "F0_pt"+fixed$(ipoint,0)+"_"+fixed$(round(100*current_point_relative_postition),0)+"%(Hz)"
	endfor
	if get_mean
		appendFile: results_file$, tab$, "mean_F0(Hz)"
	endif
	if get_median
		appendFile: results_file$, tab$, "median_F0(Hz)"
	endif
	if get_standard_deviation
		appendFile: results_file$, tab$, "std_dev_F0(Hz)"
	endif
	if get_min_max_with_time
		appendFile: results_file$, tab$, "min_F0(Hz)", tab$, "min_F0_relative_time", tab$, "max_F0(Hz)", tab$, "max_F0_relative_time"
	endif
endif
if extract_intensity
	for ipoint from 1 to n_extracted_points
		selectObject: extracted_points_relative_times
		current_point_relative_postition = Get value in cell: ipoint, 1
		appendFile: results_file$, tab$, "intensity_pt"+fixed$(ipoint,0)+"_"+fixed$(round(100*current_point_relative_postition),0)+"%(dB)"
	endfor
	if get_mean
		appendFile: results_file$, tab$, "mean_intensity(dB)"
	endif
	if get_median
		appendFile: results_file$, tab$, "median_intensity(dB)"
	endif
	if get_standard_deviation
		appendFile: results_file$, tab$, "std_dev_intensity(dB)"
	endif
	if get_min_max_with_time
		appendFile: results_file$, tab$, "min_intensity(dB)", tab$, "min_intensity_relative_time", tab$, "max_intensity(dB)", tab$, "max_intensity_relative_time"
	endif
endif
if (extract_CoG + extract_spectral_dispersion > 0) and (frequencyBandsDefinitionForCoGcomputationFile$<>"none")
	extract_energy_frequency_bands = 1
	frequencyBandsDefinitionForCoGcomputationTable = Read Table from tab-separated file: frequencyBandsDefinitionForCoGcomputationFile$
	nFrequencyBandsForCoGcomputation = Get number of rows
	for iband from 1 to nFrequencyBandsForCoGcomputation
		selectObject: frequencyBandsDefinitionForCoGcomputationTable
		currentBandStartFreqHz = Get value: iband, "minFreq"
		currentBandEndFreqHz = Get value: iband, "maxFreq"
		if get_mean
			appendFile: results_file$, tab$, "mean_intensity_band_"+fixed$(currentBandStartFreqHz,0)+"_"+fixed$(currentBandEndFreqHz,0)+"Hz(dB)"
		endif
		if get_median
			appendFile: results_file$, tab$, "median_intensity_band_"+fixed$(currentBandStartFreqHz,0)+"_"+fixed$(currentBandEndFreqHz,0)+"Hz(dB)"
		endif
		if get_standard_deviation
			appendFile: results_file$, tab$, "std_dev_intensity_band_"+fixed$(currentBandStartFreqHz,0)+"_"+fixed$(currentBandEndFreqHz,0)+"Hz(dB)"
		endif
		if get_min_max_with_time
			appendFile: results_file$, tab$, "min_intensity_band_"+fixed$(currentBandStartFreqHz,0)+"_"+fixed$(currentBandEndFreqHz,0)+"Hz(dB)"
			appendFile: results_file$, tab$, "min_intensity_band_"+fixed$(currentBandStartFreqHz,0)+"_"+fixed$(currentBandEndFreqHz,0)+"Hz_relative_time"
			appendFile: results_file$, tab$, "max_intensity_band_"+fixed$(currentBandStartFreqHz,0)+"_"+fixed$(currentBandEndFreqHz,0)+"Hz(dB)"
			appendFile: results_file$, tab$, "max_intensity_band_"+fixed$(currentBandStartFreqHz,0)+"_"+fixed$(currentBandEndFreqHz,0)+"Hz_relative_time"
		endif
	endfor
	for ipoint from 1 to n_extracted_points
		selectObject: extracted_points_relative_times
		current_point_relative_postition = Get value in cell: ipoint, 1
		for iband from 1 to nFrequencyBandsForCoGcomputation
			selectObject: frequencyBandsDefinitionForCoGcomputationTable
			currentBandStartFreqHz = Get value: iband, "minFreq"
			currentBandEndFreqHz = Get value: iband, "maxFreq"
			appendFile: results_file$, tab$, "intensity_band_"+fixed$(currentBandStartFreqHz,0)+"_"+fixed$(currentBandEndFreqHz,0)+"Hz_pt"+fixed$(ipoint,0)+"_"+fixed$(round(100*current_point_relative_postition),0)+"%(dB)"
		endfor
	endfor
	# sound pressure-related constant needed to convert energy values in frequency bands (in Pa^2.s-1) to intensity values in dB
	p0 = 2e-5
	p0squared = p0*p0
else
	extract_energy_frequency_bands = 0
endif
if extract_CoG
	for ipoint from 1 to n_extracted_points
		selectObject: extracted_points_relative_times
		current_point_relative_postition = Get value in cell: ipoint, 1
		appendFile: results_file$, tab$, "CoG_pt"+fixed$(ipoint,0)+"_"+fixed$(round(100*current_point_relative_postition),0)+"%(Hz)"
		if extract_energy_frequency_bands
			for iband from 1 to nFrequencyBandsForCoGcomputation
				selectObject: frequencyBandsDefinitionForCoGcomputationTable
				currentBandStartFreqHz = Get value: iband, "minFreq"
				currentBandEndFreqHz = Get value: iband, "maxFreq"
				appendFile: results_file$, tab$, "CoG_band_"+fixed$(currentBandStartFreqHz,0)+"_"+fixed$(currentBandEndFreqHz,0)+"Hz_pt"+fixed$(ipoint,0)+"_"+fixed$(round(100*current_point_relative_postition),0)+"%(Hz)"
			endfor
		endif
	endfor
	if get_mean
		appendFile: results_file$, tab$, "mean_CoG(Hz)"
		if extract_energy_frequency_bands
			for iband from 1 to nFrequencyBandsForCoGcomputation
				selectObject: frequencyBandsDefinitionForCoGcomputationTable
				currentBandStartFreqHz = Get value: iband, "minFreq"
				currentBandEndFreqHz = Get value: iband, "maxFreq"
				appendFile: results_file$, tab$, "mean_CoG_band_"+fixed$(currentBandStartFreqHz,0)+"_"+fixed$(currentBandEndFreqHz,0)+"Hz(Hz)"
			endfor
		endif
	endif
	if get_median
		appendFile: results_file$, tab$, "median_CoG(Hz)"
		if extract_energy_frequency_bands
			for iband from 1 to nFrequencyBandsForCoGcomputation
				selectObject: frequencyBandsDefinitionForCoGcomputationTable
				currentBandStartFreqHz = Get value: iband, "minFreq"
				currentBandEndFreqHz = Get value: iband, "maxFreq"
				appendFile: results_file$, tab$, "median_CoG_band_"+fixed$(currentBandStartFreqHz,0)+"_"+fixed$(currentBandEndFreqHz,0)+"Hz(Hz)"
			endfor
		endif
	endif
	if get_standard_deviation
		appendFile: results_file$, tab$, "std_dev_CoG(Hz)"
		if extract_energy_frequency_bands
			for iband from 1 to nFrequencyBandsForCoGcomputation
				selectObject: frequencyBandsDefinitionForCoGcomputationTable
				currentBandStartFreqHz = Get value: iband, "minFreq"
				currentBandEndFreqHz = Get value: iband, "maxFreq"
				appendFile: results_file$, tab$, "std_dev_CoG_band_"+fixed$(currentBandStartFreqHz,0)+"_"+fixed$(currentBandEndFreqHz,0)+"Hz(Hz)"
			endfor
		endif
	endif
	if get_min_max_with_time
		appendFile: results_file$, tab$, "min_CoG(Hz)", tab$, "min_CoG_relative_time", tab$, "max_CoG(Hz)", tab$, "max_CoG_relative_time"
		if extract_energy_frequency_bands
			for iband from 1 to nFrequencyBandsForCoGcomputation
				selectObject: frequencyBandsDefinitionForCoGcomputationTable
				currentBandStartFreqHz = Get value: iband, "minFreq"
				currentBandEndFreqHz = Get value: iband, "maxFreq"
				appendFile: results_file$, tab$, "min_CoG_band_"+fixed$(currentBandStartFreqHz,0)+"_"+fixed$(currentBandEndFreqHz,0)+"Hz(Hz)"
				appendFile: results_file$, tab$, "min_CoG_band_"+fixed$(currentBandStartFreqHz,0)+"_"+fixed$(currentBandEndFreqHz,0)+"Hz_relative_time"
				appendFile: results_file$, tab$, "max_CoG_band_"+fixed$(currentBandStartFreqHz,0)+"_"+fixed$(currentBandEndFreqHz,0)+"Hz(Hz)"
				appendFile: results_file$, tab$, "max_CoG_band_"+fixed$(currentBandStartFreqHz,0)+"_"+fixed$(currentBandEndFreqHz,0)+"Hz_relative_time"
			endfor
		endif
	endif
endif
if extract_spectral_dispersion
	for ipoint from 1 to n_extracted_points
		selectObject: extracted_points_relative_times
		current_point_relative_postition = Get value in cell: ipoint, 1
		appendFile: results_file$, tab$, "spectral_dispersion_pt"+fixed$(ipoint,0)+"_"+fixed$(round(100*current_point_relative_postition),0)+"%(Hz)"
		if extract_energy_frequency_bands
			for iband from 1 to nFrequencyBandsForCoGcomputation
				selectObject: frequencyBandsDefinitionForCoGcomputationTable
				currentBandStartFreqHz = Get value: iband, "minFreq"
				currentBandEndFreqHz = Get value: iband, "maxFreq"
				appendFile: results_file$, tab$, "spectral_dispersion_band_"+fixed$(currentBandStartFreqHz,0)+"_"+fixed$(currentBandEndFreqHz,0)+"Hz_pt"+fixed$(ipoint,0)+"_"+fixed$(round(100*current_point_relative_postition),0)+"%(Hz)"
			endfor
		endif
	endfor
	if get_mean
		appendFile: results_file$, tab$, "mean_spectral_dispersion(Hz)"
		if extract_energy_frequency_bands
			for iband from 1 to nFrequencyBandsForCoGcomputation
				selectObject: frequencyBandsDefinitionForCoGcomputationTable
				currentBandStartFreqHz = Get value: iband, "minFreq"
				currentBandEndFreqHz = Get value: iband, "maxFreq"
				appendFile: results_file$, tab$, "mean_spectral_dispersion_band_"+fixed$(currentBandStartFreqHz,0)+"_"+fixed$(currentBandEndFreqHz,0)+"Hz(Hz)"
			endfor
		endif
	endif
	if get_median
		appendFile: results_file$, tab$, "median_spectral_dispersion(Hz)"
		if extract_energy_frequency_bands
			for iband from 1 to nFrequencyBandsForCoGcomputation
				selectObject: frequencyBandsDefinitionForCoGcomputationTable
				currentBandStartFreqHz = Get value: iband, "minFreq"
				currentBandEndFreqHz = Get value: iband, "maxFreq"
				appendFile: results_file$, tab$, "median_spectral_dispersion_band_"+fixed$(currentBandStartFreqHz,0)+"_"+fixed$(currentBandEndFreqHz,0)+"Hz(Hz)"
			endfor
		endif
	endif
	if get_standard_deviation
		appendFile: results_file$, tab$, "std_dev_spectral_dispersion(Hz)"
		if extract_energy_frequency_bands
			for iband from 1 to nFrequencyBandsForCoGcomputation
				selectObject: frequencyBandsDefinitionForCoGcomputationTable
				currentBandStartFreqHz = Get value: iband, "minFreq"
				currentBandEndFreqHz = Get value: iband, "maxFreq"
				appendFile: results_file$, tab$, "std_dev_spectral_dispersion_band_"+fixed$(currentBandStartFreqHz,0)+"_"+fixed$(currentBandEndFreqHz,0)+"Hz(Hz)"
			endfor
		endif
	endif
	if get_min_max_with_time
		appendFile: results_file$, tab$, "min_spectral_dispersion(Hz)", tab$, "min_spectral_dispersion_relative_time", tab$, "max_spectral_dispersion(Hz)", tab$, "max_spectral_dispersion_relative_time"
		if extract_energy_frequency_bands
			for iband from 1 to nFrequencyBandsForCoGcomputation
				selectObject: frequencyBandsDefinitionForCoGcomputationTable
				currentBandStartFreqHz = Get value: iband, "minFreq"
				currentBandEndFreqHz = Get value: iband, "maxFreq"
				appendFile: results_file$, tab$, "min_spectral_dispersion_band_"+fixed$(currentBandStartFreqHz,0)+"_"+fixed$(currentBandEndFreqHz,0)+"Hz(Hz)"
				appendFile: results_file$, tab$, "min_spectral_dispersion_band_"+fixed$(currentBandStartFreqHz,0)+"_"+fixed$(currentBandEndFreqHz,0)+"Hz_relative_time"
				appendFile: results_file$, tab$, "max_spectral_dispersion_band_"+fixed$(currentBandStartFreqHz,0)+"_"+fixed$(currentBandEndFreqHz,0)+"Hz(Hz)"
				appendFile: results_file$, tab$, "max_spectral_dispersion_band_"+fixed$(currentBandStartFreqHz,0)+"_"+fixed$(currentBandEndFreqHz,0)+"Hz_relative_time"
			endfor
		endif
	endif
endif
if extract_spectral_tilt
	for ipoint from 1 to n_extracted_points
		selectObject: extracted_points_relative_times
		current_point_relative_postition = Get value in cell: ipoint, 1
		appendFile: results_file$, tab$, "spectral_tilt_pt"+fixed$(ipoint,0)+"_"+fixed$(round(100*current_point_relative_postition),0)+"%(dB)"
	endfor
	if get_mean
		appendFile: results_file$, tab$, "mean_spectral_tilt(dB)"
	endif
	if get_median
		appendFile: results_file$, tab$, "median_spectral_tilt(dB)"
	endif
	if get_standard_deviation
		appendFile: results_file$, tab$, "std_dev_spectral_tilt(dB)"
	endif
	if get_min_max_with_time
		appendFile: results_file$, tab$, "min_spectral_tilt(dB)", tab$, "min_spectral_tilt_relative_time", tab$, "max_spectral_tilt(dB)", tab$, "max_spectral_tilt_relative_time"
	endif
endif
if extract_HNR
	for ipoint from 1 to n_extracted_points
		selectObject: extracted_points_relative_times
		current_point_relative_postition = Get value in cell: ipoint, 1
		appendFile: results_file$, tab$, "HNR_pt"+fixed$(ipoint,0)+"_"+fixed$(round(100*current_point_relative_postition),0)+"%(dB)"
	endfor
	if get_mean
		appendFile: results_file$, tab$, "mean_HNR(dB)"
	endif
	if get_median
		appendFile: results_file$, tab$, "median_HNR(dB)"
	endif
	if get_standard_deviation
		appendFile: results_file$, tab$, "std_dev_HNR(dB)"
	endif
	if get_min_max_with_time
		appendFile: results_file$, tab$, "min_HNR(dB)", tab$, "min_HNR_relative_time", tab$, "max_HNR(dB)", tab$, "max_HNR_relative_time"
	endif
endif
if extract_ZCR
	for ipoint from 1 to n_extracted_points
		selectObject: extracted_points_relative_times
		current_point_relative_postition = Get value in cell: ipoint, 1
		appendFile: results_file$, tab$, "ZCR_pt"+fixed$(ipoint,0)+"_"+fixed$(round(100*current_point_relative_postition),0)+"%"
	endfor
	if get_mean
		appendFile: results_file$, tab$, "mean_ZCR"
	endif
	if get_median
		appendFile: results_file$, tab$, "median_ZCR"
	endif
	if get_standard_deviation
		appendFile: results_file$, tab$, "std_dev_ZCR"
	endif
	if get_min_max_with_time
		appendFile: results_file$, tab$, "min_ZCR", tab$, "min_ZCR_relative_time", tab$, "max_ZCR", tab$, "max_ZCR_relative_time"
	endif
endif
if extract_formants
	for ipoint from 1 to n_extracted_points
		selectObject: extracted_points_relative_times
		current_point_relative_postition = Get value in cell: ipoint, 1
		appendFile: results_file$, tab$, "F1_pt"+fixed$(ipoint,0)+"_"+fixed$(round(100*current_point_relative_postition),0)+"%(Hz)"
	endfor
	if get_mean
		appendFile: results_file$, tab$, "mean_F1(Hz)"
	endif
	if get_median
		appendFile: results_file$, tab$, "median_F1(Hz)"
	endif
	if get_standard_deviation
		appendFile: results_file$, tab$, "std_dev_F1(Hz)"
	endif
	if get_min_max_with_time
		appendFile: results_file$, tab$, "min_F1(Hz)", tab$, "min_F1_relative_time", tab$, "max_F1(Hz)", tab$, "max_F1_relative_time"
	endif
	for ipoint from 1 to n_extracted_points
		selectObject: extracted_points_relative_times
		current_point_relative_postition = Get value in cell: ipoint, 1
		appendFile: results_file$, tab$, "F2_pt"+fixed$(ipoint,0)+"_"+fixed$(round(100*current_point_relative_postition),0)+"%(Hz)"
	endfor
	if get_mean
		appendFile: results_file$, tab$, "mean_F2(Hz)"
	endif
	if get_median
		appendFile: results_file$, tab$, "median_F2(Hz)"
	endif
	if get_standard_deviation
		appendFile: results_file$, tab$, "std_dev_F2(Hz)"
	endif
	if get_min_max_with_time
		appendFile: results_file$, tab$, "min_F2(Hz)", tab$, "min_F2_relative_time", tab$, "max_F2(Hz)", tab$, "max_F2_relative_time"
	endif
	for ipoint from 1 to n_extracted_points
		selectObject: extracted_points_relative_times
		current_point_relative_postition = Get value in cell: ipoint, 1
		appendFile: results_file$, tab$, "F3_pt"+fixed$(ipoint,0)+"_"+fixed$(round(100*current_point_relative_postition),0)+"%(Hz)"
	endfor
	if get_mean
		appendFile: results_file$, tab$, "mean_F3(Hz)"
	endif
	if get_median
		appendFile: results_file$, tab$, "median_F3(Hz)"
	endif
	if get_standard_deviation
		appendFile: results_file$, tab$, "std_dev_F3(Hz)"
	endif
	if get_min_max_with_time
		appendFile: results_file$, tab$, "min_F3(Hz)", tab$, "min_F3_relative_time", tab$, "max_F3(Hz)", tab$, "max_F3_relative_time"
	endif
	for ipoint from 1 to n_extracted_points
		selectObject: extracted_points_relative_times
		current_point_relative_postition = Get value in cell: ipoint, 1
		appendFile: results_file$, tab$, "F4_pt"+fixed$(ipoint,0)+"_"+fixed$(round(100*current_point_relative_postition),0)+"%(Hz)"
	endfor
	if get_mean
		appendFile: results_file$, tab$, "mean_F4(Hz)"
	endif
	if get_median
		appendFile: results_file$, tab$, "median_F4(Hz)"
	endif
	if get_standard_deviation
		appendFile: results_file$, tab$, "std_dev_F4(Hz)"
	endif
	if get_min_max_with_time
		appendFile: results_file$, tab$, "min_F4(Hz)", tab$, "min_F4_relative_time", tab$, "max_F4(Hz)", tab$, "max_F4_relative_time"
	endif
endif
if extract_MFCC
	for ipoint from 1 to n_extracted_points
		for iMFCC from 0 to numberOfMFCCcoefficients
			selectObject: extracted_points_relative_times
			current_point_relative_postition = Get value in cell: ipoint, 1
			appendFile: results_file$, tab$, "MFCC"+fixed$(iMFCC,0)+"_pt"+fixed$(ipoint,0)+"_"+fixed$(round(100*current_point_relative_postition),0)+"%"
		endfor
	endfor
	if get_mean
		for iMFCC from 0 to numberOfMFCCcoefficients
			appendFile: results_file$, tab$, "mean_MFCC"+fixed$(iMFCC,0)
		endfor
	endif
	if get_median
		for iMFCC from 0 to numberOfMFCCcoefficients
			appendFile: results_file$, tab$, "median_MFCC"+fixed$(iMFCC,0)
		endfor
	endif
	if get_standard_deviation
		for iMFCC from 0 to numberOfMFCCcoefficients
			appendFile: results_file$, tab$, "std_dev_MFCC"+fixed$(iMFCC,0)
		endfor
	endif
	if get_min_max_with_time
		for iMFCC from 0 to numberOfMFCCcoefficients
			appendFile: results_file$, tab$, "min_MFCC"+fixed$(iMFCC,0)
			appendFile: results_file$, tab$, "min_MFCC"+fixed$(iMFCC,0)+"_relative_time"
			appendFile: results_file$, tab$, "max_MFCC"+fixed$(iMFCC,0)
			appendFile: results_file$, tab$, "max_MFCC"+fixed$(iMFCC,0)+"_relative_time"
		endfor
	endif
endif
if extract_CPP
	# CPP
	for ipoint from 1 to n_extracted_points
		selectObject: extracted_points_relative_times
		current_point_relative_postition = Get value in cell: ipoint, 1
		appendFile: results_file$, tab$, "CPP_pt"+fixed$(ipoint,0)+"_"+fixed$(round(100*current_point_relative_postition),0)+"%(Hz)"
	endfor
	if get_mean
		appendFile: results_file$, tab$, "mean_CPP(dB)"
	endif
	if get_median
		appendFile: results_file$, tab$, "median_CPP(dB)"
	endif
	if get_standard_deviation
		appendFile: results_file$, tab$, "std_dev_CPP(dB)"
	endif
	if get_min_max_with_time
		appendFile: results_file$, tab$, "min_CPP(dB)", tab$, "min_CPP_relative_time", tab$, "max_CPP(dB)", tab$, "max_CPP_relative_time"
	endif
	# CPPS
	for ipoint from 1 to n_extracted_points
		selectObject: extracted_points_relative_times
		current_point_relative_postition = Get value in cell: ipoint, 1
		appendFile: results_file$, tab$, "CPPS_pt"+fixed$(ipoint,0)+"_"+fixed$(round(100*current_point_relative_postition),0)+"%(Hz)"
	endfor
	if get_mean
		appendFile: results_file$, tab$, "mean_CPPS(dB)"
	endif
	if get_median
		appendFile: results_file$, tab$, "median_CPPS(dB)"
	endif
	if get_standard_deviation
		appendFile: results_file$, tab$, "std_dev_CPPS(dB)"
	endif
	if get_min_max_with_time
		appendFile: results_file$, tab$, "min_CPPS(dB)", tab$, "min_CPPS_relative_time", tab$, "max_CPPS(dB)", tab$, "max_CPPS_relative_time"
	endif
endif
if extract_signal_value
	for ipoint from 1 to n_extracted_points
		selectObject: extracted_points_relative_times
		current_point_relative_postition = Get value in cell: ipoint, 1
		appendFile: results_file$, tab$, "signal_value_pt"+fixed$(ipoint,0)+"_"+fixed$(round(100*current_point_relative_postition),0)+"%"
	endfor
	if get_mean
		appendFile: results_file$, tab$, "mean_signal_value"
	endif
	if get_median
		appendFile: results_file$, tab$, "median_signal_value"
	endif
	if get_standard_deviation
		appendFile: results_file$, tab$, "std_dev_signal_value"
	endif
	if get_min_max_with_time
		appendFile: results_file$, tab$, "min_signal_value", tab$, "min_signal_value_relative_time", tab$, "max_signal_value", tab$, "max_signal_value_relative_time"
	endif
endif

# The rest of the results file header will be written only when processing the first textgrid
header_written = 0
# Init ntiers to 0 before actual number of tiers is known
ntiers = 0

# fileappend "'results_file$'" 'newline$'

# Loop every selected textgrid
for itextgrid to ntextgrids
	# Get its name, display it in 'info' windows and read it
	selectObject: flist
	tg$ = Get value: itextgrid, "TextGrid"
	appendInfoLine: "Processing file ", tg$, "..."
	current_tg = Read from file: "'textgrids_folder$'/'tg$'"
	ntiers = Get number of tiers

	if header_written = 0
		# Finish writing results file header on the first loop increment
		selectObject: current_tg
		if secondary_tier>0
			# Get the name of the selected secondary tier
			selectObject: current_tg
			tiername$ = Get tier name: secondary_tier
			# If the name of the secondary tier is empty or includes only spaces, replace it with "tier" followed by its index
			if tiername$ = "" or index_regex(tiername$, "^\s*$")>0
				tiername$ = "tier" + fixed$(secondary_tier,0)
			endif
			appendFile: results_file$, tab$, tiername$, tab$, "'tiername$'_start_time", tab$, "'tiername$'_end_time"
			if extract_left_and_right_context
				appendFile: results_file$, tab$, "previous_'tiername$'", tab$, "next_'tiername$'"
			endif
		elsif secondary_tier = -1
			# Get the names of all interval tiers
			for itier from 1 to ntiers
				# Ignore it if it's the reference tier (already processed) or a point tier (no labels to extract)
				selectObject: current_tg
				interv_tier = Is interval tier: itier
				if itier<>reference_tier and interv_tier=1
					# Get tier name and write it to results file
					selectObject: current_tg
					tiername$ = Get tier name: itier
					# If the name of the tier is empty or includes only spaces, replace it with "tier" followed by its index
					if tiername$ = "" or index_regex(tiername$, "^\s*$")>0
						tiername$ = "tier" + fixed$(itier,0)
					endif
					appendFile: results_file$, tab$, tiername$, tab$, "'tiername$'_start_time", tab$, "'tiername$'_end_time"
					if extract_left_and_right_context
						appendFile: results_file$, tab$, "previous_'tiername$'", tab$, "next_'tiername$'"
					endif
				endif
			endfor
		endif
		# Append a linebreak to results file to finish writing the header
		appendFile: results_file$, newline$
		header_written = 1
	endif

	# Read corresponding sound if at least one type of acoustic analysis is selected
	if extract_F0+extract_intensity+extract_formants+extract_CoG+extract_spectral_dispersion+extract_spectral_tilt+extract_HNR+extract_ZCR+extract_MFCC+extract_CPP+extract_signal_value>0
		selectObject: flist
		snd$ = Get value: itextgrid, "wav"
		current_snd = Open long sound file: "'wavefiles_folder$'/'snd$'"
		current_snd_total_duration = Get total duration
	endif

	# Extract info from every matching interval
	selectObject: current_tg
	ninterv = Get number of intervals: reference_tier
	# Loop every interval on reference tier
	for iinterv from 1 to ninterv
		selectObject: current_tg
		label$ = Get label of interval: reference_tier, iinterv
		# Do something only if the interval label is not empty and matches the set of symbols to be processed (if defined)
		if length(label$)>0
			# Search and replace in label if requested
			if searchAndReplaceInLabels
				label$ = replace_regex$(label$, regexToReplaceInLabels$, replacementStringInLabels$, 0)
			endif
			# Check if the current label is included in the dictionary
			foundString = 0
			if filter_labels
				if dictionarySize > 0
					currentStringIndex = 1
					while foundString = 0 and currentStringIndex <= dictionarySize
						selectObject: stringsDictionary
						currentString$ = Get string: currentStringIndex
						if label$ = currentString$
							foundString = 1
						endif
						currentStringIndex = currentStringIndex + 1
					endwhile
				endif
			endif
			# If filtering is on, extract values only if the current phoneme is included
			if filter_labels = 0 or foundString = 1
				selectObject: current_tg
				#  Extract phonemic context
				if extract_left_and_right_context
					if iinterv>1
						previousLabel$ = Get label of interval: reference_tier, iinterv-1
						# Search and replace in label if requested
						if searchAndReplaceInLabels
							previousLabel$ = replace_regex$(previousLabel$, regexToReplaceInLabels$, replacementStringInLabels$, 0)
						endif
					else
						previousLabel$ = "--undefined--"
					endif
					if iinterv<ninterv
						followingLabel$ = Get label of interval: reference_tier, iinterv+1
						# Search and replace in label if requested
						if searchAndReplaceInLabels
							followingLabel$ = replace_regex$(followingLabel$, regexToReplaceInLabels$, replacementStringInLabels$, 0)
						endif
					else
						followingLabel$ = "--undefined--"
					endif
				endif
				# Extract start and end times, and calculate segment duration
				start_time = Get start time of interval: reference_tier, iinterv
				end_time = Get end time of interval: reference_tier, iinterv
				duration = end_time-start_time

				# Get the signal extract as a Sound object, and compute acoustic parameters
				if extract_F0+extract_intensity+extract_formants+extract_CoG+extract_spectral_dispersion+extract_spectral_tilt+extract_HNR+extract_ZCR+extract_MFCC+extract_signal_value>0
					# Get the start and end time of the signal extract ( Sound object including offset before and after the target interval)
					extract_start_time = start_time - offset_for_acoustic_parameters_extraction_milliseconds/1000
					extract_end_time = end_time + offset_for_acoustic_parameters_extraction_milliseconds/1000
					# Check that the extract start and end times are not off limits
					if extract_start_time < 0
						extract_start_time = 0
					endif
					if extract_end_time > current_snd_total_duration
						extract_end_time = current_snd_total_duration
					endif
					selectObject: current_snd
					current_snd_extract = Extract part: extract_start_time, extract_end_time, "yes"
					nChannels = Get number of channels
					if nChannels > 1
						current_snd_extract_tmp = current_snd_extract
						# If special value 0 is used as target channel, mix channels
						if target_channel=0
							current_snd_extract = Convert to mono
						else
						# Otherwise extract the target channel
							current_snd_extract = Extract one channel: target_channel
						endif
						removeObject: current_snd_extract_tmp
					endif

					# Get F0, intensity and formants values of the extract if requested
					# Store values extracted for each target point in TableOfReal objects
					if extract_F0
						selectObject: current_snd_extract
						current_pitch = noprogress To Pitch: timeStepF0detection, minF0, maxF0
						f0_values_extracted_points = Create TableOfReal: "f0_values", n_extracted_points, 1
						for ipoint from 1 to n_extracted_points
							selectObject: extracted_points_relative_times
							current_point_relative_position = Get value in cell: ipoint, 1
							current_point_time = start_time + (current_point_relative_position * duration)
							selectObject: current_pitch
							current_value = Get value at time: current_point_time, "Hertz", "Linear"
							selectObject: f0_values_extracted_points
							Set value: ipoint, 1, current_value
						endfor
						selectObject: current_pitch
						mean_f0 = Get mean: start_time, end_time, "Hertz"
						if get_mean
							mean_f0 = Get mean: start_time, end_time, "Hertz"
						endif
						if get_median
							median_f0 = Get quantile: start_time, end_time, 0.5, "Hertz"
						endif
						if get_standard_deviation
							std_f0 = Get standard deviation: start_time, end_time, "Hertz"
						endif
						if get_min_max_with_time
							min_f0 = Get minimum: start_time, end_time, "Hertz", "Parabolic"
							min_f0_relative_time = Get time of minimum: start_time, end_time, "Hertz", "Parabolic"
							min_f0_relative_time = (min_f0_relative_time - start_time) / duration
							max_f0 = Get maximum: start_time, end_time, "Hertz", "Parabolic"
							max_f0_relative_time = Get time of maximum: start_time, end_time, "Hertz", "Parabolic"
							max_f0_relative_time = (max_f0_relative_time - start_time) / duration
						endif
						removeObject: current_pitch
					endif
					
					if extract_intensity
						selectObject: current_snd_extract
						current_intensity = noprogress To Intensity: minF0, timeStepIntensityExtraction, subtractMeanFromIntensityValues$
						intensity_values_extracted_points = Create TableOfReal: "intensity_values", n_extracted_points, 1
						for ipoint from 1 to n_extracted_points
							selectObject: extracted_points_relative_times
							current_point_relative_position = Get value in cell: ipoint, 1
							current_point_time = start_time + (current_point_relative_position * duration)
							selectObject: current_intensity
							current_value = Get value at time: current_point_time, "Cubic"
							selectObject: intensity_values_extracted_points
							Set value: ipoint, 1, current_value
						endfor
						selectObject: current_intensity
						if get_mean
							mean_intensity = Get mean: start_time, end_time, "dB"
						endif
						if get_median
							median_intensity = Get quantile: start_time, end_time, 0.5
						endif
						if get_standard_deviation
							std_intensity = Get standard deviation: start_time, end_time
						endif
						if get_min_max_with_time
							min_intensity = Get minimum: start_time, end_time, "Cubic"
							min_intensity_relative_time = Get time of minimum: start_time, end_time, "Cubic"
							min_intensity_relative_time = (min_intensity_relative_time - start_time) / duration
							max_intensity = Get maximum: start_time, end_time, "Cubic"
							max_intensity_relative_time = Get time of maximum: start_time, end_time, "Cubic"
							max_intensity_relative_time = (max_intensity_relative_time - start_time) / duration
						endif
						removeObject: current_intensity
					endif
					
					if extract_CoG + extract_spectral_dispersion + extract_spectral_tilt > 0
						if extract_energy_frequency_bands
							for iband from 1 to nFrequencyBandsForCoGcomputation
								energy_extracted_points_band'iband' = Create TableOfReal: "energy_values_band'iband'", n_extracted_points, 1
							endfor
						endif
						if extract_CoG
							cog_values_extracted_points = Create TableOfReal: "cog_values", n_extracted_points, 1
							if extract_energy_frequency_bands
								for iband from 1 to nFrequencyBandsForCoGcomputation
									cog_values_extracted_points_band'iband' = Create TableOfReal: "cog_values_band'iband'", n_extracted_points, 1
								endfor
							endif
						endif
						if extract_spectral_dispersion
							spectral_dispersion_values_extracted_points = Create TableOfReal: "spectral_dispersion_values", n_extracted_points, 1
							if extract_energy_frequency_bands
								for iband from 1 to nFrequencyBandsForCoGcomputation
									spectral_dispersion_values_extracted_points_band'iband' = Create TableOfReal: "spectral_dispersion_values_band'iband'", n_extracted_points, 1
								endfor
							endif
						endif
						if extract_spectral_tilt
							spectral_tilt_values_extracted_points = Create TableOfReal: "spectral_tilt_values", n_extracted_points, 1
						endif
						for ipoint from 1 to n_extracted_points
							selectObject: extracted_points_relative_times
							current_point_relative_position = Get value in cell: ipoint, 1
							current_point_time = start_time + (current_point_relative_position * duration)
							# Get number of points in extract
							current_extract_start_time = current_point_time - (extractsDurationForCoGcomputationMilliseconds/2000)
							if current_extract_start_time < start_time - (offset_for_acoustic_parameters_extraction_milliseconds/1000)
								current_extract_start_time = start_time - (offset_for_acoustic_parameters_extraction_milliseconds/1000)
							endif
							current_extract_end_time = current_point_time + (extractsDurationForCoGcomputationMilliseconds/2000)
							if current_extract_end_time > end_time + (offset_for_acoustic_parameters_extraction_milliseconds/1000)
								current_extract_end_time = end_time - (offset_for_acoustic_parameters_extraction_milliseconds/1000)
							endif
							currentSoundSliceDuration = current_extract_end_time - current_extract_start_time
							selectObject: current_snd_extract
							currentSoundSlice = Extract part: current_extract_start_time, current_extract_end_time, windowShapeForCenterOfGravityComputation$, relativeWidthForCenterOfGravityComputation, "no"
							currentSoundSliceSpectrum = noprogress To Spectrum: useFFTinCenterOfGravityComputation$
							# Apply cepstral smoothing if requested
							if cepstralSmoothingBandwidthForCenterOfGravityComputationHz > 0
								selectObject: currentSoundSliceSpectrum
								currentSoundSliceSpectrumSmoothed = Cepstral smoothing: cepstralSmoothingBandwidthForCenterOfGravityComputationHz
								removeObject: currentSoundSliceSpectrum
								currentSoundSliceSpectrum = currentSoundSliceSpectrumSmoothed
							endif
							if extract_energy_frequency_bands
								for iband from 1 to nFrequencyBandsForCoGcomputation
									selectObject: frequencyBandsDefinitionForCoGcomputationTable
									currentBandStartFreqHz = Get value: iband, "minFreq"
									currentBandEndFreqHz = Get value: iband, "maxFreq"
									selectObject: currentSoundSliceSpectrum
									currentSoundSliceSpectrum_band'iband' = Copy: "soundSliceSpectrum_band'iband'"
									Filter (pass Hann band): currentBandStartFreqHz, currentBandEndFreqHz, spectralSmoothingForCoGcomputationOnBandsHertz
									selectObject: currentSoundSliceSpectrum
									current_energy_value_Pa2s = Get band energy: currentBandStartFreqHz, currentBandEndFreqHz
									current_energy_value_dB = 10*log10(current_energy_value_Pa2s/(currentSoundSliceDuration*p0squared))
									selectObject: energy_extracted_points_band'iband'
									Set value: ipoint, 1, current_energy_value_dB
								endfor
							endif
							if extract_CoG
								selectObject: currentSoundSliceSpectrum
								current_CoG_value = Get centre of gravity: powerValueInCenterOfGravityComputation
								selectObject: cog_values_extracted_points
								Set value: ipoint, 1, current_CoG_value
								if extract_energy_frequency_bands
									for iband from 1 to nFrequencyBandsForCoGcomputation
										selectObject: currentSoundSliceSpectrum_band'iband'
										current_CoG_value_band'iband' = Get centre of gravity: powerValueInCenterOfGravityComputation
										selectObject: cog_values_extracted_points_band'iband'
										Set value: ipoint, 1, current_CoG_value_band'iband'
									endfor
								endif
							endif
							if extract_spectral_dispersion
								selectObject: currentSoundSliceSpectrum
								current_spectral_dispersion_value = Get standard deviation: powerValueInCenterOfGravityComputation
								selectObject: spectral_dispersion_values_extracted_points
								Set value: ipoint, 1, current_spectral_dispersion_value
								if extract_energy_frequency_bands
									for iband from 1 to nFrequencyBandsForCoGcomputation
										selectObject: currentSoundSliceSpectrum_band'iband'
										current_spectral_dispersion_value_band'iband' = Get standard deviation: powerValueInCenterOfGravityComputation
										selectObject: spectral_dispersion_values_extracted_points_band'iband'
										Set value: ipoint, 1, current_spectral_dispersion_value_band'iband'
									endfor
								endif
							endif
							if extract_spectral_tilt
								selectObject: currentSoundSliceSpectrum
								current_spectral_tilt_value = Get band energy difference: spectralTiltFirstBandLowFrequency, spectralTiltFirstBandHighFrequency, spectralTiltSecondBandLowFrequency, spectralTiltSecondBandHighFrequency
								selectObject: spectral_tilt_values_extracted_points
								Set value: ipoint, 1, current_spectral_tilt_value
							endif
							removeObject: currentSoundSlice, currentSoundSliceSpectrum
							if extract_energy_frequency_bands
								for iband from 1 to nFrequencyBandsForCoGcomputation
									removeObject: currentSoundSliceSpectrum_band'iband'
								endfor
							endif
						endfor
						
						if get_mean
							selectObject: current_snd_extract
							currentSoundSlice = Extract part: start_time, end_time, windowShapeForCenterOfGravityComputation$, relativeWidthForCenterOfGravityComputation, "no"
							currentSoundSliceDuration = end_time - start_time
							currentSoundSliceSpectrum = noprogress To Spectrum: useFFTinCenterOfGravityComputation$
							# Apply cepstral smoothing if requested
							if cepstralSmoothingBandwidthForCenterOfGravityComputationHz > 0
								selectObject: currentSoundSliceSpectrum
								currentSoundSliceSpectrumSmoothed = Cepstral smoothing: cepstralSmoothingBandwidthForCenterOfGravityComputationHz
								removeObject: currentSoundSliceSpectrum
								currentSoundSliceSpectrum = currentSoundSliceSpectrumSmoothed
							endif
							# mean energy in each band
							if extract_energy_frequency_bands
								for iband from 1 to nFrequencyBandsForCoGcomputation
									selectObject: frequencyBandsDefinitionForCoGcomputationTable
									currentBandStartFreqHz = Get value: iband, "minFreq"
									currentBandEndFreqHz = Get value: iband, "maxFreq"
									selectObject: currentSoundSliceSpectrum
									currentSoundSliceSpectrum_band'iband' = Copy: "soundSliceSpectrum_band'iband'"
									Filter (pass Hann band): currentBandStartFreqHz, currentBandEndFreqHz, spectralSmoothingForCoGcomputationOnBandsHertz
									selectObject: currentSoundSliceSpectrum
									current_energy_value_Pa2s = Get band energy: currentBandStartFreqHz, currentBandEndFreqHz
									mean_energy_band'iband' = 10*log10(current_energy_value_Pa2s/(currentSoundSliceDuration*p0squared))
								endfor
							endif
							if extract_CoG
								selectObject: currentSoundSliceSpectrum
								mean_CoG = Get centre of gravity: powerValueInCenterOfGravityComputation
								if extract_energy_frequency_bands
									for iband from 1 to nFrequencyBandsForCoGcomputation
										selectObject: currentSoundSliceSpectrum_band'iband'
										mean_CoG_band'iband' = Get centre of gravity: powerValueInCenterOfGravityComputation
									endfor
								endif
							endif
							if extract_spectral_dispersion
								selectObject: currentSoundSliceSpectrum
								mean_spectral_dispersion = Get standard deviation: powerValueInCenterOfGravityComputation
								if extract_energy_frequency_bands
									for iband from 1 to nFrequencyBandsForCoGcomputation
										selectObject: currentSoundSliceSpectrum_band'iband'
										mean_spectral_dispersion_band'iband' = Get standard deviation: powerValueInCenterOfGravityComputation
									endfor
								endif
							endif
							if extract_spectral_tilt
								selectObject: currentSoundSliceSpectrum
								mean_spectral_tilt = Get band energy difference: spectralTiltFirstBandLowFrequency, spectralTiltFirstBandHighFrequency, spectralTiltSecondBandLowFrequency, spectralTiltSecondBandHighFrequency
							endif
							removeObject: currentSoundSlice, currentSoundSliceSpectrum
							if extract_energy_frequency_bands
								for iband from 1 to nFrequencyBandsForCoGcomputation
									removeObject: currentSoundSliceSpectrum_band'iband'
								endfor
							endif
						endif
						if extract_energy_frequency_bands
							for iband from 1 to nFrequencyBandsForCoGcomputation
								energy_values_sliding_window_band'iband' = Create Table with column names: "table_energy_values_band'iband'", 0, "times values"
							endfor
						endif
						if extract_CoG
							cog_values_sliding_window = Create Table with column names: "table_CoG_values", 0, "times values"
							if extract_energy_frequency_bands
								for iband from 1 to nFrequencyBandsForCoGcomputation
									cog_values_sliding_window_band'iband' = Create Table with column names: "table_CoG_values_band'iband'", 0, "times values"
								endfor
							endif
						endif
						if extract_spectral_dispersion
							spectral_dispersion_values_sliding_window = Create Table with column names: "table_spectral_dispersion_values", 0, "times values"
							if extract_energy_frequency_bands
								for iband from 1 to nFrequencyBandsForCoGcomputation
									spectral_dispersion_values_sliding_window_band'iband' = Create Table with column names: "table_spectral_dispersion_values_band'iband'", 0, "times values"
								endfor
							endif
						endif
						if extract_spectral_tilt
							spectral_tilt_values_sliding_window = Create Table with column names: "table_spectral_tilt_values", 0, "times values"
						endif
						if get_median+get_standard_deviation+get_min_max_with_time > 0
							# If the computation of the median, standard deviation or minimum/maximum is requested, compute values using a sliding window and store them in a Table object
							# current_point_index = 0
							current_point_index_CoG = 0
							current_point_index_spectral_dispersion = 0
							current_point_index_spectral_tilt = 0
							if extract_energy_frequency_bands
								for iband from 1 to nFrequencyBandsForCoGcomputation
									current_point_index_CoG_band'iband' = 0
									current_point_index_spectral_dispersion_band'iband' = 0
									current_point_index_energy_band'iband' = 0
								endfor
							endif
							current_point_time = start_time
							while current_point_time <= end_time
								current_extract_start_time = current_point_time - (extractsDurationForCoGcomputationMilliseconds/2000)
								if current_extract_start_time < start_time - (offset_for_acoustic_parameters_extraction_milliseconds/1000)
									current_extract_start_time = start_time - (offset_for_acoustic_parameters_extraction_milliseconds/1000)
								endif
								current_extract_end_time = current_point_time + (extractsDurationForCoGcomputationMilliseconds/2000)
								if current_extract_end_time > end_time + (offset_for_acoustic_parameters_extraction_milliseconds/1000)
									current_extract_end_time = end_time - (offset_for_acoustic_parameters_extraction_milliseconds/1000)
								endif
								selectObject: current_snd_extract
								currentSoundSlice = Extract part: current_extract_start_time, current_extract_end_time, windowShapeForCenterOfGravityComputation$, relativeWidthForCenterOfGravityComputation, "no"
								currentSoundSliceDuration = current_extract_end_time - current_extract_start_time
								currentSoundSliceSpectrum = noprogress To Spectrum: useFFTinCenterOfGravityComputation$
								# Apply cepstral smoothing if requested
								if cepstralSmoothingBandwidthForCenterOfGravityComputationHz > 0
									selectObject: currentSoundSliceSpectrum
									currentSoundSliceSpectrumSmoothed = Cepstral smoothing: cepstralSmoothingBandwidthForCenterOfGravityComputationHz
									removeObject: currentSoundSliceSpectrum
									currentSoundSliceSpectrum = currentSoundSliceSpectrumSmoothed
								endif
								if extract_energy_frequency_bands
									for iband from 1 to nFrequencyBandsForCoGcomputation
										selectObject: frequencyBandsDefinitionForCoGcomputationTable
										currentBandStartFreqHz = Get value: iband, "minFreq"
										currentBandEndFreqHz = Get value: iband, "maxFreq"
										# filtering to get CoG and/or spectral dispersion in each band
										selectObject: currentSoundSliceSpectrum
										currentSoundSliceSpectrum_band'iband' = Copy: "soundSliceSpectrum_band'iband'"
										Filter (pass Hann band): currentBandStartFreqHz, currentBandEndFreqHz, spectralSmoothingForCoGcomputationOnBandsHertz
										# energy in each band
										selectObject: currentSoundSliceSpectrum
										current_energy_value_Pa2s = Get band energy: currentBandStartFreqHz, currentBandEndFreqHz
										current_energy_value_dB = 10*log10(current_energy_value_Pa2s/(currentSoundSliceDuration*p0squared))
										selectObject: energy_values_sliding_window_band'iband'
										Append row
										current_point_index_energy_band'iband' = current_point_index_energy_band'iband' + 1
										Set numeric value: current_point_index_energy_band'iband', "times", current_point_time
										Set numeric value: current_point_index_energy_band'iband', "values", current_energy_value_dB
									endfor
								endif
								if extract_CoG
									selectObject: currentSoundSliceSpectrum
									current_CoG_value = Get centre of gravity: powerValueInCenterOfGravityComputation
									if current_CoG_value<>undefined
										selectObject: cog_values_sliding_window
										Append row
										current_point_index_CoG = current_point_index_CoG + 1
										Set numeric value: current_point_index_CoG, "times", current_point_time
										Set numeric value: current_point_index_CoG, "values", current_CoG_value
									endif
									for iband from 1 to nFrequencyBandsForCoGcomputation
										selectObject: currentSoundSliceSpectrum_band'iband'
										current_CoG_value_band'iband' = Get centre of gravity: powerValueInCenterOfGravityComputation
										if current_CoG_value_band'iband'<>undefined
											selectObject: cog_values_sliding_window_band'iband'
											Append row
											current_point_index_CoG_band'iband' = current_point_index_CoG_band'iband' + 1
											Set numeric value: current_point_index_CoG_band'iband', "times", current_point_time
											Set numeric value: current_point_index_CoG_band'iband', "values", current_CoG_value_band'iband'
										endif
									endfor
								endif
								if extract_spectral_dispersion
									selectObject: currentSoundSliceSpectrum
									current_spectral_dispersion_value = Get standard deviation: powerValueInCenterOfGravityComputation
									if current_spectral_dispersion_value<>undefined
										selectObject: spectral_dispersion_values_sliding_window
										Append row
										current_point_index_spectral_dispersion = current_point_index_spectral_dispersion + 1
										Set numeric value: current_point_index_spectral_dispersion, "times", current_point_time
										Set numeric value: current_point_index_spectral_dispersion, "values", current_spectral_dispersion_value
									endif
									for iband from 1 to nFrequencyBandsForCoGcomputation
										selectObject: currentSoundSliceSpectrum_band'iband'
										current_spectral_dispersion_value_band'iband' = Get standard deviation: powerValueInCenterOfGravityComputation
										if current_spectral_dispersion_value_band'iband'<>undefined
											selectObject: spectral_dispersion_values_sliding_window_band'iband'
											Append row
											current_point_index_spectral_dispersion_band'iband' = current_point_index_spectral_dispersion_band'iband' + 1
											Set numeric value: current_point_index_spectral_dispersion_band'iband', "times", current_point_time
											Set numeric value: current_point_index_spectral_dispersion_band'iband', "values", current_spectral_dispersion_value_band'iband'
										endif
									endfor
								endif
								if extract_spectral_tilt
									selectObject: currentSoundSliceSpectrum
									current_spectral_tilt_value = Get band energy difference: spectralTiltFirstBandLowFrequency, spectralTiltFirstBandHighFrequency, spectralTiltSecondBandLowFrequency, spectralTiltSecondBandHighFrequency
									if current_spectral_tilt_value<>undefined
										selectObject: spectral_tilt_values_sliding_window
										Append row
										current_point_index_spectral_tilt = current_point_index_spectral_tilt + 1
										Set numeric value: current_point_index_spectral_tilt, "times", current_point_time
										Set numeric value: current_point_index_spectral_tilt, "values", current_spectral_tilt_value
									endif
								endif
								removeObject: currentSoundSlice, currentSoundSliceSpectrum
								for iband from 1 to nFrequencyBandsForCoGcomputation
									removeObject: currentSoundSliceSpectrum_band'iband'
								endfor
								current_point_time = current_point_time + timestepForCoGcomputationMilliseconds/1000
							endwhile
						endif
						
						if extract_energy_frequency_bands
							for iband from 1 to nFrequencyBandsForCoGcomputation
								selectObject: energy_values_sliding_window_band'iband'
								nrows_energy_values_sliding_window_band'iband' = Get number of rows
								if get_median
									if nrows_energy_values_sliding_window_band'iband'>0
										median_energy_band'iband' = Get quantile: "values", 0.5
									else
										median_energy_band'iband' = undefined
									endif
								endif
								if get_standard_deviation
									if nrows_energy_values_sliding_window_band'iband'>0
										std_energy_band'iband' = Get standard deviation: "values"
									else
										std_energy_band'iband' = undefined
									endif
								endif
								if get_min_max_with_time
									if nrows_energy_values_sliding_window_band'iband'>0
										Sort rows: "values times"
										min_energy_band'iband' = Get value: 1, "values"
										min_energy_band'iband'_relative_time = Get value: 1, "times"
										min_energy_band'iband'_relative_time = (min_energy_band'iband'_relative_time - start_time) / duration
										# reverse sorting: replace values by their negative counterpart (so the first values in time will be selected in cas of equal values)
										Formula: "values", "- (self)"
										Sort rows: "values times"
										max_energy_band'iband' = Get value: 1, "values"
										max_energy_band'iband' = -max_energy_band'iband'
										max_energy_band'iband'_relative_time = Get value: 1, "times"
										max_energy_band'iband'_relative_time = (max_energy_band'iband'_relative_time - start_time) / duration
									else
										min_energy_band'iband' = undefined
										min_energy_band'iband'_relative_time = undefined
										max_energy_band'iband' = undefined
										max_energy_band'iband'_relative_time = undefined
									endif
								endif
								removeObject: energy_values_sliding_window_band'iband'
							endfor
						endif

						if extract_CoG
							selectObject: cog_values_sliding_window
							nrows_cog_values_sliding_window = Get number of rows
							if get_median
								if nrows_cog_values_sliding_window>0
									median_CoG = Get quantile: "values", 0.5
								else
									median_CoG = undefined
								endif
							endif
							if get_standard_deviation
								if nrows_cog_values_sliding_window>0
									std_CoG = Get standard deviation: "values"
								else
									std_CoG = undefined
								endif
							endif
							if get_min_max_with_time
								if nrows_cog_values_sliding_window>0
									Sort rows: "values times"
									min_CoG = Get value: 1, "values"
									min_CoG_relative_time = Get value: 1, "times"
									min_CoG_relative_time = (min_CoG_relative_time - start_time) / duration
									# reverse sorting: replace values by their negative counterpart (so the first values in time will be selected in cas of equal values)
									Formula: "values", "- (self)"
									Sort rows: "values times"
									max_CoG = Get value: 1, "values"
									max_CoG = -max_CoG
									max_CoG_relative_time = Get value: 1, "times"
									max_CoG_relative_time = (max_CoG_relative_time - start_time) / duration
								else
									min_CoG = undefined
									min_CoG_relative_time = undefined
									max_CoG = undefined
									max_CoG_relative_time = undefined
								endif
							endif
							removeObject: cog_values_sliding_window
							for iband from 1 to nFrequencyBandsForCoGcomputation
								selectObject: cog_values_sliding_window_band'iband'
								nrows_cog_values_sliding_window = Get number of rows
								if get_median
									if nrows_cog_values_sliding_window>0
										median_CoG_band'iband' = Get quantile: "values", 0.5
									else
										median_CoG_band'iband' = undefined
									endif
								endif
								if get_standard_deviation
									if nrows_cog_values_sliding_window>0
										std_CoG_band'iband' = Get standard deviation: "values"
									else
										std_CoG_band'iband' = undefined
									endif
								endif
								if get_min_max_with_time
									if nrows_cog_values_sliding_window>0
										Sort rows: "values times"
										min_CoG_band'iband' = Get value: 1, "values"
										min_CoG_relative_time_band'iband' = Get value: 1, "times"
										min_CoG_relative_time_band'iband' = (min_CoG_relative_time_band'iband' - start_time) / duration
										# reverse sorting: replace values by their negative counterpart (so the first values in time will be selected in cas of equal values)
										Formula: "values", "- (self)"
										Sort rows: "values times"
										max_CoG_band'iband' = Get value: 1, "values"
										max_CoG_band'iband' = -max_CoG_band'iband'
										max_CoG_relative_time_band'iband' = Get value: 1, "times"
										max_CoG_relative_time_band'iband' = (max_CoG_relative_time_band'iband' - start_time) / duration
									else
										min_CoG_band'iband' = undefined
										min_CoG_relative_time_band'iband' = undefined
										max_CoG_band'iband' = undefined
										max_CoG_relative_time_band'iband' = undefined
									endif
								endif
								removeObject: cog_values_sliding_window_band'iband'
							endfor
						endif
						if extract_spectral_dispersion
							selectObject: spectral_dispersion_values_sliding_window
							nrows_spectral_dispersion_values_sliding_window = Get number of rows
							if get_median
								if nrows_spectral_dispersion_values_sliding_window>0
									median_spectral_dispersion = Get quantile: "values", 0.5
								else
									median_spectral_dispersion = undefined
								endif
							endif
							if get_standard_deviation
								if nrows_spectral_dispersion_values_sliding_window>0
									std_spectral_dispersion = Get standard deviation: "values"
								else
									std_spectral_dispersion = undefined
								endif
							endif
							if get_min_max_with_time
								if nrows_spectral_dispersion_values_sliding_window>0
									Sort rows: "values times"
									min_spectral_dispersion = Get value: 1, "values"
									min_spectral_dispersion_relative_time = Get value: 1, "times"
									min_spectral_dispersion_relative_time = (min_spectral_dispersion_relative_time - start_time) / duration
									# reverse sorting: replace values by their negative counterpart (so the first values in time will be selected in cas of equal values)
									Formula: "values", "- (self)"
									Sort rows: "values times"
									max_spectral_dispersion = Get value: 1, "values"
									max_spectral_dispersion = -max_spectral_dispersion
									max_spectral_dispersion_relative_time = Get value: 1, "times"
									max_spectral_dispersion_relative_time = (max_spectral_dispersion_relative_time - start_time) / duration
								else
									min_spectral_dispersion = undefined
									min_spectral_dispersion_relative_time = undefined
									min_spectral_dispersion_relative_time = undefined
									max_spectral_dispersion = undefined
									max_spectral_dispersion_relative_time = undefined
									max_spectral_dispersion_relative_time = undefined
								endif
							endif
							removeObject: spectral_dispersion_values_sliding_window
							for iband from 1 to nFrequencyBandsForCoGcomputation
								selectObject: spectral_dispersion_values_sliding_window_band'iband'
								nrows_spectral_dispersion_values_sliding_window = Get number of rows
								if get_median
									if nrows_spectral_dispersion_values_sliding_window>0
										median_spectral_dispersion_band'iband' = Get quantile: "values", 0.5
									else
										median_spectral_dispersion_band'iband' = undefined
									endif
								endif
								if get_standard_deviation
									if nrows_spectral_dispersion_values_sliding_window>0
										std_spectral_dispersion_band'iband' = Get standard deviation: "values"
									else
										std_spectral_dispersion_band'iband' = undefined
									endif
								endif
								if get_min_max_with_time
									if nrows_spectral_dispersion_values_sliding_window>0
										Sort rows: "values times"
										min_spectral_dispersion_band'iband' = Get value: 1, "values"
										min_spectral_dispersion_relative_time_band'iband' = Get value: 1, "times"
										min_spectral_dispersion_relative_time_band'iband' = (min_spectral_dispersion_relative_time_band'iband' - start_time) / duration
										# reverse sorting: replace values by their negative counterpart (so the first values in time will be selected in cas of equal values)
										Formula: "values", "- (self)"
										Sort rows: "values times"
										max_spectral_dispersion_band'iband' = Get value: 1, "values"
										max_spectral_dispersion_band'iband' = -max_spectral_dispersion_band'iband'
										max_spectral_dispersion_relative_time_band'iband' = Get value: 1, "times"
										max_spectral_dispersion_relative_time_band'iband' = (max_spectral_dispersion_relative_time_band'iband' - start_time) / duration
									else
										min_spectral_dispersion_band'iband' = undefined
										min_spectral_dispersion_relative_time_band'iband' = undefined
										max_spectral_dispersion_band'iband' = undefined
										max_spectral_dispersion_relative_time_band'iband' = undefined
									endif
								endif
								removeObject: spectral_dispersion_values_sliding_window_band'iband'
							endfor
						endif
						if extract_spectral_tilt
							selectObject: spectral_tilt_values_sliding_window
							nrows_spectral_tilt_values_sliding_window = Get number of rows
							if get_median
								if nrows_spectral_tilt_values_sliding_window>0
									median_spectral_tilt = Get quantile: "values", 0.5
								else
									median_spectral_tilt = undefined
								endif
							endif
							if get_standard_deviation
								if nrows_spectral_tilt_values_sliding_window>0
									std_spectral_tilt = Get standard deviation: "values"
								else
									std_spectral_tilt = undefined
								endif
							endif
							if get_min_max_with_time
								if nrows_spectral_tilt_values_sliding_window>0
									Sort rows: "values times"
									min_spectral_tilt = Get value: 1, "values"
									min_spectral_tilt_relative_time = Get value: 1, "times"
									min_spectral_tilt_relative_time = (min_spectral_tilt_relative_time - start_time) / duration
									# reverse sorting: replace values by their negative counterpart (so the first values in time will be selected in cas of equal values)
									Formula: "values", "- (self)"
									Sort rows: "values times"
									max_spectral_tilt = Get value: 1, "values"
									max_spectral_tilt = -max_spectral_tilt
									max_spectral_tilt_relative_time = Get value: 1, "times"
									max_spectral_tilt_relative_time = (max_spectral_tilt_relative_time - start_time) / duration
								else
									min_spectral_tilt = undefined
									min_spectral_tilt_relative_time = undefined
									min_spectral_tilt_relative_time = undefined
									max_spectral_tilt = undefined
									max_spectral_tilt_relative_time = undefined
									max_spectral_tilt_relative_time = undefined
								endif
							endif
							removeObject: spectral_tilt_values_sliding_window
						endif
					endif
					
					if extract_HNR
						selectObject: current_snd_extract
						if minFrequencyUsedInBandPassFilterForHNRcomputation=0 and maxFrequencyUsedInBandPassFilterForHNRcomputation=undefined
							# Keep sound extract without prior filtering
							current_snd_extract_filterForHNR = Copy: "extract_copy_HNR"
						elsif minFrequencyUsedInBandPassFilterForHNRcomputation>0 and maxFrequencyUsedInBandPassFilterForHNRcomputation=undefined
							# High-pass filtering
							current_snd_extract_filterForHNR = Filter (stop Hann band): 0, minFrequencyUsedInBandPassFilterForHNRcomputation, smoothingUsedInBandPassFilterForHNRcomputation
						else
							# Band-pass filtering
							current_snd_extract_filterForHNR = Filter (pass Hann band): minFrequencyUsedInBandPassFilterForHNRcomputation, maxFrequencyUsedInBandPassFilterForHNRcomputation, smoothingUsedInBandPassFilterForHNRcomputation
						endif

						selectObject: current_snd_extract_filterForHNR
						current_HNR = noprogress To Harmonicity (cc): timeStepHNRextraction, minF0, silenceTresholdHNRextraction, periodsPerWindowHNRextraction
						hnr_values_extracted_points = Create TableOfReal: "hnr_values", n_extracted_points, 1
						for ipoint from 1 to n_extracted_points
							selectObject: extracted_points_relative_times
							current_point_relative_position = Get value in cell: ipoint, 1
							current_point_time = start_time + (current_point_relative_position * duration)
							selectObject: current_HNR
							current_value = Get value at time: current_point_time, interpolationMethodUsedInHNRextraction$
							selectObject: hnr_values_extracted_points
							Set value: ipoint, 1, current_value
						endfor
						if get_mean
							selectObject: current_HNR
							mean_HNR = Get mean: start_time, end_time
						endif
						if get_median
							# No "Get quantile" function defined for Harmonicity objects: export values as Table
							selectObject: current_HNR
							mat_HNR_tmp = To Matrix
							tmat_HNR_tmp = Transpose
							tor_HNR_tmp = To TableOfReal
							table_HNR_tmp = To Table: "rowLabel"
							Set column label (index): 2, "value"
							median_HNR = Get quantile: "value", 0.5
							removeObject: mat_HNR_tmp, tmat_HNR_tmp, tor_HNR_tmp, table_HNR_tmp
						endif
						if get_standard_deviation
							selectObject: current_HNR
							std_HNR = Get standard deviation: start_time, end_time
						endif
						if get_min_max_with_time
							selectObject: current_HNR
							min_HNR = Get minimum: start_time, end_time, interpolationMethodUsedInHNRextraction$
							min_HNR_relative_time = Get time of minimum: start_time, end_time, interpolationMethodUsedInHNRextraction$
							min_HNR_relative_time = (min_HNR_relative_time - start_time) / duration
							max_HNR = Get maximum: start_time, end_time, interpolationMethodUsedInHNRextraction$
							max_HNR_relative_time = Get time of maximum: start_time, end_time, interpolationMethodUsedInHNRextraction$
							max_HNR_relative_time = (max_HNR_relative_time - start_time) / duration
						endif
						removeObject: current_HNR, current_snd_extract_filterForHNR
					endif
					
					if extract_ZCR
						selectObject: current_snd_extract
						current_zcr_pp = noprogress To PointProcess (zeroes): 1, includeRaisingPartsInZeroCrossingRateComputation$, includeFallingPartsInZeroCrossingRateComputation$
						current_zcr_tt = Up to TextTier: ""
						current_zcr_tg = Into TextGrid
						removeObject: current_zcr_pp, current_zcr_tt
						
						zcr_values_extracted_points = Create TableOfReal: "zcr_values", n_extracted_points, 1
						for ipoint from 1 to n_extracted_points
							selectObject: extracted_points_relative_times
							current_point_relative_position = Get value in cell: ipoint, 1
							current_point_time = start_time + (current_point_relative_position * duration)
							# Get number of points in extract
							current_extract_start_time = current_point_time - (extractsDurationForZCRcomputationMilliseconds/2000)
							if current_extract_start_time < start_time - (offset_for_acoustic_parameters_extraction_milliseconds/1000)
								current_extract_start_time = start_time - (offset_for_acoustic_parameters_extraction_milliseconds/1000)
							endif
							current_extract_end_time = current_point_time + (extractsDurationForZCRcomputationMilliseconds/2000)
							if current_extract_end_time > end_time + (offset_for_acoustic_parameters_extraction_milliseconds/1000)
								current_extract_end_time = end_time - (offset_for_acoustic_parameters_extraction_milliseconds/1000)
							endif
							selectObject: current_zcr_tg
							current_zcr_tg_extract = Extract part: current_extract_start_time, current_extract_end_time, "no"
							nCrossingPoints = Get number of points: 1
							current_value = nCrossingPoints / (current_extract_end_time - current_extract_start_time)
							selectObject: zcr_values_extracted_points
							Set value: ipoint, 1, current_value
							removeObject: current_zcr_tg_extract
						endfor
						if get_mean
							selectObject: current_zcr_tg
							zcr_tg_extract = Extract part: start_time, end_time, "no"
							nCrossingPoints = Get number of points: 1
							mean_ZCR = nCrossingPoints / (end_time - start_time)
							removeObject: zcr_tg_extract
						endif

						# If the computation of the median, standard deviation or minumum/maximum is requested, compute values using a sliding window and store them in a Table object
						zcr_values_sliding_window = Create Table with column names: "table_ZCR_values", 0, "times values"
						if get_median+get_standard_deviation+get_min_max_with_time > 0
							current_point_index = 0
							current_point_time = start_time
							while current_point_time <= end_time
								current_extract_start_time = current_point_time - (extractsDurationForZCRcomputationMilliseconds/2000)
								if current_extract_start_time < start_time - (offset_for_acoustic_parameters_extraction_milliseconds/1000)
									current_extract_start_time = start_time - (offset_for_acoustic_parameters_extraction_milliseconds/1000)
								endif
								current_extract_end_time = current_point_time + (extractsDurationForZCRcomputationMilliseconds/2000)
								if current_extract_end_time > end_time + (offset_for_acoustic_parameters_extraction_milliseconds/1000)
									current_extract_end_time = end_time - (offset_for_acoustic_parameters_extraction_milliseconds/1000)
								endif
								selectObject: current_zcr_tg
								current_zcr_tg_extract = Extract part: current_extract_start_time, current_extract_end_time, "no"
								nCrossingPoints = Get number of points: 1
								current_value = nCrossingPoints / (current_extract_end_time - current_extract_start_time)
								if current_value<>undefined
									selectObject: zcr_values_sliding_window
									Append row
									current_point_index = current_point_index + 1
									Set numeric value: current_point_index, "times", current_point_time
									Set numeric value: current_point_index, "values", current_value
								endif
								removeObject: current_zcr_tg_extract
								current_point_time = current_point_time + timestepForZCRcomputationMilliseconds/1000
							endwhile
						endif
						
						selectObject: zcr_values_sliding_window					
						nrows_zcr_values_sliding_window = Get number of rows
						if get_median
							if nrows_zcr_values_sliding_window>0
								median_ZCR = Get quantile: "values", 0.5
							else
								median_ZCR = undefined
							endif
						endif
						if get_standard_deviation
							if nrows_zcr_values_sliding_window>0
								std_ZCR = Get standard deviation: "values"
							else
								std_ZCR = undefined
							endif
						endif
						if get_min_max_with_time
							if nrows_zcr_values_sliding_window>0
								Sort rows: "values times"
								min_ZCR = Get value: 1, "values"
								min_ZCR_relative_time = Get value: 1, "times"
								min_ZCR_relative_time = (min_ZCR_relative_time - start_time) / duration
								# reverse sorting: replace values by their negative counterpart (so the first values in time will be selected in cas of equal values)
								Formula: "values", "- (self)"
								Sort rows: "values times"
								max_ZCR = Get value: 1, "values"
								max_ZCR = -max_ZCR
								max_ZCR_relative_time = Get value: 1, "times"
								max_ZCR_relative_time = (max_ZCR_relative_time - start_time) / duration
							else
								min_ZCR = undefined
								min_ZCR_relative_time = undefined
								min_ZCR_relative_time = undefined
								max_ZCR = undefined
								max_ZCR_relative_time = undefined
								max_ZCR_relative_time = undefined
							endif
						endif
						removeObject: current_zcr_tg, zcr_values_sliding_window
					endif
					
					if extract_formants
						selectObject: current_snd_extract
						current_formant = noprogress To Formant (burg): timeStepFormantsDetection, nDetectedFormants, maxFrequencyForFormantsDetection, windowLengthFormantsDetectionSeconds, preEmphasisFormantsDetectionHertz
						f1_values_extracted_points = Create TableOfReal: "f1_values", n_extracted_points, 1
						for ipoint from 1 to n_extracted_points
							selectObject: extracted_points_relative_times
							current_point_relative_position = Get value in cell: ipoint, 1
							current_point_time = start_time + (current_point_relative_position * duration)
							selectObject: current_formant
							current_value = Get value at time: 1, current_point_time, "Hertz", "Linear"
							selectObject: f1_values_extracted_points
							Set value: ipoint, 1, current_value
						endfor
						selectObject: current_formant
						if get_mean
							mean_f1 = Get mean: 1, start_time, end_time, "Hertz"
						endif
						if get_median
							median_f1 = Get quantile: 1, start_time, end_time, "hertz", 0.5
						endif
						if get_standard_deviation
							std_f1 = Get standard deviation: 1, start_time, end_time, "Hertz"
						endif
						if get_min_max_with_time
							min_f1 = Get minimum: 1, start_time, end_time, "Hertz", "Parabolic"
							min_f1_relative_time = Get time of minimum: 1, start_time, end_time, "Hertz", "Parabolic"
							min_f1_relative_time = (min_f1_relative_time - start_time) / duration
							max_f1 = Get maximum: 1, start_time, end_time, "Hertz", "Parabolic"
							max_f1_relative_time = Get time of maximum: 1, start_time, end_time, "Hertz", "Parabolic"
							max_f1_relative_time = (max_f1_relative_time - start_time) / duration
						endif

						f2_values_extracted_points = Create TableOfReal: "f2_values", n_extracted_points, 1
						for ipoint from 1 to n_extracted_points
							selectObject: extracted_points_relative_times
							current_point_relative_position = Get value in cell: ipoint, 1
							current_point_time = start_time + (current_point_relative_position * duration)
							selectObject: current_formant
							current_value = Get value at time: 2, current_point_time, "Hertz", "Linear"
							selectObject: f2_values_extracted_points
							Set value: ipoint, 1, current_value
						endfor
						selectObject: current_formant
						if get_mean
							mean_f2 = Get mean: 2, start_time, end_time, "Hertz"
						endif
						if get_median
							median_f2 = Get quantile: 2, start_time, end_time, "hertz", 0.5
						endif
						if get_standard_deviation
							std_f2 = Get standard deviation: 2, start_time, end_time, "Hertz"
						endif
						if get_min_max_with_time
							min_f2 = Get minimum: 2, start_time, end_time, "Hertz", "Parabolic"
							min_f2_relative_time = Get time of minimum: 2, start_time, end_time, "Hertz", "Parabolic"
							min_f2_relative_time = (min_f2_relative_time - start_time) / duration
							max_f2 = Get maximum: 2, start_time, end_time, "Hertz", "Parabolic"
							max_f2_relative_time = Get time of maximum: 2, start_time, end_time, "Hertz", "Parabolic"
							max_f2_relative_time = (max_f2_relative_time - start_time) / duration
						endif

						f3_values_extracted_points = Create TableOfReal: "f3_values", n_extracted_points, 1
						for ipoint from 1 to n_extracted_points
							selectObject: extracted_points_relative_times
							current_point_relative_position = Get value in cell: ipoint, 1
							current_point_time = start_time + (current_point_relative_position * duration)
							selectObject: current_formant
							current_value = Get value at time: 3, current_point_time, "Hertz", "Linear"
							selectObject: f3_values_extracted_points
							Set value: ipoint, 1, current_value
						endfor
						selectObject: current_formant
						if get_mean
							mean_f3 = Get mean: 3, start_time, end_time, "Hertz"
						endif
						if get_median
							median_f3 = Get quantile: 3, start_time, end_time, "hertz", 0.5
						endif
						if get_standard_deviation
							std_f3 = Get standard deviation: 3, start_time, end_time, "Hertz"
						endif
						if get_min_max_with_time
							min_f3 = Get minimum: 3, start_time, end_time, "Hertz", "Parabolic"
							min_f3_relative_time = Get time of minimum: 3, start_time, end_time, "Hertz", "Parabolic"
							min_f3_relative_time = (min_f3_relative_time - start_time) / duration
							max_f3 = Get maximum: 3, start_time, end_time, "Hertz", "Parabolic"
							max_f3_relative_time = Get time of maximum: 3, start_time, end_time, "Hertz", "Parabolic"
							max_f3_relative_time = (max_f3_relative_time - start_time) / duration
						endif

						f4_values_extracted_points = Create TableOfReal: "f4_values", n_extracted_points, 1
						for ipoint from 1 to n_extracted_points
							selectObject: extracted_points_relative_times
							current_point_relative_position = Get value in cell: ipoint, 1
							current_point_time = start_time + (current_point_relative_position * duration)
							selectObject: current_formant
							current_value = Get value at time: 4, current_point_time, "Hertz", "Linear"
							selectObject: f4_values_extracted_points
							Set value: ipoint, 1, current_value
						endfor
						selectObject: current_formant
						if get_mean
							mean_f4 = Get mean: 4, start_time, end_time, "Hertz"
						endif
						if get_median
							median_f4 = Get quantile: 4, start_time, end_time, "hertz", 0.5
						endif
						if get_standard_deviation
							std_f4 = Get standard deviation: 4, start_time, end_time, "Hertz"
						endif
						if get_min_max_with_time
							min_f4 = Get minimum: 4, start_time, end_time, "Hertz", "Parabolic"
							min_f4_relative_time = Get time of minimum: 4, start_time, end_time, "Hertz", "Parabolic"
							min_f4_relative_time = (min_f4_relative_time - start_time) / duration
							max_f4 = Get maximum: 4, start_time, end_time, "Hertz", "Parabolic"
							max_f4_relative_time = Get time of maximum: 4, start_time, end_time, "Hertz", "Parabolic"
							max_f4_relative_time = (max_f4_relative_time - start_time) / duration
						endif

						removeObject: current_formant
					endif

					if extract_MFCC
						selectObject: current_snd_extract
						current_mfcc = noprogress To MFCC: numberOfMFCCcoefficients, windowLengthMFCC, timeStepMFCC, firstFilterFrequencyMFCC, distanceBetweenFiltersMFCC, maximumFrequencyMFCC
						current_mfcc_nframes = Get number of frames
						starting_frame_index_MFCC = Get frame number from time: start_time
						starting_frame_index_MFCC = round(starting_frame_index_MFCC)
						if starting_frame_index_MFCC<1
							starting_frame_index_MFCC = 1
						endif
						if starting_frame_index_MFCC>current_mfcc_nframes
							starting_frame_index_MFCC = current_mfcc_nframes
						endif
						end_frame_index_MFCC = Get frame number from time: end_time
						end_frame_index_MFCC = round(end_frame_index_MFCC)
						if end_frame_index_MFCC>current_mfcc_nframes
							end_frame_index_MFCC = current_mfcc_nframes
						endif


						mfcc_values_extracted_points = Create TableOfReal: "mfcc_values", n_extracted_points, numberOfMFCCcoefficients+1
						for ipoint from 1 to n_extracted_points
							selectObject: extracted_points_relative_times
							current_point_relative_position = Get value in cell: ipoint, 1
							current_point_time = start_time + (current_point_relative_position * duration)
							selectObject: current_mfcc
							current_frame_index_MFCC = Get frame number from time: current_point_time
							current_frame_index_MFCC = round(current_frame_index_MFCC)
							if current_frame_index_MFCC<1
								current_frame_index_MFCC = 1
							endif
							if current_frame_index_MFCC>current_mfcc_nframes
								current_frame_index_MFCC = current_mfcc_nframes
							endif
							# c0
							current_value = Get c0 value in frame: current_frame_index_MFCC
							selectObject: mfcc_values_extracted_points
							Set value: ipoint, 1, current_value
							for iMFCC from 1 to numberOfMFCCcoefficients
								selectObject: current_mfcc
								current_value = Get value in frame: current_frame_index_MFCC, iMFCC
								selectObject: mfcc_values_extracted_points
								Set value: ipoint, iMFCC+1, current_value
							endfor
						endfor
						if get_mean+get_median+get_standard_deviation+get_min_max_with_time > 0
							selectObject: current_mfcc
							tor_MFCC_tmp = To TableOfReal: "yes"
							extractedRowsRange$ = fixed$(starting_frame_index_MFCC,0)+":"+fixed$(end_frame_index_MFCC,0)
							tor_MFCC_current_snd_extract = Extract row ranges: extractedRowsRange$
							removeObject: tor_MFCC_tmp
						endif

						if get_mean
							for iMFCC from 0 to numberOfMFCCcoefficients
								selectObject: tor_MFCC_current_snd_extract
								mean_MFCC'iMFCC' = Get column mean (index): iMFCC+1
							endfor
						endif
						if get_standard_deviation
							for iMFCC from 0 to numberOfMFCCcoefficients
								selectObject: tor_MFCC_current_snd_extract
								std_MFCC'iMFCC' = Get column stdev (index): iMFCC+1
							endfor
						endif

						if get_median+get_min_max_with_time > 0
							selectObject: tor_MFCC_current_snd_extract
							table_tor_MFCC_current_snd_extract = To Table: "rowLabel"
							for iMFCC from 0 to numberOfMFCCcoefficients
								if get_median
									selectObject: table_tor_MFCC_current_snd_extract
									median_MFCC'iMFCC' = Get quantile: "c"+fixed$(iMFCC,0), 0.5
								endif
								if get_min_max_with_time
									selectObject: table_tor_MFCC_current_snd_extract
									minValue = Get minimum: "c"+fixed$(iMFCC,0)
									min_MFCC'iMFCC' = minValue
									iFrameMinValue = Search column: "c"+fixed$(iMFCC,0), "'minValue'"
									maxValue = Get maximum: "c"+fixed$(iMFCC,0)
									max_MFCC'iMFCC' = maxValue
									iFrameMaxValue = Search column: "c"+fixed$(iMFCC,0), "'maxValue'"
									selectObject: current_mfcc
									min_MFCC'iMFCC'_absolute_time = Get time from frame number: iFrameMinValue+starting_frame_index_MFCC-1
									max_MFCC'iMFCC'_absolute_time = Get time from frame number: iFrameMaxValue+starting_frame_index_MFCC-1
									min_MFCC'iMFCC'_relative_time = (min_MFCC'iMFCC'_absolute_time - start_time) / duration
									max_MFCC'iMFCC'_relative_time = (max_MFCC'iMFCC'_absolute_time - start_time) / duration
								endif
							endfor
							removeObject: table_tor_MFCC_current_snd_extract
						endif

						if get_mean+get_median+get_standard_deviation+get_min_max_with_time > 0
							removeObject: tor_MFCC_current_snd_extract
						endif
						removeObject: current_mfcc
					endif

					if extract_CPP
						# Analysis of CPP and CPPS on each target point
						cpp_values_extracted_points = Create TableOfReal: "cpp_values", n_extracted_points, 1
						cpps_values_extracted_points = Create TableOfReal: "cpp_values", n_extracted_points, 1
						for ipoint from 1 to n_extracted_points
							selectObject: extracted_points_relative_times
							current_point_relative_position = Get value in cell: ipoint, 1
							current_point_time = start_time + (current_point_relative_position * duration)
							# Get number of points in extract
							current_extract_start_time = current_point_time - (extractsDurationForCPPcomputationMilliseconds/2000)
							if current_extract_start_time < start_time - (offset_for_acoustic_parameters_extraction_milliseconds/1000)
								current_extract_start_time = start_time - (offset_for_acoustic_parameters_extraction_milliseconds/1000)
							endif
							current_extract_end_time = current_point_time + (extractsDurationForCPPcomputationMilliseconds/2000)
							if current_extract_end_time > end_time + (offset_for_acoustic_parameters_extraction_milliseconds/1000)
								current_extract_end_time = end_time - (offset_for_acoustic_parameters_extraction_milliseconds/1000)
							endif
							# Get sound extract for current frame
							selectObject: current_snd_extract
							current_cpp_snd_extract = Extract part: current_extract_start_time, current_extract_end_time, "rectangular", 1, "no"
							# Compute CPP
							selectObject: current_cpp_snd_extract
							currentFrameSpectrum = noprogress To Spectrum: "yes"
							currentFramePowerCepstrum = noprogress To PowerCepstrum
							selectObject: currentFramePowerCepstrum
							currentFrameCPP = Get peak prominence: minF0forPeakProminenceComputation, maxF0forPeakProminenceComputation, interpolationMethodPeakProminenceComputation$, trendLineQuefrencyMinValuePeakProminenceComputation, trendLineQuefrencyMaxValuePeakProminenceComputation, trendTypePeakProminenceComputation$, fitMethodPeakProminenceComputation$
							selectObject: cpp_values_extracted_points
							Set value: ipoint, 1, currentFrameCPP
							# Compute CPPS
							selectObject: current_cpp_snd_extract
							currentFramePowerCepstrogram = noprogress To PowerCepstrogram: minF0, timeStepPowerCepstrogram, maxFrequencyPowerCepstrogram, preEmphasisStartFrequencyPowerCepstrogram
							selectObject: currentFramePowerCepstrogram
							currentFrameCPPS = Get CPPS: subtractTrendBeforeSmoothingCPPS$, timeAveragingWindowCPPS, quefrencyAveragingWindowCPPS, minF0forPeakProminenceComputation, maxF0forPeakProminenceComputation, peakSearchToleranceFactorCPPS, interpolationMethodPeakProminenceComputation$, trendLineQuefrencyMinValuePeakProminenceComputation, trendLineQuefrencyMaxValuePeakProminenceComputation, trendTypePeakProminenceComputation$, fitMethodPeakProminenceComputation$
							selectObject: cpps_values_extracted_points
							Set value: ipoint, 1, currentFrameCPPS
							removeObject: current_cpp_snd_extract, currentFrameSpectrum, currentFramePowerCepstrum, currentFramePowerCepstrogram
						endfor
						# If the computation of the median, standard deviation or minumum/maximum is requested, compute values using a sliding window and store them in a Table object
						cpp_values_sliding_window = Create Table with column names: "table_CPP_values", 0, "times values"
						cpps_values_sliding_window = Create Table with column names: "table_CPP_values", 0, "times values"
						if get_mean+get_median+get_standard_deviation+get_min_max_with_time > 0
							current_point_index = 0
							current_point_time = start_time
							while current_point_time <= end_time
								current_extract_start_time = current_point_time - (extractsDurationForCPPcomputationMilliseconds/2000)
								if current_extract_start_time < start_time - (offset_for_acoustic_parameters_extraction_milliseconds/1000)
									current_extract_start_time = start_time - (offset_for_acoustic_parameters_extraction_milliseconds/1000)
								endif
								current_extract_end_time = current_point_time + (extractsDurationForCPPcomputationMilliseconds/2000)
								if current_extract_end_time > end_time + (offset_for_acoustic_parameters_extraction_milliseconds/1000)
									current_extract_end_time = end_time - (offset_for_acoustic_parameters_extraction_milliseconds/1000)
								endif
								# Get sound extract for current frame
								selectObject: current_snd_extract
								current_cpp_snd_extract = Extract part: current_extract_start_time, current_extract_end_time, "rectangular", 1, "no"
								# Compute CPP
								selectObject: current_cpp_snd_extract
								currentFrameSpectrum = noprogress To Spectrum: "yes"
								currentFramePowerCepstrum = noprogress To PowerCepstrum
								selectObject: currentFramePowerCepstrum
								currentFrameCPP = Get peak prominence: minF0forPeakProminenceComputation, maxF0forPeakProminenceComputation, interpolationMethodPeakProminenceComputation$, trendLineQuefrencyMinValuePeakProminenceComputation, trendLineQuefrencyMaxValuePeakProminenceComputation, trendTypePeakProminenceComputation$, fitMethodPeakProminenceComputation$
								if current_value<>undefined
									selectObject: cpp_values_sliding_window
									Append row
									current_point_index = current_point_index + 1
									Set numeric value: current_point_index, "times", current_point_time
									Set numeric value: current_point_index, "values", current_value
								endif
								# Compute CPPS
								selectObject: current_cpp_snd_extract
								currentFramePowerCepstrogram = noprogress To PowerCepstrogram: minF0, timeStepPowerCepstrogram, maxFrequencyPowerCepstrogram, preEmphasisStartFrequencyPowerCepstrogram
								selectObject: currentFramePowerCepstrogram
								currentFrameCPPS = Get CPPS: subtractTrendBeforeSmoothingCPPS$, timeAveragingWindowCPPS, quefrencyAveragingWindowCPPS, minF0forPeakProminenceComputation, maxF0forPeakProminenceComputation, peakSearchToleranceFactorCPPS, interpolationMethodPeakProminenceComputation$, trendLineQuefrencyMinValuePeakProminenceComputation, trendLineQuefrencyMaxValuePeakProminenceComputation, trendTypePeakProminenceComputation$, fitMethodPeakProminenceComputation$
								if current_value<>undefined
									selectObject: cpps_values_sliding_window
									Append row
									current_point_index = current_point_index + 1
									Set numeric value: current_point_index, "times", current_point_time
									Set numeric value: current_point_index, "values", current_value
								endif
								current_point_time = current_point_time + timestepForZCRcomputationMilliseconds/1000
								removeObject: current_cpp_snd_extract, currentFrameSpectrum, currentFramePowerCepstrum, currentFramePowerCepstrogram
							endwhile
						endif
						
						selectObject: cpp_values_sliding_window					
						nrows_cpp_values_sliding_window = Get number of rows
						selectObject: cpps_values_sliding_window					
						nrows_cpps_values_sliding_window = Get number of rows
						if get_mean
							if nrows_cpp_values_sliding_window>0
								selectObject: cpp_values_sliding_window	
								mean_CPP = Get mean: "values"
							else
								mean_CPP = undefined
							endif
							if nrows_cpps_values_sliding_window>0
								selectObject: cpps_values_sliding_window	
								mean_CPP = Get mean: "values"
							else
								mean_CPP = undefined
							endif
						endif
						if get_median
							if nrows_cpp_values_sliding_window>0
								selectObject: cpp_values_sliding_window	
								median_CPP = Get quantile: "values", 0.5
							else
								median_CPP = undefined
							endif
							if nrows_cpps_values_sliding_window>0
								selectObject: cpps_values_sliding_window	
								median_CPP = Get quantile: "values", 0.5
							else
								median_CPP = undefined
							endif
						endif
						if get_standard_deviation
							if nrows_cpp_values_sliding_window>0
								selectObject: cpp_values_sliding_window	
								std_CPP = Get standard deviation: "values"
							else
								std_CPP = undefined
							endif
							if nrows_cpps_values_sliding_window>0
								selectObject: cpps_values_sliding_window	
								std_CPPS = Get standard deviation: "values"
							else
								std_CPPS = undefined
							endif
						endif
						if get_min_max_with_time
							if nrows_cpp_values_sliding_window>0
								selectObject: cpp_values_sliding_window
								Sort rows: "values times"
								min_CPP = Get value: 1, "values"
								min_CPP_relative_time = Get value: 1, "times"
								min_CPP_relative_time = (min_CPP_relative_time - start_time) / duration
								# reverse sorting: replace values by their negative counterpart (so the first values in time will be selected in cas of equal values)
								Formula: "values", "- (self)"
								Sort rows: "values times"
								max_CPP = Get value: 1, "values"
								max_CPP = -max_CPP
								max_CPP_relative_time = Get value: 1, "times"
								max_CPP_relative_time = (max_CPP_relative_time - start_time) / duration
							else
								min_CPP = undefined
								min_CPP_relative_time = undefined
								min_CPP_relative_time = undefined
								max_CPP = undefined
								max_CPP_relative_time = undefined
								max_CPP_relative_time = undefined
							endif
							if nrows_cpps_values_sliding_window>0
								selectObject: cpps_values_sliding_window
								Sort rows: "values times"
								min_CPPS = Get value: 1, "values"
								min_CPPS_relative_time = Get value: 1, "times"
								min_CPPS_relative_time = (min_CPPS_relative_time - start_time) / duration
								# reverse sorting: replace values by their negative counterpart (so the first values in time will be selected in cas of equal values)
								Formula: "values", "- (self)"
								Sort rows: "values times"
								max_CPPS = Get value: 1, "values"
								max_CPPS = -max_CPPS
								max_CPPS_relative_time = Get value: 1, "times"
								max_CPPS_relative_time = (max_CPPS_relative_time - start_time) / duration
							else
								min_CPPS = undefined
								min_CPPS_relative_time = undefined
								min_CPPS_relative_time = undefined
								max_CPPS = undefined
								max_CPPS_relative_time = undefined
								max_CPPS_relative_time = undefined
							endif
						endif
						removeObject: cpp_values_sliding_window, cpps_values_sliding_window
					endif

					if extract_signal_value
						selectObject: current_snd_extract
						signal_value_extracted_points = Create TableOfReal: "signal_values", n_extracted_points, 1
						for ipoint from 1 to n_extracted_points
							selectObject: extracted_points_relative_times
							current_point_relative_position = Get value in cell: ipoint, 1
							current_point_time = start_time + (current_point_relative_position * duration)
							selectObject: current_snd_extract
							current_value = Get value at time: 1, current_point_time, "Sinc70"
							selectObject: signal_value_extracted_points
							Set value: ipoint, 1, current_value
						endfor
						if get_mean
							selectObject: current_snd_extract
							mean_signal_value = Get mean: 1, start_time, end_time
						endif
						if get_median
							# Trick to get access to the "Get quantile" function (not available for Sound objects):
							# convert to Pitch via a Matrix object
							selectObject: current_snd_extract
							signal_matrix_tmp = Down to Matrix
							signal_pitch_obj_tmp = noprogress To Pitch
							median_signal_value = Get quantile: start_time, end_time, 0.5, "Hertz"
							removeObject: signal_matrix_tmp, signal_pitch_obj_tmp
						endif
						if get_standard_deviation
							selectObject: current_snd_extract
							std_signal_value = Get standard deviation: 1, start_time, end_time
						endif
						if get_min_max_with_time
							selectObject: current_snd_extract
							min_signal_value = Get minimum: start_time, end_time, "Sinc70"
							min_signal_value_relative_time = Get time of minimum: start_time, end_time, "Sinc70"
							min_signal_value_relative_time = (min_signal_value_relative_time - start_time) / duration
							max_signal_value = Get maximum: start_time, end_time, "Sinc70"
							max_signal_value_relative_time = Get time of maximum: start_time, end_time, "Sinc70"
							max_signal_value_relative_time = (max_signal_value_relative_time - start_time) / duration
						endif
					endif

					removeObject: current_snd_extract
				endif

				# Write information to results file
				appendFile: results_file$, tg$
				if export_sound_files_name
					appendFile: results_file$, tab$, snd$
				endif
				appendFile: results_file$, tab$, label$, tab$, start_time, tab$, end_time, tab$, duration
				if extract_left_and_right_context
					appendFile: results_file$, tab$, previousLabel$, tab$, followingLabel$
				endif
				# F0
				if extract_F0
					selectObject: f0_values_extracted_points
					for ipoint from 1 to n_extracted_points
						current_point_value = Get value: ipoint, 1
						appendFile: results_file$, tab$, current_point_value
					endfor
					if get_mean
						appendFile: results_file$, tab$, mean_f0
					endif
					if get_median
						appendFile: results_file$, tab$, median_f0
					endif
					if get_standard_deviation
						appendFile: results_file$, tab$, std_f0
					endif
					if get_min_max_with_time
						appendFile: results_file$, tab$, min_f0, tab$, min_f0_relative_time, tab$, max_f0, tab$, max_f0_relative_time
					endif
				endif
				# intensity
				if extract_intensity
					selectObject: intensity_values_extracted_points
					for ipoint from 1 to n_extracted_points
						current_point_value = Get value: ipoint, 1
						appendFile: results_file$, tab$, current_point_value
					endfor
					if get_mean
						appendFile: results_file$, tab$, mean_intensity
					endif
					if get_median
						appendFile: results_file$, tab$, median_intensity
					endif
					if get_standard_deviation
						appendFile: results_file$, tab$, std_intensity
					endif
					if get_min_max_with_time
						appendFile: results_file$, tab$, min_intensity, tab$, min_intensity_relative_time, tab$, max_intensity, tab$, max_intensity_relative_time
					endif
				endif
				if extract_energy_frequency_bands
					for iband from 1 to nFrequencyBandsForCoGcomputation
						if get_mean
							appendFile: results_file$, tab$, mean_energy_band'iband'
						endif
						if get_median
							appendFile: results_file$, tab$, median_energy_band'iband'
						endif
						if get_standard_deviation
							appendFile: results_file$, tab$, std_energy_band'iband'
						endif
						if get_min_max_with_time
							appendFile: results_file$, tab$, min_energy_band'iband'
							appendFile: results_file$, tab$, min_energy_band'iband'_relative_time
							appendFile: results_file$, tab$, max_energy_band'iband'
							appendFile: results_file$, tab$, max_energy_band'iband'_relative_time
						endif
					endfor
					for ipoint from 1 to n_extracted_points
						for iband from 1 to nFrequencyBandsForCoGcomputation
							selectObject: energy_extracted_points_band'iband'
							current_point_value = Get value: ipoint, 1
							appendFile: results_file$, tab$, current_point_value
						endfor
					endfor
				endif
				# CoG
				if extract_CoG
					for ipoint from 1 to n_extracted_points
						selectObject: cog_values_extracted_points
						current_point_value = Get value: ipoint, 1
						appendFile: results_file$, tab$, current_point_value
						if extract_energy_frequency_bands
							for iband from 1 to nFrequencyBandsForCoGcomputation
								selectObject: cog_values_extracted_points_band'iband'
								current_point_value = Get value: ipoint, 1
								appendFile: results_file$, tab$, current_point_value
							endfor
						endif
					endfor
					if get_mean
						appendFile: results_file$, tab$, mean_CoG
						if extract_energy_frequency_bands
							for iband from 1 to nFrequencyBandsForCoGcomputation
								appendFile: results_file$, tab$, mean_CoG_band'iband'
							endfor
						endif
					endif
					if get_median
						appendFile: results_file$, tab$, median_CoG
						if extract_energy_frequency_bands
							for iband from 1 to nFrequencyBandsForCoGcomputation
								appendFile: results_file$, tab$, median_CoG_band'iband'
							endfor
						endif
					endif
					if get_standard_deviation
						appendFile: results_file$, tab$, std_CoG
						if extract_energy_frequency_bands
							for iband from 1 to nFrequencyBandsForCoGcomputation
								appendFile: results_file$, tab$, std_CoG_band'iband'
							endfor
						endif
					endif
					if get_min_max_with_time
						appendFile: results_file$, tab$, min_CoG, tab$, min_CoG_relative_time, tab$, max_CoG, tab$, max_CoG_relative_time
						if extract_energy_frequency_bands
							for iband from 1 to nFrequencyBandsForCoGcomputation
								appendFile: results_file$, tab$, min_CoG_band'iband'
								appendFile: results_file$, tab$, min_CoG_relative_time_band'iband'
								appendFile: results_file$, tab$, max_CoG_band'iband'
								appendFile: results_file$, tab$, max_CoG_relative_time_band'iband'
							endfor
						endif
					endif
				endif
				# spectral dispersion
				if extract_spectral_dispersion
					for ipoint from 1 to n_extracted_points
						selectObject: spectral_dispersion_values_extracted_points
						current_point_value = Get value: ipoint, 1
						appendFile: results_file$, tab$, current_point_value
						if extract_energy_frequency_bands
							for iband from 1 to nFrequencyBandsForCoGcomputation
								selectObject: spectral_dispersion_values_extracted_points_band'iband'
								current_point_value = Get value: ipoint, 1
								appendFile: results_file$, tab$, current_point_value
							endfor
						endif
					endfor
					if get_mean
						appendFile: results_file$, tab$, mean_spectral_dispersion
						if extract_energy_frequency_bands
							for iband from 1 to nFrequencyBandsForCoGcomputation
								appendFile: results_file$, tab$, mean_spectral_dispersion_band'iband'
							endfor
						endif
					endif
					if get_median
						appendFile: results_file$, tab$, median_spectral_dispersion
						if extract_energy_frequency_bands
							for iband from 1 to nFrequencyBandsForCoGcomputation
								appendFile: results_file$, tab$, median_spectral_dispersion_band'iband'
							endfor
						endif
					endif
					if get_standard_deviation
						appendFile: results_file$, tab$, std_spectral_dispersion
						if extract_energy_frequency_bands
							for iband from 1 to nFrequencyBandsForCoGcomputation
								appendFile: results_file$, tab$, std_spectral_dispersion_band'iband'
							endfor
						endif
					endif
					if get_min_max_with_time
						appendFile: results_file$, tab$, min_spectral_dispersion, tab$, min_spectral_dispersion_relative_time, tab$, max_spectral_dispersion, tab$, max_spectral_dispersion_relative_time
						if extract_energy_frequency_bands
							for iband from 1 to nFrequencyBandsForCoGcomputation
								appendFile: results_file$, tab$, min_spectral_dispersion_band'iband'
								appendFile: results_file$, tab$, min_spectral_dispersion_relative_time_band'iband'
								appendFile: results_file$, tab$, max_spectral_dispersion_band'iband'
								appendFile: results_file$, tab$, max_spectral_dispersion_relative_time_band'iband'
							endfor
						endif
					endif
				endif
				# spectral tilt
				if extract_spectral_tilt
					selectObject: spectral_tilt_values_extracted_points
					for ipoint from 1 to n_extracted_points
						current_point_value = Get value: ipoint, 1
						appendFile: results_file$, tab$, current_point_value
					endfor
					if get_mean
						appendFile: results_file$, tab$, mean_spectral_tilt
					endif
					if get_median
						appendFile: results_file$, tab$, median_spectral_tilt
					endif
					if get_standard_deviation
						appendFile: results_file$, tab$, std_spectral_tilt
					endif
					if get_min_max_with_time
						appendFile: results_file$, tab$, min_spectral_tilt, tab$, min_spectral_tilt_relative_time, tab$, max_spectral_tilt, tab$, max_spectral_tilt_relative_time
					endif
				endif
				# HNR
				if extract_HNR
					selectObject: hnr_values_extracted_points
					for ipoint from 1 to n_extracted_points
						current_point_value = Get value: ipoint, 1
						appendFile: results_file$, tab$, current_point_value
					endfor
					if get_mean
						appendFile: results_file$, tab$, mean_HNR
					endif
					if get_median
						appendFile: results_file$, tab$, median_HNR
					endif
					if get_standard_deviation
						appendFile: results_file$, tab$, std_HNR
					endif
					if get_min_max_with_time
						appendFile: results_file$, tab$, min_HNR, tab$, min_HNR_relative_time, tab$, max_HNR, tab$, max_HNR_relative_time
					endif
				endif
				# ZCR
				if extract_ZCR
					selectObject: zcr_values_extracted_points
					for ipoint from 1 to n_extracted_points
						current_point_value = Get value: ipoint, 1
						appendFile: results_file$, tab$, current_point_value
					endfor
					if get_mean
						appendFile: results_file$, tab$, mean_ZCR
					endif
					if get_median
						appendFile: results_file$, tab$, median_ZCR
					endif
					if get_standard_deviation
						appendFile: results_file$, tab$, std_ZCR
					endif
					if get_min_max_with_time
						appendFile: results_file$, tab$, min_ZCR, tab$, min_ZCR_relative_time, tab$, max_ZCR, tab$, max_ZCR_relative_time
					endif
				endif
				# formants
				if extract_formants
					selectObject: f1_values_extracted_points
					for ipoint from 1 to n_extracted_points
						current_point_value = Get value: ipoint, 1
						appendFile: results_file$, tab$, current_point_value
					endfor
					if get_mean
						appendFile: results_file$, tab$, mean_f1
					endif
					if get_median
						appendFile: results_file$, tab$, median_f1
					endif
					if get_standard_deviation
						appendFile: results_file$, tab$, std_f1
					endif
					if get_min_max_with_time
						appendFile: results_file$, tab$, min_f1, tab$, min_f1_relative_time, tab$, max_f1, tab$, max_f1_relative_time
					endif

					selectObject: f2_values_extracted_points
					for ipoint from 1 to n_extracted_points
						current_point_value = Get value: ipoint, 1
						appendFile: results_file$, tab$, current_point_value
					endfor
					if get_mean
						appendFile: results_file$, tab$, mean_f2
					endif
					if get_median
						appendFile: results_file$, tab$, median_f2
					endif
					if get_standard_deviation
						appendFile: results_file$, tab$, std_f2
					endif
					if get_min_max_with_time
						appendFile: results_file$, tab$, min_f2, tab$, min_f2_relative_time, tab$, max_f2, tab$, max_f2_relative_time
					endif

					selectObject: f3_values_extracted_points
					for ipoint from 1 to n_extracted_points
						current_point_value = Get value: ipoint, 1
						appendFile: results_file$, tab$, current_point_value
					endfor
					if get_mean
						appendFile: results_file$, tab$, mean_f3
					endif
					if get_median
						appendFile: results_file$, tab$, median_f3
					endif
					if get_standard_deviation
						appendFile: results_file$, tab$, std_f3
					endif
					if get_min_max_with_time
						appendFile: results_file$, tab$, min_f3, tab$, min_f3_relative_time, tab$, max_f3, tab$, max_f3_relative_time
					endif

					selectObject: f4_values_extracted_points
					for ipoint from 1 to n_extracted_points
						current_point_value = Get value: ipoint, 1
						appendFile: results_file$, tab$, current_point_value
					endfor
					if get_mean
						appendFile: results_file$, tab$, mean_f4
					endif
					if get_median
						appendFile: results_file$, tab$, median_f4
					endif
					if get_standard_deviation
						appendFile: results_file$, tab$, std_f4
					endif
					if get_min_max_with_time
						appendFile: results_file$, tab$, min_f4, tab$, min_f4_relative_time, tab$, max_f4, tab$, max_f4_relative_time
					endif
				endif
				# MFCC
				if extract_MFCC
					selectObject: mfcc_values_extracted_points
					for ipoint from 1 to n_extracted_points
						for iMFCC from 0 to numberOfMFCCcoefficients
							current_point_value = Get value: ipoint, iMFCC+1
							appendFile: results_file$, tab$, current_point_value
						endfor
					endfor
					if get_mean
						for iMFCC from 0 to numberOfMFCCcoefficients
							appendFile: results_file$, tab$, mean_MFCC'iMFCC'
						endfor
					endif
					if get_median
						for iMFCC from 0 to numberOfMFCCcoefficients
							appendFile: results_file$, tab$, median_MFCC'iMFCC'
						endfor
					endif
					if get_standard_deviation
						for iMFCC from 0 to numberOfMFCCcoefficients
							appendFile: results_file$, tab$, std_MFCC'iMFCC'
						endfor
					endif
					if get_min_max_with_time
						for iMFCC from 0 to numberOfMFCCcoefficients
							appendFile: results_file$, tab$, min_MFCC'iMFCC'
							appendFile: results_file$, tab$, min_MFCC'iMFCC'_relative_time
							appendFile: results_file$, tab$, max_MFCC'iMFCC'
							appendFile: results_file$, tab$, max_MFCC'iMFCC'_relative_time
						endfor
					endif
				endif
				# CPP and CPPS
				if extract_CPP
					selectObject: cpp_values_extracted_points
					for ipoint from 1 to n_extracted_points
						current_point_value = Get value: ipoint, 1
						appendFile: results_file$, tab$, current_point_value
					endfor
					if get_mean
						appendFile: results_file$, tab$, mean_CPP
					endif
					if get_median
						appendFile: results_file$, tab$, median_CPP
					endif
					if get_standard_deviation
						appendFile: results_file$, tab$, std_CPP
					endif
					if get_min_max_with_time
						appendFile: results_file$, tab$, min_CPP, tab$, min_CPP_relative_time, tab$, max_CPP, tab$, max_CPP_relative_time
					endif
					selectObject: cpps_values_extracted_points
					for ipoint from 1 to n_extracted_points
						current_point_value = Get value: ipoint, 1
						appendFile: results_file$, tab$, current_point_value
					endfor
					if get_mean
						appendFile: results_file$, tab$, mean_CPPS
					endif
					if get_median
						appendFile: results_file$, tab$, median_CPPS
					endif
					if get_standard_deviation
						appendFile: results_file$, tab$, std_CPPS
					endif
					if get_min_max_with_time
						appendFile: results_file$, tab$, min_CPPS, tab$, min_CPPS_relative_time, tab$, max_CPPS, tab$, max_CPPS_relative_time
					endif
				endif
				# Raw signal values
				if extract_signal_value
					selectObject: signal_value_extracted_points
					for ipoint from 1 to n_extracted_points
						current_point_value = Get value: ipoint, 1
						appendFile: results_file$, tab$, current_point_value
					endfor
					if get_mean
						appendFile: results_file$, tab$, mean_signal_value
					endif
					if get_median
						appendFile: results_file$, tab$, median_signal_value
					endif
					if get_standard_deviation
						appendFile: results_file$, tab$, std_signal_value
					endif
					if get_min_max_with_time
						appendFile: results_file$, tab$, min_signal_value, tab$, min_signal_value_relative_time, tab$, max_signal_value, tab$, max_signal_value_relative_time
					endif
				endif

				# Clean-up: remove temporary storage of values extracted on each point
				if extract_F0
					removeObject: f0_values_extracted_points
				endif
				if extract_intensity
					removeObject: intensity_values_extracted_points
				endif
				if extract_energy_frequency_bands
					for iband from 1 to nFrequencyBandsForCoGcomputation
						removeObject: energy_extracted_points_band'iband'
					endfor
				endif
				if extract_CoG
					removeObject: cog_values_extracted_points
					if extract_energy_frequency_bands
						for iband from 1 to nFrequencyBandsForCoGcomputation
							removeObject: cog_values_extracted_points_band'iband'
						endfor
					endif
				endif
				if extract_spectral_dispersion
					removeObject: spectral_dispersion_values_extracted_points
					if extract_energy_frequency_bands
						for iband from 1 to nFrequencyBandsForCoGcomputation
							removeObject: spectral_dispersion_values_extracted_points_band'iband'
						endfor
					endif
				endif
				if extract_spectral_tilt
					removeObject: spectral_tilt_values_extracted_points
				endif
				if extract_HNR
					removeObject: hnr_values_extracted_points
				endif
				if extract_ZCR
					removeObject: zcr_values_extracted_points
				endif
				if extract_formants
					removeObject: f1_values_extracted_points, f2_values_extracted_points, f3_values_extracted_points, f4_values_extracted_points
				endif
				if extract_MFCC
					removeObject: mfcc_values_extracted_points
				endif
				if extract_CPP
					removeObject: cpp_values_extracted_points, cpps_values_extracted_points
				endif
				if extract_signal_value
					removeObject: signal_value_extracted_points
				endif

				# Extract labels from other tiers and append information to the results file
				# Get interval midpoint (used as reference to extract information from other tiers)
				mid_point = start_time + duration/2
				if secondary_tier>0
					# Get the corresponding label on the selected secondary tier
					selectObject: current_tg
					intervtmp = Get interval at time: secondary_tier, mid_point
					secondaryTierlabel$ = Get label of interval: secondary_tier, intervtmp
					# Search and replace in label if requested
					if searchAndReplaceInLabels
						secondaryTierlabel$ = replace_regex$(secondaryTierlabel$, regexToReplaceInLabels$, replacementStringInLabels$, 0)
					endif
					secondaryTierStartTime = Get start time of interval: secondary_tier, intervtmp
					secondaryTierEndTime = Get end time of interval: secondary_tier, intervtmp
					appendFile: results_file$, tab$, secondaryTierlabel$, tab$, secondaryTierStartTime, tab$, secondaryTierEndTime
					if extract_left_and_right_context
						if intervtmp-1 > 0
							previousLabelSecondaryTier$ = Get label of interval: secondary_tier, intervtmp-1
							# Search and replace in label if requested
							if searchAndReplaceInLabels
								previousLabelSecondaryTier$ = replace_regex$(previousLabelSecondaryTier$, regexToReplaceInLabels$, replacementStringInLabels$, 0)
							endif
						else
							previousLabelSecondaryTier$ = "--undefined--"
						endif
						nIntervalsSecondaryTier = Get number of intervals: secondary_tier
						if intervtmp+1 <= nIntervalsSecondaryTier
							followingLabelSecondaryTier$ = Get label of interval: secondary_tier, intervtmp+1
							# Search and replace in label if requested
							if searchAndReplaceInLabels
								followingLabelSecondaryTier$ = replace_regex$(followingLabelSecondaryTier$, regexToReplaceInLabels$, replacementStringInLabels$, 0)
							endif
						else
							followingLabelSecondaryTier$ = "--undefined--"
						endif
						appendFile: results_file$, tab$, previousLabelSecondaryTier$, tab$, followingLabelSecondaryTier$
					endif
				elsif secondary_tier = -1
					# Get the corresponding labels on all interval tiers
					# Loop every tier
					for itier from 1 to ntiers
						# Ignore it if it's the reference tier (already processed) or a point tier (no labels to extract)
						selectObject: current_tg
						interv_tier = Is interval tier: itier
						if itier<>reference_tier and interv_tier=1
							selectObject: current_tg
							# Get label at reference tier's current interval midpoint and append it to results file
							intervtmp = Get interval at time: itier, mid_point
							secondaryTierlabel$ = Get label of interval: itier, intervtmp
							# Search and replace in label if requested
							if searchAndReplaceInLabels
								secondaryTierlabel$ = replace_regex$(secondaryTierlabel$, regexToReplaceInLabels$, replacementStringInLabels$, 0)
							endif
							secondaryTierStartTime = Get start time of interval: itier, intervtmp
							secondaryTierEndTime = Get end time of interval: itier, intervtmp
							appendFile: results_file$, tab$, secondaryTierlabel$, tab$, secondaryTierStartTime, tab$, secondaryTierEndTime
							if extract_left_and_right_context
								if intervtmp-1 > 0
									previousLabelSecondaryTier$ = Get label of interval: itier, intervtmp-1
									# Search and replace in label if requested
									if searchAndReplaceInLabels
										previousLabelSecondaryTier$ = replace_regex$(previousLabelSecondaryTier$, regexToReplaceInLabels$, replacementStringInLabels$, 0)
									endif
								else
									previousLabelSecondaryTier$ = "--undefined--"
								endif
								nIntervalsSecondaryTier = Get number of intervals: itier
								if intervtmp+1 <= nIntervalsSecondaryTier
									followingLabelSecondaryTier$ = Get label of interval: itier, intervtmp+1
									# Search and replace in label if requested
									if searchAndReplaceInLabels
										followingLabelSecondaryTier$ = replace_regex$(followingLabelSecondaryTier$, regexToReplaceInLabels$, replacementStringInLabels$, 0)
									endif
								else
									followingLabelSecondaryTier$ = "--undefined--"
								endif
								appendFile: results_file$, tab$, previousLabelSecondaryTier$, tab$, followingLabelSecondaryTier$
							endif
						endif
					endfor
				endif

				# Append a line break to the results file before proceeding to the next interval
				appendFile: results_file$, newline$
			endif
		endif
	endfor
	# Clean-up: remove current textgrid, pitch, intensity, formant and sound objects
	removeObject: current_tg
	if extract_F0+extract_intensity+extract_formants+extract_CoG+extract_spectral_dispersion+extract_spectral_tilt+extract_HNR+extract_ZCR+extract_MFCC+extract_signal_value>0
		removeObject: current_snd
	endif
endfor

appendInfoLine: newline$, "Processed ", ntextgrids, " files."

# Clean-up: remove lists of textgrids and vowels
removeObject: flist, extracted_points_relative_times
if filter_labels = 1
	removeObject: stringsDictionary
endif
if extract_energy_frequency_bands
	removeObject: frequencyBandsDefinitionForCoGcomputationTable 
endif
