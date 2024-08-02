# Praat_acoustic_analysis
Set of parameterizable Praat scripts for acoustic speech analysis on a set of matched sound and textgrid files.

The main script extracts selected parameters on every label matching a set defined in an external text file with their start and end time and duration. If no external text file is specified, all non-empty intervals are processed.
Aligned labels on other tiers can be optionally extracted (either one selected secondary tier or all others), as well as previous and next labels on each target tier.

The following acoustic parameters can be extracted on request on a variable number of user-defined points, plus mean value on the whole interval, min, max with time and standard deviation on request:
- F0
- Intensity
- Formants 1 to 4
- Spectral center of gravity (CoG), optionally limited to frequency bands, with energy in each band (may be used for the extraction of the level of energy in frequency bands)
- Spectral dispersion
- Spectral tilt (difference of energy between frequency bands)
- Harmonics-to-noise ratio (HNR)
- Zero-crossing rate (ZCR)
- Cepstral Peak Prominence (CPP) + smoothed version CPPS
- MFCC coefficients
- Signal value (may be useful when applied to physiological signals)

Selected features and target extraction points are defined in an external text file.
The script assumes that matched textgrid and sound files have the same name, and that all textgrids
have the same structure.
