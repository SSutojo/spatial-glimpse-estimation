A spatial glimpse estimation algorithm

------------------------------------------------------------------------------------------

The spatial_glimpse_estimation algorithm divides a binaural mixture of multiple speakers 
into spectro-temporal glimpses that can be assigned to different positions in the azimuth 
plane. By combining glimpses that originate from the same azimuth position, a binary mask 
is estimated that can be used for source separation.

------------------------------------------------------------------------------------------
1. Processing an individual sentence mixture (example)

Download the HRTF database (filename: KEMAR.h5) from the following link:
https://zenodo.org/record/1226873#.YnzFu-jMLcs
Save the file KEMAR.h5 to the subdirectory \stimuli\HRTF\MMHR-HRTF

Add all paths by running the script "collect_paths".
To perform the source separation on an arbitrary mixed signal, use the script 
"process_individual_sentence.m". An example mixture (variable: "noisy_mix") and a set of 
suitable model parameters (struct: "params") are pre-defined in the script.
By running this script, you can make sure that all necessary functions are available or 
may need to be compiled for your system (see Sec. 4). If all necessary functions are 
available the script performs the source separation on the given sentence mixture and 
generates a number of plots, showing the estimated features, feature contrasts, binary 
masks, etc.

If you want to process a different signal mixture, replace the variable "noisy_mix" (in 
line 49) with your own binaural signal (dimensions: signal length X number of channels). 
If the separate source signals are also available, the variable "signal_components" (in 
line 49) can also be replaced (dimensions: signal length X number of channels X number of 
separate sources).
Potentially, certain parameters such as the sampling frequency 
have to be adapted when using a different signal (see Sec. 2. for parameter settings).

The separation algorithm is executed in spatial_glimpses_main and can be divided into 4 
main stages:
-	pre-processing (division into time frames and frequency subbands) yields a 
	time-frequency representation (T-F units) for the feature extraction
-	feature extraction (extraction of various features, such as azimuth or periodicity 
	estimates per individual T-F units)
-	local grouping (estimation of feature contrasts between neighboring T-F units and 
	based on this, formation of spectro-temporal glimpses as clusters of T-F units with 
	homogeneous features)
-	global grouping (identification of promiment spatial features in each glimpse and 
	assignment of spatial glimpses to previously detected source locations), yields 
	estimation of binary mask per source

The outputs of each separate stage can be found in the results structure "main_out".
------------------------------------------------------------------------------------------
2. parameter settings

The parameters are controlled with the def_..._params functions and with update_params_new 
which allows to update individual parameters according to table of new parameters, read 
from a …_settings_... .txt file (see directory: evaluation>specs_eval_condidions).

A complete set of all parameters is set using the def_..._params and init_..._params 
functions. While def_... sets independent parameters, init_... derives further parameters 
that are dependent on the settings in def_… (such as filterbank parameters).

The parameter structs resulting from init_… are named according to the individual stages 
of the segmentation algorithm (e.g. prep_params contains als parameters needed during 
pre-processing) and are summarized in a single struct called params.
Any parameter can be changed by replacing it directly in the def_..._params function or 
by using the update_params_new function. This function takes two cell arrays as input, 
one of which contains the new parameter values, while the other defines the exact 
fieldnames within the params struct which should be changed. All parameters that are not 
updated here, remain the same as defined by def_..._params.
Note, that all parameters which are dependent on other parameters should not be updated 
to avoid inconsistencies (see def_..._params functions for all independent parameters

------------------------------------------------------------------------------------------
3. available contrast estimation model

The best performing contrast estimation model according to the corresponding publication
is included in this package. Its parameter settings can be used with:
[eval_settings, var_names] = load_model_settings('model_settings'); 
eval_settings              = eval_settings(1,:);   

To use a single (untrained) contrast features (here azimuth correlation) use:
[eval_settings, var_names] = load_model_settings('model_settings'); 
eval_settings              = eval_settings(2,:);   

The included contrast estimation model and localizer are trained with the delivered HRTF 
database.

------------------------------------------------------------------------------------------
4. compile -mex functions

To compile missing –mex functions (e.g. xcorrNormMEX) for your system, use the command:
mex –largeArrayDims xcorrNormMEX.cpp
(the new function should appear in the current working directory)

------------------------------------------------------------------------------------------
5. Reference

Developed with Matlab 2018b

For citing this data, please use "Segmentation of Multitalker Mixtures Based on Local 
Feature Contrasts and Auditory Glimpses", by S. Sutojo, T. May and S. van de Par, IEEE 
Transactions on Audio, Speech, and Language Processing, vol.30, 2022.
https://doi.org/10.1109/TASLP.2022.3155285

This work includes external ressources:

HRTF database
 	J. Thiemann and S. van de Par,"A multiple model high-resolution head-related 
	impulse response database for aided and unaided ears", 
        EURASIP J. Adv. Signal Process., no. 9, pp, 1-9, 2019.
	link: https://zenodo.org/record/1226873#.YnzFu-jMLcs   
	(downlad KEMAR.h5 database) , retrieved May 12, 2022. 
	

auditory modeling toolbox
 	The pre-processing was performed with the invertible gammatone filterbank 
	(hohmann2002, V. Hohmann,"Frequency analysis and synthesis using a Gammatone 
	filterbank", Acta Acust. united with Acustica, vol. 88, no. 3, pp.433-442, 
	2002) implemented in the auditory modeling toolbox version 1.1.0 (Majdak, P., 
	Hollomey, C., and Baumgartner, R. "AMT 1.x: A toolbox for reproducible research 
	in auditory modeling," Acta Acust. vol. 6, no.19, 2022)
        link: https://amtoolbox.org/download.php , retrieved May 12, 2022.
 	
Netlab - Pattern analysis toolbox
        Ian Nabney (2022)
        link: https://www.mathworks.com/matlabcentral/fileexchange/2654-netlab), 
	MATLAB Central File Exchange. Retrieved May 12, 2022.


You can also use the "noisy_mix" function to create signal mixtures using the 
externatal speech database TIMIT acoustic-phonetic continuous speech corpus
	J. S. Garofolo et al., “TIMIT acoustic-phonetic continuous speech corpus LDC93S1,” Web Download, Philadelphia, PA, USA, Linguistic Data Consortium, 1993.
        link: https://catalog.ldc.upenn.edu/LDC93s1

To generate sepakers mixtures with the function 'noisy_mix', set the input variable 
'mode' to either 'training' or 'testing' (depending which part of the speech corpus 
should be used) and define the speaker IDs and sentence IDs with the input variable 
'mixIDs' (see function for further details and example).

