function [ground_truth] = pre_process_ground_truth(sig_components, params)
%APPLY FILTERBANK AND SEGMENTATION INTO TIME FRAMES TO THE SEPARATE SIGNAL
%COMPONENTS TO GET THE IDEAL BINARY MASKS ETC.
% 
% INPUTS:   sig_components:     time signals of the separate speakers/sources
%                               (full sig.length X 2 earChans X number of components)
%           prep_params:        structure with parameters
%           opts{1}:            defines the type of binary mask that is
%                               applied 
%           opts{2}:            the binary mask that is supplied (number of time frames X number of filterbands)
% 
% OUTPUTS:  ground_truth. ...
%           ... t_f_units       for each source: output of filterbank, segmented into time (#of samples in time frame X #time frames X #filterbands)
%           ... power_mats      for each source: power per t-f-unit (#time frames X #filterbands)
%           ... num_components  number of components/sources
%           ... sig_components  separate time (stereo) signals of the sources
%           ... ideal_source_map|map in which each t-f-unit is labeled with the ID of the dominant source (#time frames X #filterbands)
%           ... IBMs            the ideal binary mask for each separate component (#time frames X #filterbands X #sources)
%           ... ideal_contours  contours around the glimpses/ local clusters (2*#timeframes-1 X 2*#filterbands-1)
%           ... local_SNRs      the local SNRs in each t-f-unit for all separate sources (#time frames X #filterbands X #sources)

num_components = size(sig_components,3);    % get the number of components and save it to the ground truth struct

%% pre-processing for each components separately
for ii = 1:num_components                   
    curr_sig = sig_components(:,:,ii); 
    
    switch params.fb_type                   
        case 'am_t'                                                             % the filterbank from the auditory modelling toolbox
                           
            if size(curr_sig,1)~=2                                              % check the dimensions of the input matrix and change them if necessary
                curr_sig = curr_sig.';
            end               
            earL_mix = curr_sig(1,:);
            earR_mix = curr_sig(2,:);
            analyzer = hohmann2002(params.fs, params.f_low, params.base_frequency, params.f_high, params.filters_per_ERB);
            [Ldecomp_sig, analyzer1] = hohmann2002_process(analyzer, earL_mix.');% apply filterbank (decomposition into subbands)
            [Rdecomp_sig, analyzer2] = hohmann2002_process(analyzer, earR_mix.'); 
            Ldecomp_sig = real(Ldecomp_sig);                                  % important for further processing: only use real part of filter outputs
            Rdecomp_sig = real(Rdecomp_sig);     
            decomp_curr_sig = cat(3,Ldecomp_sig,Rdecomp_sig);                   % input shape for featureextract_block should be (full sig.length X bands X earChans)
                
        case 'may'                                                              % filterbank from the may localizer             
            x = genAudio(curr_sig, params.fs);                             % generate an audio object
            [decomp_curr_sig, ~] = featureextract_auditory(params, x, false);          % the second input should be an audio object 
                
        otherwise
            error('undefined filterbank option')
    end
    FEXT.hopsize        = params.hopsize;
    FEXT.blocksize      = params.blocksize_min;
    FEXT.blkwindow      = params.blkwindow_min;
    [ curr_decomp_seg ] = featureextract_block( FEXT, decomp_curr_sig ); % segmentation into time frames
        
    % calculate the power matrix for each component
    curr_squared = (abs(curr_decomp_seg)).^2;                                   % square of the magnitude
    curr_power_mat      = squeeze(sum(curr_squared,1));                                       % sum over the whole window length

    ground_truth.t_f_units{ii}  = curr_decomp_seg;                              
    ground_truth.power_mats{ii} = curr_power_mat;
    ground_truth.decomp_sig{ii} = decomp_curr_sig;
end


%% sum each sources' power_mat over both ear Channels and combine in one matrix
power_mat_sep = zeros(size(curr_power_mat,1),size(curr_power_mat,2),num_components); % #time frames X #Filterbands X #components/speakers

for jj = 1:num_components
    power_mat_sep(:,:,jj) = sum(ground_truth.power_mats{1,jj},3);               % power mat for the separate components (but summed over both ear channels)
end


%% get ideal mask

[ideal_source_map, IBMs, local_SNRs]              = get_ideal_mask(power_mat_sep, 'n', 1);   % the bin mask option is removed

nFilts  = size(local_SNRs,2);   % filterbands
nFrames = size(local_SNRs,1);   % time frames

%% ideal contrast values

sigm_SNRs = 1./(1+exp(-local_SNRs)); % use sigmoid func to derive sth like a source presence probability from the local SNRs

source_probs = sigm_SNRs./sum(sigm_SNRs,3); % normalize to get a probability

ideal_contrasts = ones((nFrames*2-1),(nFilts*2-1));

for aa = 1:nFrames-1        % within filterband transitions
    for bb = 1:nFilts
        ideal_contrasts(aa*2,bb*2-1)= sum(source_probs(aa,bb,:).*source_probs(aa+1,bb,:)); 
    end
end
for aa = 1:nFrames          % across filterband transitions
    for bb = 1:nFilts-1
        ideal_contrasts(aa*2-1,bb*2)= sum(source_probs(aa,bb,:).*source_probs(aa,bb+1,:));
    end
end

ideal_contrasts = 1 - ideal_contrasts;

%ideal_contrasts    = bilinear_interpolation(ideal_contrasts, 'favor_high_contrasts','vertices_only');
ideal_contrasts    = bilinear_interpolation(ideal_contrasts, 'favor_high_contrasts','vertices_and_TFunits');

%% ideal contours

ideal_contours = zeros(size(ideal_contrasts));
for bb = 2:2:size(ideal_contrasts,1)
    for cc = 1:2:size(ideal_contrasts,2)
        ideal_contours(bb,cc) = (ideal_contrasts(bb,cc)> 0.5);
    end
end
for dd = 1:2:size(ideal_contrasts,1)
    for ee = 2:2:size(ideal_contrasts,2)
        ideal_contours(dd,ee) = (ideal_contrasts(dd,ee)> 0.5);
    end
end
ideal_contours = bilinear_interpolation(ideal_contours,'favor_high_contrasts','vertices_only');
% or :
%[ideal_contours]                = get_contours(ideal_source_map);
%% ideal glimpses (based on ideal contours)

[regions, ~,~,~]    = regiongrow(ideal_contours, 0, 0.2);
ideal_glimpses      = zeros(nFrames, nFilts);
for bb = 1:nFrames 
    for cc = 1:nFilts
        ideal_glimpses(bb,cc)=regions((2*bb)-1,(2*cc)-1);
    end
end


%% 
ground_truth.ideal_contrasts    = ideal_contrasts;
ground_truth.num_components     = num_components;
ground_truth.sig_components     = sig_components;
ground_truth.ideal_source_map   = ideal_source_map;       % labels each t-f-unit with an integer that represents the locally dominant source
ground_truth.IBMs               = IBMs;
ground_truth.ideal_contours     = ideal_contours;
ground_truth.local_SNRs         = local_SNRs;
ground_truth.ideal_glimpses     = ideal_glimpses;    % after filling the ideal contours 
ground_truth.power_mat_sep      = power_mat_sep;


