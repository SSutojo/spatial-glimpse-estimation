function [pre_processed] = pre_process(sig_mix, prep_params)
%APPLY FILTERBANK AND SEGMENTATION INTO TIME FRAMES 
%
% 1.    a gammatone filterbank is applied to the signal
% 2.    optionally, a haircell model or so can be applied here
% 3.    the signal is segmented into time frames
%
% INPUTS  : stereo mix      - time signal (2 ear Channels X full signal length) 
%           prep_params     - structure with needed parameters and FB
%
% OUTPUT  : decomp_sig_seg (blocksize X num_win X filter_bands X earChans)
%           power_mats (num_win X filterbands X earChans)
%
% left ear is the first layer
nChans = prep_params.nChans;

%***************************** gammatone filterbank **********************
switch prep_params.fb_type
    case 'am_t'                                     % the filterbank from the auditory modelling toolbox
        
        
        if size(sig_mix,1)~=nChans                    % check the dimensions of the input matrix and change them if necessary
            sig_mix = sig_mix.';
        end
        
        decomp_sig      = zeros(size(sig_mix,2),prep_params.nFilts,nChans);
        decomp_sig_comp = zeros(size(sig_mix,2),prep_params.nFilts,nChans);
        for ii = 1:nChans
            mix_ii                  = sig_mix(ii,:);% for channel ii
            %[decomp_sig_ii, ~ ]     = gfb_analyzer_process(prep_params.analyzer, mix_ii);
            [decomp_sig_ii, ~ ]     = hohmann2002_process(prep_params.analyzer, mix_ii.');
            decomp_sig(:,:,ii)      = real(decomp_sig_ii);
            decomp_sig_comp(:,:,ii) = decomp_sig_ii;
            
        end
        
        pre_processed.decomp_sig_comp        = decomp_sig_comp;
        
    case 'may'                                      % the filterbank from Tobias May's localizer
        
        x = genAudio(sig_mix, prep_params.fs);   % desired input is an audio object instead of a stereo signal
        [decomp_sig, ~] = featureextract_auditory(prep_params, x, false); % last input is 'false', means that the haircell model 
                                                                                    % isn't applied here (done later, see below)
                          
    otherwise 
        error('undefined filterbank option')
end

% testsig = squeeze(sum(decomp_sig,2));
% soundsc(testsig,prep_params.fs)

hairc_dec = prep_params.haircell(decomp_sig);       % max(afout, 0); application of haircell model to the filterbands
%***************************** segmentation into time frames *************
FEXT.hopsize        = prep_params.hopsize;
FEXT.blocksize      = prep_params.blocksize_min;
FEXT.blkwindow      = prep_params.blkwindow_min;
[ decomp_sig_seg ]  = featureextract_block( FEXT, decomp_sig );
[ hairc_dec_seg]    = featureextract_block( FEXT, hairc_dec );

%***************************** calculate summed powers ***************
decomp_seg_squared  = (abs(decomp_sig_seg)).^2;      % calculate energy of the time signal in each filterband+time-frame tile
power_mat           = sum(decomp_seg_squared,1);              % sum over whole window length
power_mat           = squeeze(power_mat);


pre_processed.decomp_sig        = decomp_sig;       % filterbank outputs, no time segmentation
pre_processed.hairc_dec         = hairc_dec;        % filterbank outputs, haircell model, no time segmentation

pre_processed.decomp_sig_seg    = decomp_sig_seg;   % finestructure (direct outputs of the filterbank)
pre_processed.hairc_dec_seg     = hairc_dec_seg;    % filterbank outputs after haircell model
pre_processed.power_mat         = power_mat;        % power in each time-frequency unit 
pre_processed.cfs               = prep_params.cfs;  % safe the vector with center frequencies of the filterbank to the decomposed signal
end




