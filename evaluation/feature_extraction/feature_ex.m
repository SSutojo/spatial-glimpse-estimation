function [feature_mats] = feature_ex(pre_processed ,feature_types, feature_params, prep_params)
%MAIN FUNCTION FOR FEATURE EXTRACTION
% 
% INPUTS:
%      feature_types - defines which types of features should be extracted (ITD,ILD,pitch)
%      decomp_sig_seg - full signal decomposed into frequency bands and
%                       time segments (DIMS: blocksize X num_win X filter_bands X earChans)
%      feature_params, prep_params - structure with all the parameters
%
% OUTPUTS:
%      feature_mats   - a structure that contains a field/matrix for each
%                       extracted feature type
%
% NOTE
% the shape of the feature matrices should be (N x F) with N=number of time
% frames and F=number of filterbands (for ILDs and ITDs)
% 
% hopsize needs to be consistent across all features and ground truth, that
% is why it is taken from prep-Params

% get some parameters
decomp_sig_seg              = pre_processed.decomp_sig_seg;             % outputs of the FB segmented into time frames -> no haircell model
Ldecomp_sig_seg             = squeeze(decomp_sig_seg(:,:,:,1));         % left ear

num_win                     = size(Ldecomp_sig_seg,2);                  % number of time frames
bands                       = size(Ldecomp_sig_seg,3);                  % number of frequency bands 



num_features                = length(feature_types);                    % number of feature types that should be extracted
for aa = 1:num_features
    
    switch upper(feature_types{aa})
        case 'PERIODICITY'                                              % calculate normalized auto correlation for each t-f-unit to estimate periodicity , f0
            % check if the blocksize is a whole-number multiple of the hopsize and determine the needed amount of zero-padding at start and end
            if mod(feature_params.pitch_blocksize,prep_params.hopsize) ~= 0 
                error('all blocksizes must be whole-number multiples of the hopsize')
            else
                zero_padding        = zeros(((feature_params.pitch_blocksize/2) - prep_params.hopsize),bands,2); % for zero_padding at beginning and end of the filterbank outputs before time segmentation
            end  
            FEXT_pitch.hopsize        = prep_params.hopsize;
            FEXT_pitch.blocksize      = feature_params.pitch_blocksize;            
            
            % window, according to type and blocksize in samples          
            if strcmp(feature_params.pitch_wintype,'rect')==1    
                FEXT_pitch.blkwindow       = ones(feature_params.pitch_blocksize, 1);
            elseif strcmp(feature_params.pitch_wintype,'hann')==1
                FEXT_pitch.blkwindow       = hann(feature_params.pitch_blocksize);
            elseif strcmp(feature_params.pitch_wintype,'shan')==1
                FEXT_pitch.blkwindow       = sqrt(hann(feature_params.pitch_blocksize));
            else
                error('define a valid window type')
            end
            
            % determine if haircell model should be used before the periodicity estimation (also in pre_processed.hairc_dec_seg)        
            if strcmp(feature_params.model,'haircell_yes') == 1         
                fb_out              = cat(1, zero_padding, pre_processed.hairc_dec, zero_padding);               
                dec_seg             = featureextract_block(FEXT_pitch, fb_out);
            elseif strcmp(feature_params.model,'waveform')==1           % the direct output of the filterbank (finestructure, also in pre_processed.decomp_sig_seg)
                fb_out              = cat(1, zero_padding, pre_processed.decomp_sig, zero_padding); 
                dec_seg             = featureextract_block(FEXT_pitch, fb_out);          
            end
            dec_seg                 = squeeze(sum(dec_seg,4));              % sum over both ears and squeeze. dimension(window length X number of time frames X filterbands) 

            % calculate normalized autocorrelation NAC in each T-F unit (divide by NAC of the window)
            NAC_norm                = calcACorr(FEXT_pitch.blkwindow, feature_params.maxLags, feature_params.pitch_scale);
            complete_NAC_mat        = zeros((feature_params.maxLags-feature_params.minLags+1), num_win, bands );
            for kk = 1:bands                                            % over all gammatone filters
                for jj = 1:num_win                                      % over all time frames
                    t_f_unit                    = dec_seg(:,jj,kk);                 
                    [NAC_vec,~]                 = calcACorr(t_f_unit, feature_params.maxLags, feature_params.pitch_scale);
                    complete_NAC_mat(:,jj,kk)   = NAC_vec(feature_params.minLags:end) ./ NAC_norm(feature_params.minLags:end);
                end                 
            end
            feature_mats.complete_NAC_mat       = complete_NAC_mat;
            
            
        case 'AZIMUTH'                                                  % extract spatial cues, azimuth probabilities
            % check if the blocksize is a whole-number multiple of the hopsize and determine the needed amount of zero-padding at start and end
            if mod(feature_params.azm_blocksize,prep_params.hopsize) ~= 0 
                 error('all blocksizes must be whole-number multiples of the hopsize')
            else
                 zero_padding       = zeros(((feature_params.azm_blocksize/2) - prep_params.hopsize),bands,2); % for zero_padding at beginning and end of the filterbank outputs before time segmentation
            end     
            FEXT_azm.hopsize        = prep_params.hopsize;
            FEXT_azm.blocksize      = feature_params.azm_blocksize;
            
            % window, according to type and blocksize in samples 
            if strcmp(feature_params.azm_wintype,'rect')==1     
                FEXT_azm.blkwindow       = ones(feature_params.azm_blocksize, 1);
            elseif strcmp(feature_params.azm_wintype,'hann')==1
                FEXT_azm.blkwindow       = hann(feature_params.azm_blocksize);
            elseif strcmp(feature_params.azm_wintype,'shan')==1
                FEXT_azm.blkwindow       = sqrt(hann(feature_params.azm_blocksize));
            else
                error('define a valid window type')
            end
            
            % determine if haircell model should be used
            FEXT_azm.HaircellModel  = feature_params.HaircellModel;            
            FEXT_azm.lagsITD        = feature_params.lagsITD;
            FEXT_azm.maxDelay       = feature_params.maxDelay;
            FEXT_azm.fs             = prep_params.fs;        
            fb_out                  = cat(1, zero_padding, pre_processed.decomp_sig, zero_padding); 
            dec_seg                 = featureextract_block(FEXT_azm, fb_out);

            % get ITDs and ILDs            
            [ ILDs ]                = featureextract_ILD(FEXT_azm, dec_seg);           
            [ ITDs,xcorrs,~]        = featureextract_ITD(FEXT_azm, dec_seg);
            
            % use GMM classifier in struct C to map ITD/ILD combination to azimuth plane
            gmmCL                   = feature_params.C.classifier.gmmFinal;            
            nFrames                 = size(ILDs, 1);    % # of time frames
            nChannels               = length(gmmCL);    % # of frequency channels in pre-trained GMMs (should be consistent with feature matrix)
            nChannels_features      = size(ILDs,2);     % # of frequency channels in pre-processed mat/ feature mat (should by consistent with GMMS)
            
            % nr. of filterchannels in feature mat and pre-trained GMMS must be equal since GMM should be trained on the same t-f- unit size
            if nChannels ~= nChannels_features 
                warning('The numbers of filter channels in the feature matrix and in the pre-trained GMM are not equal')
            end
            
            nClasses                = size(gmmCL{1},2); % each class is an azimuth position
            az_probs1               = zeros(nClasses, nFrames, nChannels);
            az_probs                = zeros(nClasses, nFrames, nChannels);
            log_az_probs            = zeros(nClasses, nFrames, nChannels);
            for ff = 1:nChannels
                fvect                       = [ ILDs(:, ff) ITDs(:, ff) ];
                for ii = 1:nClasses
                    az_probs1(ii, :, ff)    = gmmprob(gmmCL{ff}(ii), fvect);   
                end
                for jj = 1:nClasses
                    az_probs(jj,:,ff)       = az_probs1(jj,:,ff)./(sum(az_probs1(:,:,ff),1)+eps);
                    log_az_probs(jj, :, ff) = log(az_probs(jj, :, ff)+eps); % apply log to the az_probs and add a small number to avoid NaNs
                end
            end
            feature_mats.az_probs       = az_probs;
            feature_mats.log_az_probs   = log_az_probs;
            feature_mats.ILDs           = ILDs;
            feature_mats.ITDs           = ITDs;
            feature_mats.xcorrs         = xcorrs;
    end
end

%% power_mat (braucht man immer daher kein case)
if mod(feature_params.pow_blocksize,prep_params.hopsize) ~= 0
    error('all blocksizes must be whole-number multiples of the hopsize')
else
    zero_padding = zeros(((feature_params.pow_blocksize/2) - prep_params.hopsize),bands,2); % for zero_padding at beginning and end of the filterbank outputs before time segmentation
end

FEXT_pow.hopsize    = prep_params.hopsize;
FEXT_pow.blocksize  = feature_params.pow_blocksize;

% window, according to type and blocksize in samples
if strcmp(feature_params.pow_wintype,'rect')==1    
    FEXT_pow.blkwindow       = ones(feature_params.pow_blocksize, 1);
elseif strcmp(feature_params.pow_wintype,'hann')==1
    FEXT_pow.blkwindow       = hann(feature_params.pow_blocksize);
elseif strcmp(feature_params.pow_wintype,'shan')==1
    FEXT_pow.blkwindow       = sqrt(hann(feature_params.pow_blocksize));
else
    error('define a valid window type')
end

fb_out                          = cat(1, zero_padding, pre_processed.decomp_sig, zero_padding);
dec_seg                         = featureextract_block(FEXT_pow, fb_out); % segmentation into time frames
dec_seg_squared                 = (abs(dec_seg)).^2;      % calculate energy of the time signal in each filterband+time-frame tile
power_mat                       = squeeze(sum(dec_seg_squared,1));              % sum over whole window length

feature_mats.power_mat          = power_mat;
feature_mats.summed_power_mat   = squeeze(sum(power_mat,3));  % also sum over both ear channels
feature_mats.time_sig           = squeeze(sum(dec_seg,4)); 
feature_mats.fileIdx_ac         = [1 nFrames];
feature_mats.fileIdx_inb        = [1 nFrames-1];

end









