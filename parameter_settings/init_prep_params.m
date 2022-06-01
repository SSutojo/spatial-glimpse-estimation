function [prep_params] = init_prep_params(p_prep)
%INIT_PREP_PARAMS initialize the parameter struct for pre-processing, based
%on the independent parameters (set in def_prep_params), further (dependent) parameters are derived 


prep_params = p_prep;

if prep_params.fs == 16000
    if strcmp(prep_params.HRIR_type,'KEMAR')==1
        prep_params.loadHRIRfun = @(az, el) loadHRIR('KEMAR.h5', az, el, '16k' ); 
    else
        error('define HRIR load function')
    end
else
    error('define HRIR load function')
end

%% temp segmentation
prep_params.hopsize             = prep_params.hopsize_ms * prep_params.fs * 1e-3; % in samples
prep_params.blocksize_min       = prep_params.Lwin_min * prep_params.fs * 1e-3;        
if strcmp(prep_params.window_type,'rect')
    prep_params.blkwindow_min       = ones(prep_params.blocksize_min, 1);  
elseif strcmp(prep_params.window_type,'hann')
    prep_params.blkwindow_min       = hann(prep_params.blocksize_min); 
elseif strcmp(prep_params.window_type,'shan') % sqrt_hann
    prep_params.blkwindow_min       = sqrt(hann(prep_params.blocksize_min)); 
else
    error('define a valid window type')
end 

%% filterbank
prep_params.erbstep             = 1/prep_params.filters_per_ERB;

switch prep_params.fb_type
    case 'am_t'
        prep_params.desired_delay_sec   = 0.01;  % before 0.004
        prep_params.base_frequency      = prep_params.f_low+1; % not used later, just defined to avoid errors        
        prep_params.cfs                 = ERB_center_freqs(prep_params.f_low, prep_params.f_high, 1/prep_params.filters_per_ERB); % center frequencies in Hz
        prep_params.bw                  = ((prep_params.cfs./9.26449) + 24.7); % bandwidth in dependence of the center frequency in Hz
        prep_params.start_freqs         = prep_params.cfs      - prep_params.bw.*0.5;
        prep_params.start_freqs(end+1)  = prep_params.cfs(end) + prep_params.bw(end).*0.5; % append the highest frequency (as end freq)
        prep_params.nFilts              = length(prep_params.cfs); % number of gammatone filters
        
       % prep_params.analyzer            = gfb_analyzer_new(prep_params.fs, prep_params.f_low, prep_params.base_frequency, prep_params.f_high, prep_params.filters_per_ERB); % create analyzer struct              
        prep_params.analyzer            = hohmann2002(prep_params.fs, prep_params.f_low, prep_params.base_frequency, prep_params.f_high, prep_params.filters_per_ERB); % create analyzer struct              
        
    case 'may'
        HzToERBRate                     = @(Hz) 21.4*log10(Hz.*0.00437 + 1); 
        ERBRateToHz                     = @(erbf) (10.^(erbf/21.4)-1)/4.37e-3;
        erbLo                           = ceil(HzToERBRate(prep_params.f_low))+ (prep_params.erbstep/2); % the lowest filter should not extend below f_low  
        erbHi                           = floor(HzToERBRate(prep_params.f_high))- (prep_params.erbstep/2); % the highest erb
%         erbLo                           = ceil(erbrate(prep_params.f_low))+ (prep_params.erbstep/2); % the lowest filter should not extend below f_low  
%         erbHi                           = floor(erbrate(prep_params.f_high))- (prep_params.erbstep/2); % the highest erb
        prep_params.cfs                 = ERBRateToHz(erbLo : prep_params.erbstep : erbHi); % vector with ERB spaced center frequencies for filterbands
%         prep_params.cfs                 = inverbrate(erbLo : prep_params.erbstep : erbHi); % vector with ERB spaced center frequencies for filterbands
        
        prep_params.bw                  = ((prep_params.cfs/9.26449) + 24.7); % calculate bandwidth in dependence of the center frequency in Hz (not kHz)
        B                               = 2 * pi * prep_params.bw * 1.019;
        prep_params.delays              = round(round(prep_params.fs * 3./B));
        prep_params.start_freqs         = prep_params.cfs - prep_params.bw.*0.5;
        prep_params.start_freqs(end+1)  = prep_params.cfs(end) + prep_params.bw(end).*0.5; % append the highest frequency (as end freq)
        
        prep_params.nFilts              = length(prep_params.cfs);
        for n = 1:nChans
            prep_params.gtfb_state{n}   = zeros(8, prep_params.nFilts); % maybe the 8 should be replaced. Does this determine the number of channels? (first input was an 8)
        end
        prep_params.gtfb                = @(x,n,STATE) MayGTFB(x,n,STATE); 
               
        md = max(prep_params.delays);       % initialize delay alignmenty memory
        for mm = 1:prep_params.nFilts
            prep_params.delaymem{mm} = zeros( md-prep_params.delays(mm), 1, nChans);
        end
        
    otherwise 
        error('undefined filterbank')
end 

%% haircell model

switch prep_params.haircell_type
    case 'none'
        prep_params.haircell            = @(x) x; % no haircell model is applied
    case 'halfwave_smooth'        
        [prep_params.C,prep_params.D]   = butter(prep_params.smooth_order,prep_params.smooth_freq/(0.5*prep_params.fs),'low'); % design a lowpass butterworth filter
        prep_params.haircell            = @(x) filter(prep_params.C,prep_params.D,max(x,0)); % halfwave rectification and smoothing with a low-pass filter
    case 'halfwave'
        prep_params.haircell            = @(x) max(x,0); % halfwave rectification
    otherwise 
        error('unknown haircell model')
end 



function [out, prep_params] = MayGTFB(x, n, prep_params) % replaces the old Filterbank
    [out, prep_params.gtfb_state{n}] = gammatone2MEX(x, prep_params.fs, ...
        prep_params.cfs(1), prep_params.cfs(end), ...
        [prep_params.nFilts 1:prep_params.nFilts], 0, 1, 0, prep_params.gtfb_state{n});    
end

function hrir = loadHRIR(varargin)
    hrir = loadHRIRnear(varargin{:});
end




end

