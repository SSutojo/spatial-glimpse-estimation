function [ FEXT ] = featureextract_init( nChans, varargin )
%FEATUREEXTRACT_INIT initialize feature extractor
%   Sets up the structure to extract features

% set up a list of ERB spaced center frequencies as default for the
% One-Zero filterbank
HzToERBRate = @(Hz) 21.4*log10(Hz*0.00437 + 1.0);
ERBRateToHz = @(erbf) (10.^(erbf/21.4)-1)/4.37e-3;
erbLo = HzToERBRate(80);
erbHi = HzToERBRate(5000);       
fCenter(1)    = ERBRateToHz(erbLo);
fCenter(2:31) = ERBRateToHz(erbLo + (1:30) .* (erbHi-erbLo) / (32-1.0));
fCenter(32)   = ERBRateToHz(erbHi);

p = inputParser;
p.addRequired('nChans', @isscalar);
p.addParameter('fb', 'may', @isstr);
p.addParameter('fs', 16000, @isscalar);
p.addParameter('blocksize', 320, @isscalar);
p.addParameter('window', []);
p.addParameter('hopsize', 160, @isscalar);
p.addParameter('cfs', fCenter, @isvector);
p.addParameter('maxITDsec', 1.1e-3, @isscalar);
p.addParameter('maxTDOAsec', 1e-4, @isscalar);
p.addParameter('minXCreliable', 0.3, @isscalar);
p.addParameter('removeEndpoints', true, @islogical);
p.addParameter('HaircellModel', 'halfwave', @isstr);
p.parse(nChans, varargin{:});

FEXT = p.Results;
FEXT.nFilts = length(FEXT.cfs);

fCentErb = ((fCenter/9.26449) + 24.7);
B        = 2 * pi * fCentErb * 1.019;
FEXT.delays    = round(round(FEXT.fs * 3./B));

FEXT.maxDelay = ceil(p.Results.maxITDsec * FEXT.fs);
FEXT.lagsITD = floor(-p.Results.maxITDsec * FEXT.fs):1:FEXT.maxDelay;

switch p.Results.HaircellModel
    case 'none'
        FEXT.haircell = @(x) x;
    case 'halfwave'
        FEXT.haircell = @(x) max(x,0);
    case 'roman'
        FEXT.haircell = @(x) sqrt(max(x,0));
    otherwise % TBD: 'haircell' or 'envelope'
        error('Haircell model %s not implemented.', p.Results.NeuralModel);
end

switch p.Results.fb
    case 'may'
        for n=1:nChans
            FEXT.gtfb_state{n} = zeros(8, FEXT.nFilts);
        end
        FEXT.gtfb = @(x, n, STATE) MayGTFB(x, n, STATE);
    case 'OZGTFB'
        for n=1:nChans
            FEXT.gtfb_state{n} = zeros(2, 4, FEXT.nFilts);
        end
        FEXT.gtfb = @(x, n) OZGTFB(FEXT.fs, FEXT.cfs, x, FEXT.gtfb_state{n});
    case 'FFT'
        FEXT.gtfb = @(x, n) error('No GTFB, use FFT!');
        FEXT.fftsize = 2^nextpow2(FEXT.blocksize);
        % calculate the GTFB band to which the FFT bins belong
        F = (FEXT.fftsize/2)+1;
        fb = linspace(0, FEXT.fs/2, F);
        fb_idx = zeros(1,F);
        fb_ss = zeros(2,FEXT.nFilts);
        fCp = [0 fCenter FEXT.fs/2];
        for f=1:F
            [~, fb_idx(f)] = min(abs(fCp-fb(f)));
        end
        for f=1:FEXT.nFilts
            fb_ss(1,f) = find(fb_idx==(f+1),1);
            fb_ss(2,f) = find(fb_idx==(f+1),1,'last');
        end
        FEXT.fb_ss = fb_ss-1;
    otherwise
        error('FB %s not implemented.', p.Results.fb);
end

if isempty(p.Results.window)
    FEXT.blkwindow = ones(FEXT.blocksize, 1);
elseif isvector(p.Results.window) && (length(p.Results.window)==p.Results.blocksize)
    FEXT.blkwindow = p.Results.window;
else
    error('Invalid window specification.');
end

% initialize delay alignmenty memory
md = max(FEXT.delays);
for mm = 1:FEXT.nFilts
    FEXT.delaymem{mm} = zeros( md-FEXT.delays(mm), 1, nChans);
end

end

function [ out, FEXT ] = MayGTFB( x, n, FEXT )
    [out, FEXT.gtfb_state{n}] = gammatone2MEX(x, FEXT.fs, ...
        FEXT.cfs(1), FEXT.cfs(end), ...
        [FEXT.nFilts 1:FEXT.nFilts], 0, 1, 0, FEXT.gtfb_state{n}); 
end
