function [noisyMix, sigComponents] = noisy_mix(numSpeakers,azmVec, SNR, prep_params, mixIDs, mode, varargin)
%NOISY_MIX: create a spatialized mixture of TIMIT speakers/sentences, 
% optional noise (different noise types available) and optional harmonic 
% tone complex (can also be amplitude modulated, see options below)
%
% INPUTS: 
%   numSpeakers     - number of speakers
%   azmVec          - azimuth positions of speakers
%   SNR             - SNR of noise in dB (noise power/power of individual speaker)
%   prep_params     - struct with parameter settings
%   mixIDs          - cell array with speaker and sentence IDs (filenames)
%                     define the speaker IDs first, then the corresponding sentence IDs
%                     example for 3 speakers: mixIDs{1,1}='FCMH1';mixIDs{1,2}='MDAW1';mixIDs{1,3}='FJSJ0';mixIDs{1,4}='SX143'	;mixIDs{1,5}='SX283'; mixIDs{1,6}='SI2114';
%   mode            - either 'training' or 'testing' (different datasets), or 'example'
%   varargin        - settings if harmonic tone complex is also added
% 
% OUTPUTS: 
%   noisyMix        - stereo signal (binaural depending on used HRIRs)
%   sigComponents   - matrix with stereo signals of individual components
%
%


%% load mono signals
if isempty(varargin)
    htc{1}    = false;      rawSigs = cell(numSpeakers,1);   endpts = cell(numSpeakers,1);
    room.name = 'anechoic';
    
elseif ~isempty(varargin)
    if size(varargin,2) == 1
        room = varargin{1,1};
        op.spat_mode = 'hrtf';
        op.hrtf_database = 'mitkemar.sofa';
        op.fs = 44100;
        SrcRecDist  = 1.7;
        htc{1}=false;  rawSigs = cell(numSpeakers,1);   endpts = cell(numSpeakers,1);
    elseif size(varargin,2) == 2
        
        room = varargin{1,1};
        op.spat_mode = 'hrtf';
        op.hrtf_database = 'mitkemar.sofa';
        op.fs = 44100;
        %SrcRecDist  = 1.7;
        SrcRecDist  = 1.5;
        htc = varargin{1,2};
        if htc{1}==true
            rawSigs = cell(numSpeakers+1,1);
            endpts = cell(numSpeakers+1,1);
            f0          = htc{2};
            partial1    = 1;
            partial2    = htc{3};
            Fs          = 16000;
            azmVec      = [azmVec ; htc{4}];
        else
            rawSigs = cell(numSpeakers,1);
            endpts = cell(numSpeakers,1);
        end
    end
end

if strcmp(mode,'training')==1
    for aa = 1:numSpeakers
        path11              = [pwd,filesep,'stimuli',filesep,'TIMIT',filesep,'TRAIN',filesep,'DR8',filesep,cell2mat(mixIDs(1,aa))];
        addpath(path11)
        if prep_params.CARL == true  % if this runs on the cluster, a slightly different wacReadTimit function is needed
            [sig1,~,endpts1]    = wavReadTimit_mac([ cell2mat(mixIDs(1,numSpeakers+aa)),'.WAV']);
        else
            [sig1,~,endpts1]    = wavReadTimit([ cell2mat(mixIDs(1,numSpeakers+aa)),'.WAV']);
        end
        sig1                = resample(sig1,prep_params.fs,16000);
        rawSigs{aa,1}       = sig1;
        endpts{aa,1}        = ceil(endpts1.*(prep_params.fs/16000));
        rmpath(path11)
    end
    
elseif strcmp(mode,'testing')==1
    for aa = 1:numSpeakers
        path11              = [pwd,filesep,'stimuli',filesep,'TIMIT',filesep,'TEST',filesep,'DR8',filesep,cell2mat(mixIDs(1,aa))];
        addpath(path11)
        if prep_params.CARL == true
            [sig1,~,endpts1]    = wavReadTimit_mac([ cell2mat(mixIDs(1,numSpeakers+aa)),'.WAV']);
        else
            [sig1,~,endpts1]    = wavReadTimit([ cell2mat(mixIDs(1,numSpeakers+aa)),'.WAV']);
        end
        sig1                = resample(sig1,prep_params.fs,16000);
        rawSigs{aa,1}      = sig1;
        endpts{aa,1}        = ceil(endpts1.*(prep_params.fs/16000));
        rmpath(path11)
    end
    
elseif strcmp(mode,'example')==1
    for aa = 1:numSpeakers
        rawSigs{aa,1} = audioread(['sent',num2str(aa),'.wav']);
        endpts{aa,1}  = [1, 1; length(rawSigs{aa,1}-1) , length(rawSigs{aa,1}-1)];
    end
    
else
    disp('define if the training or testing set should be used')
end
%% generate harmonic tone complex if desired
if htc{1}==true
    sig_length      = length(sig1)/Fs;
    [htc_time_sig, ~]        = harmonic_tone_complex(f0, partial1, partial2, sig_length, Fs);       
    if htc{5} == true
        % create amplitude modulation and apply to the htc time signal
        AM_freq         = htc{6};  % frequency of the amplitude modulation
        AM_t            = 1/Fs.*[1:length(htc_time_sig)];
        AM_vec          = sin((2*pi*AM_freq).*AM_t);
        htc_time_sig    = htc_time_sig .* AM_vec;
    end
    rawSigs{end,1}  = htc_time_sig;
    endpts{end,1}   = [1 ,1; length(htc_time_sig)-1, length(htc_time_sig)-1]; 
end

%% cut, normalize, and spatialize speech signals
numSources      = size(rawSigs,1);
cutSigs         = cell(numSources,1);
sigLengths      = zeros(numSources,1);
for cc = 1:numSources 
    rawSig          = rawSigs{cc,1};
    endpoints       = endpts{cc,1};
    cutSigs{cc,1}   = rawSig(endpoints(1,2):endpoints(end,1));          % cut out silent parts in beginning and end of the sentence
    sigLengths(cc)  = endpoints(end,1)-endpoints(1,2)+1;
end
cutSigs{1,1}        = cutSigs{1,1}./rms(cutSigs{1,1});          % this will be the reference signal to adjust the other speakers

spatSigs            = cell((numSources),1);
elVec               = zeros(numSources,1);                     % elevation is always zero
fs_string           = [num2str(prep_params.fs/1000),'k']; 

spatSigLengths      = zeros(numSources,1);
for dd = 1:numSources      % the hrirs are adapted to the timit sentences 16 kHz
    if strcmp(room.name,'anechoic')==1 % use no additional razr room simulation
        currIr          = loadHRIRnear('KEMAR.h5',azmVec(dd),elVec(dd),fs_string);% there should be an azimuth available for each speaker  currIr  
    elseif strcmp(room.name,'anechoic')==0 % use razr to simulate room and generate the impulse response
        [SrcPosX, SrcPosY]  = Azm2CartRazr(room.recpos(1), room.recpos(2), SrcRecDist, azmVec(dd));
        room.srcpos         = [SrcPosX, SrcPosY, room.recpos(3)];
        ir = razr(room, op);
        currIr = resample(ir.sig,16000,44100);  % razr works best with 44.1 kHz sampling rate. Downsample to match TIMIT sentences
    end

    monoSig         = cutSigs{dd,1};
    earL            = conv(monoSig, currIr(:,1),'full');            % spatialize   earL            = conv(monoSig, currIr(:,1),'same'); 
    earR            = conv(monoSig, currIr(:,2),'full'); 
    earL = earL(1:length(monoSig));  % cut to the length of the original mono Signal
    earR = earR(1:length(monoSig));
    spatSigs{dd,1}  = [earL(:) earR(:)];
    spatSigLengths(dd)  = size(earL,2);
end 


%% adjust speech signals to the same level, make them same length, add silence
[maxLength, ~]  = max(sigLengths);                  % find the longest and the shortest sentence
[minLength, ~]  = min(sigLengths);


silence         = zeros(prep_params.fs*0.02,2);     % silence         = zeros(prep_params.fs*0.1,2); or short fs*0.04

for aa = 2:numSources    % adjust SNR, match them to speaker 1 (speaker 1 is already rms normalized)
    refSig          = spatSigs{1,1}; % reference signal is speaker 1
    currSig         = spatSigs{aa,1};
    summed_ref      = sum(refSig(1:minLength,:),2);   % important: sum over both ears before calculating SNR
    summed_curr     = sum(currSig(1:minLength,:),2);
    speakerFact     = sqrt(sum(summed_ref.^2) ./ sum(summed_curr.^2));
    
    spatSigs{aa,1}  = currSig * speakerFact;
end

% adding zeros to the short signals for equal length
for aa = 1:numSpeakers    
    currSig         = spatSigs{aa,1};
    addLength       = maxLength - length(currSig);      % additional length
    paddedSig       = [currSig ; zeros(addLength,2)];   % zero pad every signal     
    spatSigs{aa,1}  = [silence; paddedSig ; silence];
end

%% generate noise     
L               = size(spatSigs{1,1},1);            % desired length of noise (in samples)
noisetype       = prep_params.noise;

switch  noisetype
    case 'diffuse'
        b           = 1;                                  % filter option
        noise_obj   = spat_genDiffNoise(prep_params, L, b);
    case 'diffuse_pink'              
        noise_obj   = spat_genDiffNoisePINK(prep_params, L); %  pink diffuse noise   
    case 'white'
        noise_obj   =  genAudio(randn(L,prep_params.nChans), prep_params.fs);
    case 'colored'
        if strcmp(prep_params.noise,'pink')
            beta = -1; % 1/f noise
        else
            beta = -2; %1/f^2 noise
        end
        % create frequency vector folded at center:
        freqs = [linspace(0,.5,floor(L/2)),fliplr(linspace(0,.5,ceil(L/2)))]';
        freqs = (freqs.^2).^(beta/2);
        freqs(freqs==inf) = 0; %to prevent inf errors
        freqs = [freqs freqs];
        phase_vals = rand(L,prep_params.nChans);
        % now apply an inverse fft:
        cNoise = real(ifft(sqrt(freqs).*(cos(2*pi*phase_vals)+1i*sin(2*pi*phase_vals))));        
        noise_obj =  genAudio(cNoise, prep_params.fs);
    case 'icra5'
        [icra5 icra_fs] = audioread( [ pwd,filesep,'stimuli',filesep,'icra5_0.25_fullscale.wav']);
        % *************downsample to 16000 kHz
        icra5_resampled  = resample(icra5,prep_params.fs,icra_fs);
        icra5_bit_nChannel = zeros(L,prep_params.nChans);
        for ch = 1:prep_params.nChans
            iStart  = randi( size(icra5_resampled,1)-L);   
            icra5_bit = icra5_resampled(iStart:iStart+L-1,1);
            icra5_bit_nChannel(:,ch) =  icra5_bit ;
        end
        noise_obj = genAudio(icra5_bit, prep_params.fs);
    otherwise
        error('unsupported noise type')
end

%% adjust SNR between speakers and noise   
speech_segment      = spatSigs{1,1}(length(silence):minLength,:); % use a non-silent segments for the SNR calculation 
noise_sig           = noise_obj.data;
noise_segment       = noise_sig(length(silence):minLength,:);

summed_speech_seg   = sum(speech_segment,2); % sum over both ears before caluclating the SNR, that makes a bit of a difference
summed_noise_seg    = sum(noise_segment,2);
current_snr         = sqrt(sum(summed_speech_seg.^2) ./ sum(summed_noise_seg.^2));

snrFactor           = 10^(-SNR/10);    % inverse value of the factor that needs to be multiplied with the current SNR to yield the desired SNR in dB
noise_fact          = current_snr * snrFactor;

noise_sig                   = noise_sig.* noise_fact;
spatSigs{(numSources+1),1}  = noise_sig;
%% arrange outputs

noisyMix = zeros(L,2);
sigComponents  = zeros(L,2,(numSources+1));
for aa = 1:(numSources+1)    
    noisyMix                = noisyMix + spatSigs{aa,1};            % add up all sources
    sigComponents(:,:,aa)   = spatSigs{aa,1};                       % save all sources separately
end




