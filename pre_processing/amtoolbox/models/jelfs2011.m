function [benefit, weighted_SNR, weighted_bmld] = jelfs2011(target,interferer,varargin)
%JELFS2011  Binaural speech advantage
%
%   Usage:  [benefit weighted_SNR weighted_bmld] = jelfs2011(target,interferer,fs)
%  
%   Input parameters:
%     target           : Binaural target impulse response (or stimulus)
%     interfererer     : Binaural interferer impulse response (or stimulus)
%                        Multiple interfering impulse responses MUST be
%                        concatenated, not added.
%
%   Output parameters:
%     benefit       : spatial release from masking (SRM)in dB
%     weighted_SNR  : component of SRM due to better-ear listening (dB)
%     weighted_bmld : component of SRM due to binaural unmasking (dB)
%    
%   JELFS2011(target,interferer,fs) computes the increase in speech
%   intelligibility of the target when the target and interferer are
%   spatially separated. They are preferably represented by their impulse
%   responses, but can be represented by noise recordings of equivalent
%   spectral shape emitted from the same source locations (using the same
%   noise duration for target and interferer). The impulse responses are
%   assumed to be sampled at a sampling frequency of fs Hz.  If the
%   modelled sources differ in spectral shape, this can be simulated by
%   pre-filtering the impulse responses.
%
%   [benefit, weighted_SNR, weighted_bmld]=JELFS2011(...) additionally
%   returns the benefit from the SII weighted SNR and the SII weighted BMLD.
%
%   If target or interferer are cell-arrays, the HRTF data will be loaded. The first
%   argument in the cell-array is the azimuth angle, and the second
%   parameter is the database type. The elevation is set to zero.
%   function. 
%
%   Example:
%   --------
%
%   The following code will load HRIRs from the 'kemar' database and
%   compute the binaural speech intelligibility advantage for a target
%   at 0 degrees and interferers at 300 and 90 degrees:
%
%     jelfs2011({0,'kemar'},{[330 90],'kemar'})
%
%   See also: culling2004 
% 
%   References:
%     J. Culling, S. Jelfs, and M. Lavandier. Mapping Speech Intelligibility
%     in Noisy Rooms. In Proceedings of the 128th convention of the Audio
%     Engineering Society, Convention paper 8050, 2010.
%     
%     S. Jelfs, J. Culling, and M. Lavandier. Revision and validation of a
%     binaural model for speech intelligibility in noise. Hearing Research,
%     2011.
%     
%     M. Lavandier, S. Jelfs, J. Culling, A. Watkins, A. Raimond, and
%     S. Makin. Binaural prediction of speech intelligibility in reverberant
%     rooms with multiple noise sources. J. Acoust. Soc. Am.,
%     131(1):218--231, 2012.
%     
%
%   Url: http://amtoolbox.org/amt-1.1.0/doc/models/jelfs2011.php

% Copyright (C) 2009-2021 Piotr Majdak, Clara Hollomey, and the AMT team.
% This file is part of Auditory Modeling Toolbox (AMT) version 1.1.0
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

%   #StatusDoc: Good
%   #StatusCode: Perfect
%   #Verification: Unknown
%   #Requirements: MATLAB

%   AUTHOR: Peter Sondergaard
%   AUTHOR: Matthieu Lavandier
%   adapted for AMT by Clara Hollomey (2021)
  
%auditory filterbank either via filterbank or via single gammatone
  definput.flags.verify = {'filterbank', 'single'};
  definput.flags.ears={'both','left','right'};
  definput.keyvals.fs=[];
  definput.keyvals.pad=1024;
  [flags,kv,fs]=ltfatarghelper({'fs'},definput,varargin);
  
  % If target or interferer are cell arrays, load HRTFs.
  if iscell(target)
    X=SOFAload(fullfile(amt_basepath,'hrtf',mfilename,[target{2} '.sofa']));
    idx=find(X.SourcePosition(:,1)==target{1} & X.SourcePosition(:,2)==0);
    target=squeeze(X.Data.IR(idx,:,:))';
    target=postpad(target,size(target,1)+kv.pad);
    fs=X.Data.SamplingRate;
  end
  
  if iscell(interferer)
    azims=numel(interferer{1});
    X=SOFAload(fullfile(amt_basepath,'hrtf',mfilename,[interferer{2} '.sofa']));
    for ii=1:azims
      idx(ii)=find(X.SourcePosition(:,1)==mod(interferer{1}(ii),360) & X.SourcePosition(:,2)==0);
    end
    interferer=shiftdim(X.Data.IR(idx,:,:),2);
    interferer=postpad(interferer,size(interferer,1)+kv.pad);
    fs2=X.Data.SamplingRate;
    
    if fs2~=fs
      error('%s: Mis-match between target and interferer sampling rate.',upper(mfilename));
    end;
    % Old code compatibility
    if ndims(interferer)==3
      s=size(interferer);
      interferer=reshape(interferer,s(1)*s(2),s(3));
      interferer=interferer/sqrt(azims);
    end;
  end;  
  
  if isempty(fs)
    error('%s: You must specify the sampling rate, fs.',upper(mfilename));
  end;
  
if flags.do_filterbank  
  % Make sure that there is at least 1 erb per channel, and get
  % the gammatone filters.
  nchannels=ceil(freqtoerb(fs/2));
  fc=erbspace(26.0508,fs/2,nchannels);
  [b,a] = gammatone(fc,fs,'complex', 'n', 5);
  %gammatone_c(target_in(:,1),fs,fc(n));
  effective_SNR = zeros(nchannels,1);
  bmld_prediction = zeros(nchannels,1);
  
  targ_f = 2*real(ufilterbankz(b,a,target));
  int_f  = 2*real(ufilterbankz(b,a,interferer));
  
  for n = 1:nchannels
    % Calculate the effect of BMLD
    if flags.do_both
      % cross-correlate left and right signal in channel n for both the
      % target and the inteferer
      [phase_t(n), coher_t(n)] = do_xcorr(targ_f(:,n,1),targ_f(:,n,2),fs,fc(n)); 
      [phase_i(n), coher_i(n)] = do_xcorr( int_f(:,n,1), int_f(:,n,2),fs,fc(n)); 
      bmld_prediction(n) = culling2004(coher_i(n),phase_t(n),phase_i(n),fc(n));
    end 
    % Calculate the effect of better-ear SNR
    left_SNR  = sum(targ_f(:,n,1).^2) / sum(int_f(:,n,1).^2);
    right_SNR = sum(targ_f(:,n,2).^2) / sum(int_f(:,n,2).^2);
    
    if flags.do_both
      SNR = max(left_SNR,right_SNR);
    end;
    
    if flags.do_left
      SNR = left_SNR;
    end;
    
    if flags.do_right
      SNR = right_SNR;
    end
    
    % combination
    effective_SNR(n) = 10*log10(SNR);
  end
  effective_SNR = effective_SNR';
  bmld_prediction = bmld_prediction';
else
    nerbs = 1:0.5:round(f2erbrate(fs/2));
    for n = 1:length(nerbs)
        % get filter cf
        fc(n) = round(erbrate2f(nerbs(n)));         
        % filter target and interferer separately
        targ_left = auditoryfilterbank(target(:,1),fs,fc(n), 'lavandier2022');    
        targ_right = auditoryfilterbank(target(:,2),fs,fc(n), 'lavandier2022');   
        int_left = auditoryfilterbank(interferer(:,1),fs,fc(n), 'lavandier2022');       
        int_right = auditoryfilterbank(interferer(:,2),fs,fc(n), 'lavandier2022');
        % BMLD
        [int_phase, int_coherence] = do_xcorr(int_left,int_right,fs,fc(n)); % cross-correlate
        [target_phase] = do_xcorr(targ_left,targ_right,fs,fc(n));    
        bmld_prediction(n) = (bmld(int_coherence,target_phase,int_phase,fc(n)))';
        % better-ear SNR in dB based on energy (independent of 0 padding of
        % BRIRs) energy=10*Log10(sum(sig.*sig))
        left_SNR = 10*log10(sum(targ_left.^2)/sum(int_left.^2));
        right_SNR = 10*log10(sum(targ_right.^2)/sum(int_right.^2));
        effective_SNR(n) = max(left_SNR,right_SNR);
    end
    %integration accross frequency using SII weightings
    %weightings = calc_SII_weightings(fc);
    %weightings = f2siiweightings(fc);
    %weighted_bmld = sum(bmld_prediction'.*weightings',2);
    %weighted_better_ear = sum(better_ear_prediction.*weightings',2);

    %twoears_benefit = weighted_better_ear + weighted_bmld;
    
end


  % Calculate the SII weighting
  weightings = f2siiweightings(fc);
  
  if flags.do_both
    weighted_bmld = sum(bmld_prediction.*weightings');
  else
    weighted_bmld = 0;
  end
  
  weighted_SNR = sum(effective_SNR.*weightings');
  
  benefit = weighted_SNR + weighted_bmld;
end

  
function [phase, coherence] = do_xcorr(left, right, fs, fc)
    [iacc, lags] = xcorr(left,right,round(fs/(fc*2)),'coeff'); %round(fs/(fc*2)) is for conformity with Durlach's 1972 formulation which allows time delays up to 
                                                               %+/- half the period of the channel centre frequency.
    [coherence, delay_samp] = max(iacc);
    phase = fc*2*pi*lags(delay_samp)/fs;
end

