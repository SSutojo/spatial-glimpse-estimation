function [StimSpecsTrain, StimSpecsVal] = evalStimSpecs(SetSize)
%STIMSPECS Summary of this function goes here
%   Detailed explanation goes here

if strcmp(SetSize,'small') == 1
    f = 1;      % yields a total of 512 mixes (=training files)  1
    v = 1/8;    % 64 mixes     1/8
elseif strcmp(SetSize,'medium') == 1
    f = 2;      % 1024 mixes
    v = 1/4;    % 128 mixes
elseif strcmp(SetSize,'large') == 1
    f = 3;      % 4096 mixes   f=8 
    v = 1/2;    % 512 mixes       v = 1
elseif strcmp(SetSize,'test') == 1
    f = 1/8;
    v = 1/8;    
else
    error('define a valid SetSize')
end


% make 50 mixes that are all FMFM, then remove azimuth data later

%% list for training/validation Set
StimSpecsTrain{1,1} = 10;      StimSpecsTrain{1,2} = 'FMFMF';       StimSpecsTrain{1,3} = 10;       StimSpecsTrain{1,4} = 3;
StimSpecsTrain{2,1} = 10;      StimSpecsTrain{2,2} = 'MFMFM';       StimSpecsTrain{2,3} = 10;       StimSpecsTrain{2,4} = 3;
StimSpecsTrain{3,1} = 10;      StimSpecsTrain{3,2} = 'FMFMF';       StimSpecsTrain{3,3} = 10;       StimSpecsTrain{3,4} = 3;
StimSpecsTrain{4,1} = 10;      StimSpecsTrain{4,2} = 'MFMFM';       StimSpecsTrain{4,3} = 10;       StimSpecsTrain{4,4} = 3;
StimSpecsTrain{5,1} = 30;      StimSpecsTrain{5,2} = 'FMFMF';       StimSpecsTrain{5,3} = 10;       StimSpecsTrain{5,4} = 3;
StimSpecsTrain{6,1} = 30;      StimSpecsTrain{6,2} = 'MFMFM';       StimSpecsTrain{6,3} = 10;       StimSpecsTrain{6,4} = 3;

StimSpecsVal = [];
% create a really small set for debugging
%StimSpecsTrain{1,1} = 20;      StimSpecsTrain{1,2} = 'FMFM';       StimSpecsTrain{1,3} = 10;       StimSpecsTrain{1,4} = 3;



end

