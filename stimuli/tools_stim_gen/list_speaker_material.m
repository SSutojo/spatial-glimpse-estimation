function speaker_material = list_speaker_material(sub_dir, speaker_ID, sentence_type)
% to be used with the timit database
% returns all files that belong to the given speaker
% timit wavread function is needed to open these files
% 
% INPUTS :  sub_dir - subdirectory in which the speech files are stored
%           speaker_ID - string with the speaker identity (see timit readme)
%           sentence_type - also see the timit readme. 4 inputs are
%           possible here: SA (dialect sentence), SX (ponetically compact
%           sentences), SI (phonetically diverse sentences), 'all' will
%           list all three sentence types e.g. all sentences in the folder
% OUTPUT :  speaker_material - struct with all time signals for this ID 
%
%
% 1. finds directory whtat holds the speaker material
% 2. lists all wavefiles

curr_dir = pwd; % path of the current directory
speaker_dir = fullfile(curr_dir, sub_dir, speaker_ID); % path for the speaker's speech material
cd(speaker_dir); % go to the directory with speaker material

if strcmp('SA',sentence_type) || strcmp('SX',sentence_type) || strcmp('SI',sentence_type) == 1
    
    wav_list = ls([sentence_type,'*.wav']); % list all .wav files of desired sentence type that are located in this directory
    
    for ii = 1:size(wav_list,1)
        filename = num2str(wav_list(ii,:));
        filename = strtrim(filename); % remove white space from filename
        addpath(fullfile(curr_dir,'TIMITutils'));
        [x] = wavReadTimit(filename);
        speaker_material(1,ii).time_signal = x; % number of fields is equal to the number of sentences
    end
    cd(curr_dir); % return to previous working directory
    
    
elseif strcmp('all',sentence_type)
    list1 = ls(['SA','*.wav']); % list all three sentence types
    list2 = ls(['SX','*.wav']);
    list3 = ls(['SI','*.wav']);
    wav_list = [cellstr(list1) ; cellstr(list2); cellstr(list3)];
    wav_list = char(wav_list);
    
    for jj = 1:size(wav_list,1)
        filename = num2str(wav_list(jj,:));
        filename = strtrim(filename); 
        addpath(fullfile(curr_dir,'TIMITutils'));
        [x] = wavReadTimit(filename);
        speaker_material(1,jj).time_signal = x;
    end
    cd(curr_dir); 
else    
    error('invalid sentence type')
end 

end




