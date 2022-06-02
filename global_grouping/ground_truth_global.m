function [ est_source_map, glimpseLabels, glimpsePairLabels, glimpsePairSim] = ground_truth_global( glimpse_map, ground_truth, params )
%GLOBAL_GROUPING assign the local clusters to same sources
%   get an estimate of the source map (est_source_map labels each t-f-tile
%   with an integer associated with either a source or background)
%   
% INPUTS:   - clustering_out: This is the output strucutre of the local clustering stage
%           - feature_mats: structure that holds different features types
%             feature matrices have the dimensions (length_of_feature_vector X num_win X num_filterbands)
%           - mode{1}: the type of feature that is used for global grouping
%             ('azimuth_probs','log_azimuth_probs','periodicity_degree')
%             here, the mode could possibly be set to a feature combination
%             but this mode has to be installed
% OUTPUTS:  - est_source_map: a map of size ....  in which every
%             time-frequency unit is labeled with an integer that indicates 
%             the dominant source
%           - peak_features : the most pronounced features e.g. the most
%           frequently estimated azimuths (most likely speaker positions or
%           the most pronounced/ frequently occurring fundamental
%           frequencies


%% get the local clusters
ideal_source_map                = ground_truth.ideal_source_map;
num_glimpses                    = max(glimpse_map(:)); % get the number of separate clusters

%% for each glimpse, average azimuth probabilities
% glimpse_azs                 = zeros(num_glimpses,4);   
% num_sources                 = max(ideal_source_map(:));

glimpse_sourceIDs           = zeros(num_glimpses,4);
glimpse_sourceSNRs          = zeros(num_glimpses,ground_truth.num_components);
for zz = 1:num_glimpses
    [rowz,colz]             = find(glimpse_map == zz);              % get row/column indices of the elements in the regarded glimpse
    size_zz                 = length(rowz);                         % number of t-f-units units in the regarded glimpse 
    source_IDs              = zeros(size_zz,1);
    for uu = 1:size_zz
        source_IDs(uu) = ideal_source_map(rowz(uu),colz(uu));   % affilitation of each t-f unit in the glimpse with a source
    end
    
    
    for vv = 1:ground_truth.num_components
        noise_sum = 0;
        sig_sum   = 0;
        for uu = 1:size_zz
            noise_mat = ground_truth.power_mat_sep(rowz(uu),colz(uu),setdiff(1:end,vv));
            noise_sum = noise_sum + sum(noise_mat,3);
            sig_sum   = sig_sum + ground_truth.power_mat_sep(rowz(uu),colz(uu),vv);
        end
        glimpse_sourceSNRs(zz,vv) = 10*log10(sig_sum/noise_sum);
    end
    
    [n,bin]  = hist(source_IDs,unique(source_IDs));      
    highest_bins = find(n==max(n)); % 
    
    if length(highest_bins)==1
        source = bin(highest_bins);  % most frequent source in the glimpse
        occ_part =  max(n)/size_zz;
    elseif length(highest_bins)>1   % glimpse cannot be assigned, therefore pick the first of the highest bins (Alternative is unten auskommentier)
        source = bin(highest_bins(1));    
        occ_part =  max(n)/size_zz;    
    else
        disp('warning: empty glimpse')
    end
      
    glimpse_sourceIDs(zz,1)       = zz;                      % 1st col.= number of the glimpse 
    glimpse_sourceIDs(zz,2)       = source;               % 2nd col.= ID of the dominant source (=source that dominates most t-f units in glimpse)
    glimpse_sourceIDs(zz,3)       = occ_part;                      % 3rd col.= percentage of  t-f units in glimpse occupied by dominant source
    glimpse_sourceIDs(zz,4)       = size_zz;                  % 4th col.= number of t-f-units in the glimpse    
end

%%
glimpse_map_mat             = zeros(size(glimpse_map,1),size(glimpse_map,2),num_glimpses);

for pp = 1:num_glimpses                                             % for each glimpse: check if it is near an overall_feature_peak    
    curr_sourceID           = glimpse_sourceIDs(pp,2);
    glimpse_map_mat(:,:,pp) = (glimpse_map==pp).*curr_sourceID;

end

est_source_map              = sum(glimpse_map_mat,3); 
glimpseLabels               = glimpse_sourceIDs(:,2);


%% ideal glimpse Pair Labels
glimpsePairLabels               = zeros(((num_glimpses*(num_glimpses-1))/2),1);
counter                 = 1;
for cc = 1:(num_glimpses-1)
    for dd = (cc+1):num_glimpses
        if glimpseLabels(cc) == glimpseLabels(dd)
            glimpsePairLabels(counter,1) = 1;
        elseif glimpseLabels(cc) ~= glimpseLabels(dd)
            glimpsePairLabels(counter,1) = 0;
        end
        counter = counter + 1;
    end
end


%% ideal glimpse Pair Similarity
glimpsePairSim    = zeros(((num_glimpses*(num_glimpses-1))/2),1);
sigm_SNRs               = 1./(1+exp(-glimpse_sourceSNRs)); % use sigmoid func to derive sth like a source presence probability from the local SNRs
source_probs            = sigm_SNRs./sum(sigm_SNRs,2); % normalize to get a probability
counter = 1;
for cc = 1:(num_glimpses-1)
    for dd = (cc+1):num_glimpses
        
        glimpsePairSim(counter,1) = sum(source_probs(cc,:).*source_probs(dd,:));
        counter = counter + 1;
    end
end


end


