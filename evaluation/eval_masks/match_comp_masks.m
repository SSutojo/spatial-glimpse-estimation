function [ matched_est_map, matched_est_IBMs, tru_IBMs, eval_results] = match_comp_masks( est_source_map, tru_source_map)
%MATCH_COMP_MASKS matches and compares two source maps
%   est_source_map and ideal_source_map are both of same size (number of
%   time frames X number of filter bands). Each t-f unit is labeled with an
%   integer that represents the dominant source.
%   In a first step the source maps are matched (e.g. in the estimated
%   source map, the label 1 might represent the same source as the label 3
%   in the ideal source map). This is done by assigning each estimated
%   source the ideal source with wich it has the largest overlap.
%   INPUTS:
%           est_source_map  - estimated source map
%           tru_source_map  - ideal/true source map (obtained from oracle
%           information)
%   OUTPUTS:
%           matched_est_map - the estimated source map after matching
%           matched_est_IBMs- in each layer is the IBM per estimated source
%           tru_IBMs        - in each layer is the IBM per true source
%           per_corr_per_source - vector with correctly identified
%                             t-f-units per true source
%           eval_results    - struct with evaluation results (percent correctly classified t-f units, balanced accuracy,
%                           vector with correctly identified t-f-units per true source)
% ****
% EXTRA CASES for cases in which the number of estimated
% sources is not equal to the true number of sources


frames                      = size(est_source_map,1);
bands                       = size(est_source_map,2);

est_ID_vec                  = unique(est_source_map);                   % different IDs available in est_source_map
est_num_sources             = numel(est_ID_vec);                        % estimated number of sources

tru_ID_vec                  = unique(tru_source_map);                   % different IDs available in tru_source_map
tru_num_sources             = numel(tru_ID_vec);                      % true number of sources

%% match estimated source maps and ground truth source map
% the source maps should have the maximum possible overlap assign the estimated IBM to the ground truth IBM that covers most of it
% create a binary version of each source map (each source is a layer in which the associated t-f units are set to ones)
% for every estimated source, calculate the overlap with a tru source
est_IBMs                    = zeros(size(est_source_map,1),size(est_source_map,2),est_num_sources);
tru_IBMs                    = zeros(size(tru_source_map,1),size(tru_source_map,2),tru_num_sources);
for aa = 1:est_num_sources
    aa_ID                   = est_ID_vec(aa);
    est_IBMs(:,:,aa)        = est_source_map == aa_ID;                      % set all pixels of class aa to 1 and everything else to 0
end
for bb = 1:tru_num_sources
    bb_ID                   = tru_ID_vec(bb);
    tru_IBMs(:,:,bb)        = tru_source_map == bb_ID;
end

% calculate the overlap between each class in the est_source_map and the tru_source map 
overlap_mat                 = zeros(est_num_sources, tru_num_sources, 3);  % overlap_mat has 3 layers
                                                                    % 1st layer: overlap of class n in est_map and class n in tru_map,
                                                                    % 2nd layer: size of class n in est_source_map,
                                                                    % 3rd layer: size of class n in tru_source_map
for cc = 1:est_num_sources                                          % each source is treated as a class
    combined_mat            = repmat(est_IBMs(:,:,cc),[1 1 tru_num_sources]).*tru_IBMs;
    overlap_mat(cc,:,1)     = sum(sum(combined_mat,1),2);                                % overlap of the class cc in est_map and tru_map
    overlap_mat(cc,:,2)     = repmat(sum(sum(est_IBMs(:,:,cc))), [1 tru_num_sources]);   % size of class cc in est_source_map
    overlap_mat(cc,:,3)     = sum(sum(tru_IBMs(:,:,:),1),2);                             % size of classes 1 to end in tru_source_map
end
intersection_mat            = overlap_mat(:,:,1);       % # of overlapping datapoints between each class in est_source_map and tru_source_map

eval_results.est_num_sources = est_num_sources;
eval_results.tru_num_sources = tru_num_sources;

%%

if est_num_sources == tru_num_sources
    
    matched_IDs                 = zeros(est_num_sources,2); % column vector with the ID of the est. source in 1st col. and the matched tru source ID in 2nd col.
    for dd = 1:est_num_sources
        [~ , I1]                = max(max(intersection_mat,[],2));
        [~ , I2]                = max(intersection_mat(I1,:));
        
        matched_IDs(dd,:)       = [est_ID_vec(I1),tru_ID_vec(I2)]; % 1st col. is the index of the est. source,  2nd column the true source with largest overlap
        intersection_mat(I1,:)  = ones.*(-1);           % set the used rows and columns to minus 1 so they won't be used again
        intersection_mat(:,I2)  = ones.*(-1);
    end
    
    % switch the IDs in the est_source_map to the matched IDs
    gl_mat                  = zeros(size(est_source_map,1),size(est_source_map,2),est_num_sources);
    for ee = 1:est_num_sources
        gl_mat(:,:,ee)      = est_source_map == matched_IDs(ee,1);           % set the t-f units to 1 that belong to the ID in the estimated mask
        gl_mat(:,:,ee)      = gl_mat(:,:,ee).*matched_IDs(ee,2);             % multiply the ones with the matched ID
    end
    matched_est_map = sum(gl_mat,3);
    %% calculate the percent correct classified (overlap of ideal and estimated mask after matching)        
    % count the number of correclty classified time-frequency units
    correct_units               = sum(sum(matched_est_map==tru_source_map));     
    % calculate as percentage of all t-f-units
    all_units                   = bands*frames;                      
    per_corr                    = correct_units/(all_units/100);    
    %% generate further outputs  
    matched_est_IBMs                    = zeros(size(est_IBMs));        % the matched_est_IBMs instead of est_IBMs
    for ff = 1: size(matched_IDs,1)
        true_ID                         = matched_IDs(ff,2);
        true_ID_ind                     = find(tru_ID_vec==true_ID);  % index of true ID within the origninal ID vector
        matching_est_ID                 = matched_IDs(ff,1);
        matching_est_ID_ind             = find(est_ID_vec==matching_est_ID);      
        matched_est_IBMs(:,:,true_ID_ind)= est_IBMs(:,:,matching_est_ID_ind);      
    end       
    per_corr_per_source                 = zeros(tru_num_sources,1);     % percent correct per per true source
    bacc_per_source                     = zeros(tru_num_sources,1);     % balanced accuracy for each source (tp/p + tn/n)/2
    for gg = 1:tru_num_sources 
        positives                       = sum(sum(tru_IBMs(:,:,gg)));   % total number of positives for this source (= 1s in the IBM)
        negatives                       = all_units-positives;     % total number of negatives for this source (= 0s in the IBM)
        true_pos                        = sum(sum(tru_IBMs(:,:,gg)==matched_est_IBMs(:,:,gg)));  % number of overlapping ones in true and est IBM
        true_neg                        = length(find((tru_IBMs(:,:,gg)+matched_est_IBMs(:,:,gg))==0));% alle Nullen, die auch als Nullen erkannt wurden       
        bacc_per_source(gg)             = (true_pos/positives + true_neg/negatives)/2;
        per_corr_per_source(gg)         = true_pos/all_units*100;        
    end 
    
    eval_results.bacc_per_source             = bacc_per_source;   
    eval_results.per_corr_per_source        = per_corr_per_source;
    eval_results.per_corr = per_corr;
    
elseif est_num_sources < tru_num_sources
%     disp('number of estimated sources is smaller than number of true sources')
    % assign each est_source to the tru_sources with the largest overlap
    matched_IDs                 = zeros(est_num_sources,2); % column vector with the ID of the est. source in 1st col. and the matched tru source ID in 2nd col.
    for dd = 1:est_num_sources
        [~ , I1]                = max(max(intersection_mat,[],2));
        [~ , I2]                = max(intersection_mat(I1,:));
        
        matched_IDs(dd,:)       = [est_ID_vec(I1),tru_ID_vec(I2)]; % 1st col. is the index of the est. source,  2nd column the true source with largest overlap
        intersection_mat(I1,:)  = ones.*(-1);           % set the used rows and columns to minus 1 so they won't be used again
        intersection_mat(:,I2)  = ones.*(-1);
    end
    
    % switch the IDs in the est_source_map to the matched IDs
    gl_mat = zeros(size(est_source_map,1),size(est_source_map,2),est_num_sources);
    for ee = 1:est_num_sources
        gl_mat(:,:,ee) = est_source_map == matched_IDs(ee,1);           % set the t-f units to 1 that belong to the ID in the estimated mask
        gl_mat(:,:,ee) = gl_mat(:,:,ee).*matched_IDs(ee,2);             % multiply the ones with the matched ID
        
    end
    matched_est_map = sum(gl_mat,3);
    
    %% calculate the percent correct classified (overlap of ideal and estimated mask after matching)
    % count the number of correclty classified time-frequency units
    correct_units               = sum(sum(matched_est_map==tru_source_map));
    
    % calculate as percentage of all t-f-units
    all_units                   = bands*frames;
    per_corr                    = correct_units/(all_units/100);
    
    %% generate further outputs
    matched_est_IBMs                    = zeros(size(est_IBMs));        % the matched_est_IBMs instead of est_IBMs
    for ff = 1: size(matched_IDs,1)
        true_ID                         = matched_IDs(ff,2);
        true_ID_ind                     = find(tru_ID_vec==true_ID);  % index of true ID within the origninal ID vector
        matching_est_ID                 = matched_IDs(ff,1);
        matching_est_ID_ind             = find(est_ID_vec==matching_est_ID);
        matched_est_IBMs(:,:,true_ID_ind)= est_IBMs(:,:,matching_est_ID_ind);
    end
    
    per_corr_per_source                 = zeros(tru_num_sources,1);     % percent correct per per true source
    bacc_per_source                     = zeros(tru_num_sources,1);     % balanced accuracy for each source (tp/p + tn/n)/2
    for gg = 1:est_num_sources
        positives                       = sum(sum(tru_IBMs(:,:,gg)));   % total number of positives for this source (= 1s in the IBM)
        negatives                       = all_units-positives;     % total number of negatives for this source (= 0s in the IBM)
        true_pos                        = sum(sum(tru_IBMs(:,:,gg)==matched_est_IBMs(:,:,gg)));  % number of overlapping ones in true and est IBM
        true_neg                        = length(find((tru_IBMs(:,:,gg)+matched_est_IBMs(:,:,gg))==0));% alle Nullen, die auch als Nullen erkannt wurden
        
        bacc_per_source(gg)             = (true_pos/positives + true_neg/negatives)/2;
        per_corr_per_source(gg)         = true_pos/all_units*100;
    end
    
    eval_results.bacc_per_source             = bacc_per_source;
    eval_results.per_corr_per_source         = per_corr_per_source;
    eval_results.per_corr                     = per_corr;
    
elseif est_num_sources > tru_num_sources
%     disp('number of estimated sources is larger than number of true sources')
    matched_IDs                 = zeros(tru_num_sources,2); % column vector with the ID of the est. source in 1st col. and the matched tru source ID in 2nd col.
    % theres not a matched true source for every estimated source
    for dd = 1:tru_num_sources
        [~ , I1]                = max(max(intersection_mat,[],2));
        [~ , I2]                = max(intersection_mat(I1,:));
        
        matched_IDs(dd,:)       = [est_ID_vec(I1),tru_ID_vec(I2)]; % 1st col. is the index of the est. source,  2nd column the true source with largest overlap
        intersection_mat(I1,:)  = ones.*(-1);           % set the used rows and columns to minus 1 so they won't be used again
        intersection_mat(:,I2)  = ones.*(-1);
    end

    gl_mat = zeros(size(est_source_map,1),size(est_source_map,2),tru_num_sources);
    for ee = 1:tru_num_sources
        gl_mat(:,:,ee) = est_source_map == matched_IDs(ee,1);           % set the t-f units to 1 that belong to the ID in the estimated mask
        gl_mat(:,:,ee) = gl_mat(:,:,ee).*matched_IDs(ee,2);             % multiply the ones with the matched ID   
    end    
    matched_est_map = sum(gl_mat,3);
    
    %% calculate the percent correct classified (overlap of ideal and estimated mask after matching)
    % count the number of correclty classified time-frequency units
    correct_units               = sum(sum(matched_est_map==tru_source_map));
    
    % calculate as percentage of all t-f-units
    all_units                   = bands*frames;
    per_corr                    = correct_units/(all_units/100);
    
    %% generate further outputs
    matched_est_IBMs                    = zeros(size(est_IBMs));        % the matched_est_IBMs instead of est_IBMs
    for ff = 1:tru_num_sources
        true_ID                         = matched_IDs(ff,2);
        true_ID_ind                     = find(tru_ID_vec==true_ID);  % index of true ID within the origninal ID vector
        matching_est_ID                 = matched_IDs(ff,1);
        matching_est_ID_ind             = find(est_ID_vec==matching_est_ID);
        matched_est_IBMs(:,:,true_ID_ind)= est_IBMs(:,:,matching_est_ID_ind);
    end
    
    per_corr_per_source                 = zeros(tru_num_sources,1);     % percent correct per per true source
    bacc_per_source                     = zeros(tru_num_sources,1);     % balanced accuracy for each source (tp/p + tn/n)/2
    for gg = 1:tru_num_sources
        positives                       = sum(sum(tru_IBMs(:,:,gg)));   % total number of positives for this source (= 1s in the IBM)
        negatives                       = all_units-positives;     % total number of negatives for this source (= 0s in the IBM)
        true_pos                        = sum(sum(tru_IBMs(:,:,gg)==matched_est_IBMs(:,:,gg)));  % number of overlapping ones in true and est IBM
        true_neg                        = length(find((tru_IBMs(:,:,gg)+matched_est_IBMs(:,:,gg))==0));% alle Nullen, die auch als Nullen erkannt wurden
        
        bacc_per_source(gg)             = (true_pos/positives + true_neg/negatives)/2;
        per_corr_per_source(gg)         = true_pos/all_units*100;
    end    
    
    eval_results.bacc_per_source             = bacc_per_source;
    eval_results.per_corr_per_source        = per_corr_per_source;
    eval_results.per_corr                     = per_corr;
    
    
end

end

