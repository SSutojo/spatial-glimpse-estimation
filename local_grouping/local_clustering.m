function [clustering_out] = local_clustering(feature_mats, params)
%LOCALLY GROUP THE TIME-FREQUENCY UNITS ACCORDING TO THE FEATURES IN EACH UNIT 
%PRIOR TO THE GROUPING, A SIMILARITY VALUE (SV) IS CALCULATED FOR EACH
%TRANSITION BETWEEN TWO ADJACENT TIME-FREQUENCY UNITS
% 
% INPUTS:   features_all.feature_mats: feature matrices that are used for
%                       estimation of local contrasts and glimpses
%                       dimensions (length(feature_vec) X number_of_time_frames X number_of_filterbands)
%           params:     parameter settings (see "set_local_params")
%
% OUTPUTS:  clustering_out. ...
%           .SV_mat:    matrix with an SV for every transition between two
%                       adjacent t-f-units. size: num_win*2-1,nFilts*2-1,num_SV_types
%           .contrast_map: matrix with a contrast value (range:0-1) for every transition
%                       and interpolated values at vertices and T-F units,
%                       this matrix can be used for visualization
%                       (1=high contrast, 0=no contrast)
%                       size: num_win*2-1,nFilts*2-1,num_SV_types
%           .contours:  matrix with 1s at each contour and 0s at each joint
%                       transition (same size as SV_mat). Basically, the
%                       SV_mat after applying a criterion
%           .regions_est: matrix, same size as .contours but the regions enclosed by the contours
%                       are filled with ongoing integers
%           .glimpse_map: matrix displaying the local clusters, each
%                       t-f-unit is labeled with an integer that labels the
%                       local cluster it belongs to (like a glimpse)
%           .bin_mask:  if a binary mask was applied it is also stored here
%                       so that it can be checked or plottet in the
%                       analysis
%           .SV_vert:   all across-frequency-band SV values
%           .SV_horz:   all within-frequency-band SV values
%           .SV_types:  types of relative features that are used for contrast estimation


% get number of time frames and filter bands based on the first field (feature type)
all_fields  = fieldnames(feature_mats);
first_field = all_fields{1,1};
long_name   = ['feature_mats.',first_field];
if strcmp(first_field,'power_mat')==1  % in this case the dimensions of the feature_mat are (num_win X nFilts)
    num_win = size(eval(long_name),1);
    nFilts = size(eval(long_name),2);
else                                        % in all other cases the dims are (length(feature_vec) X num_win X nFilts)
    num_win = size(eval(long_name),2);
    nFilts = size(eval(long_name),3);
end

%% calculate the relative features
bin_mask                                    = ones(num_win,nFilts); % currently, no binary mask is used, so this is set to 1 for every T-F unit, but option can be useful when VAD or sth should be integrated here
mask_opts{1}                                = 'lc_WITHOUT_bin_mask';
mask_opts{2}                                = bin_mask;   

[ SV_mat_raw, SV_horz_raw, SV_vert_raw ]    = combined_SV_vec( feature_mats, params.local_params, mask_opts, params.local_params.SV_types); % params.local_params.SV_types defines which SVs (features and correlation types) are used for local clustering

%% combine different SV types with a pre-trained model (get model name and load model if it already exists)
if strcmp(params.local_params.model_confs.model_type,'none') == 1
    SV_mat                      = SV_mat_raw;
    SV_horz                     = SV_horz_raw;
    SV_vert                     = SV_vert_raw;
    clustering_out.model_name   = SVtypes2name(params.local_params.SV_types);
else
    %************** get name of desired model and load it from file *****************************************************  
    hash_struct.feature_params  = params.feature_params;
    hash_struct.SV_types        = params.local_params.SV_types;
    hash_struct.model_confs     = params.local_params.model_confs;
    [model_name]                = DataHash(hash_struct);    %    model_name =[model_name,'balanced'];

    C                   = load([pwd,filesep,'local_grouping',filesep,'contrast_estimation',filesep,'trained_models',filesep,model_name,'.mat']);
    modelFinal          = C.modelFinal;
    
    
    %% these should be the same for all cases
    % get number of filterbands, time frames and SV (relative feature) types
    nFilts      = size(SV_horz_raw,2);  nFrames = size(SV_vert_raw,1);  nSVtypes = size(SV_vert_raw,3);        
    
    fileIdx_ac  = feature_mats.fileIdx_ac;
    fileIdx_inb = feature_mats.fileIdx_inb;
    num_files   = size(fileIdx_ac,1);
    
    SV_mat          = ones(size(SV_mat_raw,1),size(SV_mat_raw,2));      % matrix for all similarity values (SVs), will later be used to get contrast map
    SV_vert         = zeros(nFrames   , nFilts-1);                      % across band similarity values
    SV_horz         = zeros(nFrames-1 , nFilts  );                      % within band similarity values

    %% different options for model range (broadband either joint model for across and within subband, or separately, or subband specific models)
    % _______________________________________________________________________________________________________________________________        
    % model includes all filterbands at the same time (broadband) and combines across and within filterband transitions    
    if strcmp(params.local_params.model_confs.model_range,'broadband_all') == 1         
        if strcmp(params.local_params.model_confs.model_type,'dnn') == 1
            %************** prepare the features for classification ****************************************************** 
            ac_SVs                  = zeros((nFilts-1)*nSVtypes,nFrames); % across freq band
            inb_SVs                 = zeros(nFilts*nSVtypes,nFrames-1);   % within frequency band similarity values
            for aa = 1:nFilts-1
                ac_SVs(((aa-1)*nSVtypes+1):(aa*nSVtypes),:)   = squeeze(SV_vert_raw(:,aa,:)).';
            end
            for aa = 1:nFilts
                inb_SVs(((aa-1)*nSVtypes+1):(aa*nSVtypes),:)  = squeeze(SV_horz_raw(:,aa,:)).'; % the entry from frame x, band y and feature type z can be found in the dimension: frame x, new_dim (y-1)*num_feats + z
            end
            all_SVs                 = [];       % stack across and within band
            for aa = 1:num_files
                startIdx_ac         = fileIdx_ac(aa,1);
                endIdx_ac           = fileIdx_ac(aa,2);
                startIdx_inb        = fileIdx_inb(aa,1);
                endIdx_inb          = fileIdx_inb(aa,2);
                temp_features       = cat(1,ac_SVs(:,startIdx_ac:endIdx_ac-1),inb_SVs(:,startIdx_inb:endIdx_inb));  % discards the across band SV value from the last frame of each file
                all_SVs             = cat(2,all_SVs,temp_features); % first half of the dims are the across band SV vals/labels, the second half are the within band SVs/labels
            end
            FEATURES_TEST           = postProcessFeatures(all_SVs.',modelFinal);    % context etc. 
            %************** classify and write results into the SV_mat ***************************************************** 
            predicted               = nntest(modelFinal.DNN, FEATURES_TEST);
            SV_all                  = 1-predicted;
            for bb = 1:(nFilts-1)
                for aa = 1:(nFrames-1)            % for the broadband_all mode the last frame is not considered because . The first part of the vector is the across band SV values
                    SV_mat((aa*2-1),bb*2) = SV_all(aa,bb);
                    SV_vert(aa,bb)        = SV_all(aa,bb);
                end
            end
            for bb = nFilts:(2*nFilts-1)
                for aa = 1:(nFrames-1)
                    SV_mat(aa*2,(bb-nFilts+1)*2-1)   = SV_all(aa,bb);
                    SV_horz(aa,(bb-nFilts+1))        = SV_all(aa,bb);
                end
            end
        else
            error('define valid model type')
        end
    % _______________________________________________________________________________________________________________________________        
    % model includes all filterbands at the same time (broadband) but treats across and within filterband transitions separately (sep)
    elseif strcmp(params.local_params.model_confs.model_range,'broadband_sep') == 1       
        if strcmp(params.local_params.model_confs.model_type,'dnn') == 1
            inBandDNN       = modelFinal.inBand;
            acBandDNN       = modelFinal.acBand;
            %************** prepare the features for classification ****************************************************** 
            ac_SVs          = zeros((nFilts-1)*nSVtypes,nFrames); % across frequency band similarity values
            inb_SVs         = zeros(nFilts*nSVtypes,nFrames-1);   % within frequency band similarity values
            for aa = 1:nFilts-1
                ac_SVs(((aa-1)*nSVtypes+1):(aa*nSVtypes),:)   = squeeze(SV_vert_raw(:,aa,:)).';
            end
            for aa = 1:nFilts
                inb_SVs(((aa-1)*nSVtypes+1):(aa*nSVtypes),:)  = squeeze(SV_horz_raw(:,aa,:)).'; % the entry from frame x, band y and feature type z can be found in the dimension: frame x, new_dim (y-1)*num_feats + z
            end
            [FEATURES_TESTac] = postProcessFeatures(ac_SVs.',acBandDNN); % input must be (nFrames X nFeats(e.g.nBands))
            [FEATURES_TESTin] = postProcessFeatures(inb_SVs.',inBandDNN);
            %************** classify and write results into the SV_mat ***************************************************** 
            SV_vert1    = nntest(acBandDNN.DNN, FEATURES_TESTac);
            SV_horz1    = nntest(inBandDNN.DNN, FEATURES_TESTin);
            SV_vert     = 1-SV_vert1;       % to make this an SV and not contrast val
            SV_horz     = 1-SV_horz1;
            for cc = 1:(nFilts-1)
                for dd = 1:num_win
                    SV_mat(dd*2-1,cc*2) = SV_vert(dd,cc);
                end
            end
            for ee = 1:nFilts
                for ff = 1:num_win-1
                    SV_mat(ff*2,ee*2-1) = SV_horz(ff,ee);
                end
            end
        else
            error('define valid model type')
        end
    % _______________________________________________________________________________________________________________________________        
    % model for each filterband (subband) separately        
    elseif strcmp(params.local_params.model_confs.model_range,'subband_sep') == 1   
        if strcmp(params.local_params.model_confs.model_type,'dnn') == 1
            inBandDNN = modelFinal.inBand;
            acBandDNN = modelFinal.acBand;
            %************** prepare the features for classification and classify ************************************** 
            for bb= 1:nFilts                                        % within band transitions
                SV_vec                  = SV_mat_raw(2:2:end,bb*2-1,:);
                SV_vec                  = (squeeze(SV_vec)).';
                FEATURES1               = postProcessFeatures(SV_vec.',inBandDNN);
                inBand_contour_probs    = nntest(inBandDNN.DNN{1,bb}, FEATURES1);
                for win = 1:size(inBand_contour_probs,1)
                    SV_horz(win,bb)         = 1-inBand_contour_probs(win,1);
                    SV_mat(win*2,bb*2-1)    = 1-inBand_contour_probs(win,1);
                end
            end
            for bb= 1:nFilts-1                                  % across band transitions
                SV_vec                  = SV_mat_raw(1:2:end,bb*2,:);
                SV_vec                  = (squeeze(SV_vec)).';
                FEATURES1               = postProcessFeatures(SV_vec.',acBandDNN);
                acBand_contour_probs    = nntest(inBandDNN.DNN{1,bb}, FEATURES1);
                for win = 1:size(acBand_contour_probs,1)
                    SV_vert(win,bb)         = 1-acBand_contour_probs(win,1);
                    SV_mat(win*2-1,bb*2)    = 1-acBand_contour_probs(win,1);
                end
            end
        else
            error('define valid model type')
        end
    else
        error('define valid model range in local_params.model_confs')
    end
    clustering_out.model_name   = model_name;
end

%% **************** Contrast Map *********************************************
% convert SV_mat into contrast map, use bilinear interpolation for the missing values at edges/vertices, 
% also fill in interpolated values at placeholders for T-F units(this is just for better visualization and doesn't affect glimpse formation) 
inv_SV_mat      = 1-SV_mat;     % contrast value is 1-similarity value
contrast_map    = bilinear_interpolation(inv_SV_mat, 'favor_high_contrasts','vertices_and_TFunits'); % weighting option = 'none', 'favor_low_contrasts', 'favor_high_contrasts'

%% **************** Glimpse Map ************************************************
% derive contours of glimpses from the contrast map (different methods are available)
if strcmp(params.local_params.contour_method,'regiongrow')==1 % use a threshold and then regiongrow (params.threshold_factor affects the size of the local clusters)
    %******** 1. set a threshold to contrast map
    thresh_regiongrow   = params.local_params.regiongrow.contrast_thresh;   
    contrast_map        = squeeze(contrast_map(:,:,1)); 
    contours            = zeros(size(contrast_map));
    for aa = 2:2:size(contours,1)
        for bb = 1:size(contours,2)
            contours(aa,bb) = (contrast_map(aa,bb) >= thresh_regiongrow);
        end
    end
    for aa = 1:size(contours,1)
        for bb = 2:2:size(contours,2)
            contours(aa,bb) = (contrast_map(aa,bb) >= thresh_regiongrow);
        end
    end
    %******** 2. fill the resulting contours with regiongrow
    [regions_est, ~, ~, ~]  = regiongrow(contours, 0, 0.2); % fill the contours with numbers.
    %******** 3. get glimpse map by removing the contours from regions_est and just keeping the glimpse assignment for each t-f-unit
    glimpse_map             = zeros(num_win,nFilts); 
    for bb = 1:num_win % time frames
        for cc = 1:nFilts % filterbands
            glimpse_map(bb,cc) = regions_est((2*bb)-1,(2*cc)-1);
        end
    end  
    
elseif strcmp(params.local_params.contour_method,'boundary_expansion')==1 % boundaries are expanded, before using regiongrow, to close small gaps in the boundaries    
    %******** 1. use a contrast map in which only values at the vertices are interpolated
    contrast_map_boundary_ex  = bilinear_interpolation(inv_SV_mat, 'favor_high_contrasts','vertices_only'); 
    
    %******** 2. set a threshold to contrast map and expand resulting boundaries by Nexpansions pixels in each direction
    thresh_vals         = params.local_params.rep_bound.contrast_threshs;  % threshold to initially convert contrast map into binary map, this can also be a vector of different thresholds!
    Nexpansions         = params.local_params.rep_bound.Nexpansions;  % number of expansions, must be an integer    
    [regions_est,~,~,~] = SoftBoundaryExp(contrast_map_boundary_ex,thresh_vals,Nexpansions); % can be done for different criteria
    
    %******** 3. get resulting glimpse map and assign individual labels to the leftover (single) pixels
    glimpse_map         = get_glimpse_map(regions_est);
    [ro,co]             = find(glimpse_map==0);  % find all zeros
    newIDs              = [max(glimpse_map(:)):1:(max(glimpse_map(:))+length(ro)-1)]; % give each t-f-unit that is labeled with zero, an individual integer as label    
    for idx = 1:length(newIDs) 
        glimpse_map(ro(idx),co(idx))=newIDs(idx);
    end
    contours            = get_contours(glimpse_map);
    
elseif strcmp(params.local_params.contour_method,'DilationErosion') == 1 % repeated boundary/regions expansion for different Nexpansions (details see function)
    %******** 1. use a contrast map in which only values at the vertices are interpolated
    contrast_map_boundary_ex      = bilinear_interpolation(inv_SV_mat, 'favor_high_contrasts','vertices_only');
    
    %******** 2. set a threshold to contrast map and expand resulting boundaries by Nexpansions pixels in each direction (a range of Nexpansion values are possible)
    contrast_thresh     = params.local_params.soft_bound.contrast_thresh;      
    expVec              = params.local_params.soft_bound.start_end_iteration;   
    [regions_est,~,~,~] = RepBoundaryExp(contrast_map_boundary_ex,contrast_thresh,expVec);
    
    %******** 3. get resulting glimpse map and assign individual labels to the leftover (single) pixels
    glimpse_map         = get_glimpse_map(regions_est);
    [ro,co]             = find(glimpse_map==0);  % find all zeros, these are glimpses that consist of only one T-F unit   
    newIDs              = [max(glimpse_map(:)):1:(max(glimpse_map(:))+length(ro)-1)]; % give each t-f-unit that is labeled with zero, an individual integer as label
    for idx = 1:length(newIDs) 
        glimpse_map(ro(idx),co(idx)) = newIDs(idx);
    end
    contours            = get_contours(glimpse_map);
    
elseif strcmp(params.local_params.contour_method,'gbSuperpixels')==1 % use graph-based superpixels (see function for more details)
    [contours,glimpse_map, regions_est] = gb_superpixels(contrast_map, params.local_params.gb_superpixels.tau); 

elseif strcmp(params.local_params.contour_method,'slic_superpixels')==1 % slic_superpixels, gradient-ascent-based method (see function for more details)
    [contours, glimpse_map, ~]  = slic_superpixels(contrast_map, params.local_params.slic_superpixels.slic_compactness, params.local_params.slic_superpixels.slic_stepsize);
    regions_est                 = glimpse_map;

else
    error('undefined contour option')
end

%% **************** make sure that the glimpse IDs are a continuous row of integers **************
uniqueIDs   = unique(glimpse_map);
oldIDs      = sort(uniqueIDs,'ascend');  % these IDs potentially jump over some integers 

for bb = 1:numel(uniqueIDs)
    [rowb,colb]  = find(glimpse_map==oldIDs(bb));  % find elements of the glimpse
    for cc = 1:length(rowb)
        glimpse_map(rowb(cc),colb(cc))=bb;
    end
end

%% ****************
clustering_out.SV_mat       = SV_mat;               % matrix that hold the similarity values between beighboring t-f-units
clustering_out.contrast_map = contrast_map;         % SV_values converted into contrast values and interpolated vertices/t-f-units
clustering_out.regions_est  = regions_est;          % matrix of glimpse labels WITH sorrounding contours/boundaries
clustering_out.glimpse_map  = glimpse_map;          % matrix of glimpse labels only
clustering_out.contours     = contours;             % matrix with 1s to label contours

clustering_out.bin_mask     = bin_mask;
clustering_out.SV_vert      = SV_vert;
clustering_out.SV_horz      = SV_horz;
clustering_out.SV_types     = params.local_params.SV_types;


end




