function [evaluation_data] = collect_eval_data_nSp(evaluated_stage, stimuli, params)
%COLLECT DATA SETS FOR EVALUATION 
% uses speech material from timit database to evaluate different stages of
% the speaker separation algorithm
% 
% IN
%   evaluated_stage - possibility to only evaluate the glimpse formation
%                     stage: 'local_groupin_only', or both glimpse
%                     formation and glimpse labeling:
%                     'local_and_global_grouping' (this is more time consuming)
%   stimuli         - struct that contains details of the stimuli that are
%                     used as evaluation set
%                     .mix_IDs = combination of timit speaker ID and sentence ID according to the table 'StimTableEvaluation.mat'
%                     .azimuths = vector with azimuth positions of each speaker in each mix
%                     .snr_vec = different snr values (snr between speech and noise)
%                     .nSP_vec = different speaker numbers
%   params          - struct that contains prep_params (parameters for the pre-processing), feature_params (parameters for the feature
%                     extraction), and local_params (parameters for the local clustering)
% 
% OUT  
%   evaluation_data ... see detailed descriptions at end of function


%% get model_name for the contrast estimation model (local contrast)
params.local_params.SV_types                = SVs_control_order(params.local_params.SV_types);
model_confs                                 = params.local_params.model_confs;
model_confs.SV_types                        = params.local_params.SV_types;
model_confs.params.feature_params           = params.feature_params;

hash_struct.feature_params                  = params.feature_params;
hash_struct.SV_types                        = params.local_params.SV_types;
hash_struct.model_confs                     = params.local_params.model_confs;
[model_name]                                = DataHash(hash_struct);



%% if the model exists, generate evaluation data (stimuli, features, ground_truth)
if isfile(['local_grouping',filesep,'contrast_estimation',filesep,'trained_models',filesep,model_name,'.mat']) == true | strcmp(params.local_params.model_confs.model_type,'none') == 1
    % *********************** settings and speech material for evaluated conditions
    snr_vec                 = stimuli.snr_vec;
    nSp_vec                 = stimuli.nSp_vec;
    mix_IDs                 = stimuli.mix_IDs;
    rooms_vec               = stimuli.rooms_vec;
    num_bands               = params.prep_params.nFilts; 
    mix_table               = load('StimTableEvaluation.mat');    % contains speaker and sentence IDs for the TIMIT mixtures
    SpeakerMix_table        = mix_table.StimTableEval;

    % *********************** load features for each mixID-angle-SNR combination from file and do local clustering ************ 
    % init result arrays
    all_contours_bacc               = [];
    filebased_acc                   = [];

    % contrast metrics
    filebased_rmse           = zeros(length(snr_vec),length(mix_IDs)); % root-mean-square error in each file
    filebased_vaf            = zeros(length(snr_vec),length(mix_IDs)); % variance accounted for
    filebased_mse            = zeros(length(snr_vec),length(mix_IDs));
    filebased_mae            = zeros(length(snr_vec),length(mix_IDs));
    
    filebased_ROC            = cell(length(snr_vec),length(mix_IDs));
    filebased_hit_fa         = zeros(length(snr_vec),length(mix_IDs));  % hit-false alarm rate after applying 0.5 as threshold to true&estimate
    filebased_per_corr       = zeros(length(snr_vec),length(mix_IDs));
    
    filebased_jaccard        = zeros(length(snr_vec),length(mix_IDs));
    filebased_labeled_acc    = zeros(length(snr_vec),length(mix_IDs));
    filebased_nEstGlimpses   = zeros(length(snr_vec),length(mix_IDs));
    filebased_nTrueGlimpses  = zeros(length(snr_vec),length(mix_IDs));
    filebased_nPixels        = zeros(length(snr_vec),length(mix_IDs));
    
    % IBM metrics
    filebased_pcorr                     = zeros(length(snr_vec),length(mix_IDs)); % percent correct of estimated source_map
    filebased_IBM_bacc                  = zeros(length(snr_vec),length(mix_IDs)); % accuracy for an IBM
    filebased_IBM_pcorr                 = zeros(length(snr_vec),length(mix_IDs));
    filebased_global_bestpcorr          = zeros(length(snr_vec),length(mix_IDs));
    filebased_target_rec                = zeros(length(snr_vec),length(mix_IDs));
    
    filebased_baseline_IBM_bacc     = zeros(length(snr_vec),length(mix_IDs));
    filebased_baseline_IBM_pcorr    = zeros(length(snr_vec),length(mix_IDs));
    filebased_baseline_target_rec   = zeros(length(snr_vec),length(mix_IDs));
    
    filebased_global_ROC            = cell(length(snr_vec),length(mix_IDs));
    
    all_inb_SV_vecs     = cell(num_bands,1);   % should have the size {number of filterbands} and for each band ((time_frames-1) or valid time frames)
    all_inb_c_vecs      = cell(num_bands,1);   % corresponding labels to the above cell array
    all_inb_contrasts   = cell(num_bands,1);
    all_ac_SV_vecs      = cell(num_bands-1,1); % should have the size {number of filterbands-1} and for each band ((time_frames) or valid time frames)
    all_ac_c_vecs       = cell(num_bands-1,1); % corresponding labels to the above cell array
    all_ac_contrasts    = cell(num_bands-1,1);
    
    
    
    % *********************** generate mix, pre-process, extract features, and ground truth
    for ff = 1:numel(rooms_vec)
        room.name = rooms_vec{ff};
        if strcmp(rooms_vec{ff},'anechoic') == 1
            room.name = rooms_vec{ff};
        elseif strcmp(rooms_vec{ff},'anechoic') == 0
            clear ans
            room_script = ['get_room_',rooms_vec{ff}];
            run(room_script);
            room = ans;         
        end
        
        for cc = 1:length(nSp_vec)
            nSp     = nSp_vec(cc);
            
            for bb = 1:length(snr_vec)
                SNR     = snr_vec(bb);
                
                for aa = 1:length(mix_IDs)
                    % ********** generate stimulus
                    mix_ID              = mix_IDs(aa);
                    tmpMix              = table2array(SpeakerMix_table(mix_ID,[2:(nSp+1),7:(nSp+6)]));
                    azmVec              = zeros(nSp,1);
                    for zz = 1:nSp
                        azmVec(zz,1) = stimuli.azimuths(aa,zz);
                    end
                    [stereo_mix,sig_components] = noisy_mix(nSp, azmVec, SNR, params.prep_params, tmpMix,'testing',room);
                    
                    % ********** pre-process mixture and extract features
                    [pre_processed]     = pre_process(stereo_mix, params.prep_params);
                    
                    feature_types{1,1}  = 'AZIMUTH';                feature_types{1,2}  = 'PERIODICITY';
                    [feature_mats]      = feature_ex(pre_processed, feature_types, params.feature_params, params.prep_params);
                    summed_power_mat    = feature_mats.power_mat;
                    
                    % ********** pre-process individual components and get ground truth (requires separate signal components)
                    [ground_truth]      = pre_process_ground_truth(sig_components, params.prep_params);
                    
                    % ********** do the local clustering/ glimpse estimation
                    [clustering_out]    = local_clustering(feature_mats, params);
                    
                    % ******** evaluate the output of local clustering as compared to the available ground_truth
                    for xx = 1:num_bands                                        % within filterband transitions
                        SV_vec = clustering_out.SV_mat(2:2:end,xx*2-1,:);
                        c_vec  = ground_truth.ideal_contours(2:2:end,xx*2-1);
                        contr  = ground_truth.ideal_contrasts(2:2:end,xx*2-1);
                        all_inb_SV_vecs{xx} = [all_inb_SV_vecs{xx}; SV_vec]; % collect data over all mix IDs
                        all_inb_c_vecs{xx}  = [all_inb_c_vecs{xx};c_vec];
                        all_inb_contrasts{xx}= [all_inb_contrasts{xx};contr];
                    end
                    for yy = 1:num_bands-1                                      % across filterband transitions
                        SV_vec = clustering_out.SV_mat(1:2:end,yy*2,:);
                        c_vec  = ground_truth.ideal_contours(1:2:end,yy*2);
                        contr  = ground_truth.ideal_contrasts(1:2:end,yy*2);
                        all_ac_SV_vecs{yy}  = [all_ac_SV_vecs{yy}; SV_vec];
                        all_ac_c_vecs{yy}   = [all_ac_c_vecs{yy}; c_vec];
                        all_ac_contrasts{yy}= [all_ac_contrasts{yy}; contr];
                    end
                    
                    % evaluate the binary labels/contours
                    [Cout_struct]               = evaluate_contours(ground_truth.ideal_contours,clustering_out.contours);
                    filebased_hit_fa(bb,aa)     = Cout_struct.hfar;
                    filebased_per_corr(bb,aa)   = Cout_struct.per_corr;
                    all_contours_bacc           = [all_contours_bacc ; Cout_struct.bacc];
                    filebased_acc               = [filebased_acc ; Cout_struct.acc];
                    
                    % evaluate the soft labels/contrasts
                    [out_struct]                = evaluate_contrasts(ground_truth.ideal_contrasts,clustering_out.contrast_map);
                    filebased_ROC{bb,aa}        = out_struct.ROC_data;
                    filebased_vaf(bb,aa)        = out_struct.vaf;
                    filebased_mse(bb,aa)        = out_struct.mse;
                    filebased_mae(bb,aa)        = out_struct.mae;
                    filebased_rmse(bb,aa)       = out_struct.rmse;
                    
                    % evaluate the glimpse_maps
                    glimpse_opts{1} = ground_truth.ideal_source_map;
                    glimpse_opts{2} = 1; % target source for labeled accuracy
                    [out_GL]        = evaluate_glimpses(ground_truth.ideal_glimpses,clustering_out.glimpse_map,{'MEAN_JACCARD','W_MEAN_JACCARD','LABELED_ACC'},glimpse_opts);
                    
                    filebased_jaccard(bb,aa)        = out_GL.mean_jaccard;
                    filebased_labeled_acc(bb,aa)    = out_GL.labeled_acc;
                    filebased_nEstGlimpses(bb,aa)   = out_GL.nEstGlimpses;
                    filebased_nTrueGlimpses(bb,aa)  = out_GL.nTrueGlimpses;
                    filebased_nPixels(bb,aa)        = out_GL.nPixels;
                    
                    
                    
                    if strcmp(evaluated_stage,'local_and_global_grouping')==1   % do global grouping WITH and WITHOUT prior local clustering (second version is the baseline system)
                        params.global_params.est_num_sources    = nSp+1;        % assuming the true number of sources (number of speakers + background noise)
                        
                        %% baseline system 1 (no glimpses, individual T-F units are assigned to azm peaks)
                        params_baseline                             = params;
                        params_baseline.global_params.feature_types = {'peak_azm_val','direct'};
                        params_baseline.global_params.group_by      = 'absolute_glimpse_features';
                        params_baseline.global_params.grouping_method = 'histogram1D';
                        
                        mini_glimpses                               = reshape([1:numel(summed_power_mat(:,:,1))],[size(summed_power_mat,1),size(summed_power_mat,2)]); % reshape a vector with integers from 1 to number of t-f-units
                        baseline_in.glimpse_map                     = mini_glimpses+1;      % in this glimpse map, each t-f-unit is its own glimpse
                        
                        [baseline_out]                              = global_grouping(baseline_in.glimpse_map, feature_mats, params_baseline);
                        [baseline_source_map, ~, ~, baseline_eval]  = match_comp_masks(baseline_out.est_source_map, ground_truth.ideal_source_map);
                        [baseline_rec_target]                    = sum_reconstructed_target(ground_truth.ideal_source_map, baseline_source_map, feature_mats.summed_power_mat);
                        
                        %% baseline system 2 (glimpses are assigned based on ground truth T-F unit labels) i.e. best possible glimpse assignment
                        [gt_glimpse_assignment, glimpseLabels]      = ground_truth_global( clustering_out.glimpse_map, ground_truth, params );
                        [~, ~, ~, git_global_out]                   = match_comp_masks( gt_glimpse_assignment, ground_truth.ideal_source_map);
                        
                        %% glimplse assignment based on feature estimates
                        [global_out]                                = global_grouping(clustering_out.glimpse_map, feature_mats, params);
                        
                        % averaging over all sources
                        [est_source_map,est_matched_IBMs,~,res_glob]= match_comp_masks(global_out.est_source_map, ground_truth.ideal_source_map);
                        [rec_target]                             = sum_reconstructed_target(ground_truth.ideal_source_map, est_source_map, feature_mats.summed_power_mat);
                        
                        % considering only a single target source
                        [glimpse_IBM1_bacc]                         = evaluate_IBM(ground_truth.IBMs(:,:,1),est_matched_IBMs(:,:,1));
                        
                        %% ROC analysis of glimpse similarity + clustering are used to achieve glimpse assignment
                        if strcmp(params.global_params.grouping_method,'linkage') == 1 || strcmp(params.global_params.grouping_method,'k_means') || strcmp(params.global_params.grouping_method,'k_means_weight') || strcmp(params.global_params.grouping_method,'k_meansSARI')
                            [glimpse_sizes]         = average_glimpse_features(feature_mats,clustering_out.glimpse_map,{'glimpse_size'},params);
                            big_glimpse_IDs         = find(glimpse_sizes{1,1}>=params.global_params.min_glimpse_size); % consider only glimpses that are above a minimum size
                            
                            num_glimpses            = numel(big_glimpse_IDs);
                            glimpsePairSims         = global_out.glimpse_pair_sim;   % similarity values of glimpse pairs
                            glimpsePairLabels       = zeros(((num_glimpses*(num_glimpses-1))/2),1);  % based on ground truth assignment of glimpses
                            counter                 = 1;
                            for ee = 1:(num_glimpses-1)
                                for dd = (ee+1):num_glimpses
                                    if glimpseLabels(big_glimpse_IDs(ee)) == glimpseLabels(big_glimpse_IDs(dd))
                                        glimpsePairLabels(counter,1) = 1;
                                    elseif glimpseLabels(big_glimpse_IDs(ee)) ~= glimpseLabels(big_glimpse_IDs(dd))
                                        glimpsePairLabels(counter,1) = 0;
                                    end
                                    counter = counter + 1;
                                end
                            end
                            allJointPairs       = glimpsePairSims(find(glimpsePairLabels));         % class with connected pairs
                            allSepPairs         = glimpsePairSims(find(glimpsePairLabels*(-1)+1));  % class with separated pairs
                            allJointPairs       = allJointPairs(1:100:end);       % these are way too many datapoints, so only take every 100 th point
                            allSepPairs         = allSepPairs(1:100:end);
                            
                            s_data              = unique(sort([allSepPairs; allJointPairs]));          % Sorted data points
                            s_data(isnan(s_data)) = [];                                 % Delete NaN values
                            d_data              = diff(s_data);                         % Difference between consecutive points
                            if(isempty(d_data))
                                ROCglobal.param.AROC = 0.5;
                            else
                                ROCglobal           = roc_curve(allSepPairs,allJointPairs, 0, 0);
                            end
                        else
                            ROCglobal.param.AROC = 0.5;   % no ROC analysis for setting 'histogram1D' as it doesn't use glimpse similarities
                        end
                        
                        %%
                        filebased_target_rec(bb,aa)         = rec_target;           % reconstructed target energy (when applying the estimated binary masks to the mixed signal)
                        filebased_pcorr(bb,aa)              = res_glob.per_corr;
                        filebased_IBM_bacc(bb,aa)           = glimpse_IBM1_bacc;  % accuracy for an IBM
                        filebased_IBM_pcorr(bb,aa)          = res_glob.per_corr(1);  % percent correct for the one
                        filebased_global_bestpcorr(bb,aa)   = git_global_out.per_corr;  % percent correct when applying the best possible labeling on est. glimpses
                        
                        filebased_baseline_IBM_pcorr(bb,aa) = baseline_eval.per_corr;
                        filebased_baseline_target_rec(bb,aa)= baseline_rec_target; %
                        filebased_global_ROC{bb,aa}         = ROCglobal;
                        
                    end
                    
                    
                    clear stereo_mix decomp_seg_sigs feature_mats pre_processed features_all ground_truth_data sig_components ground_truth
                    
                end
            end
        end
    end
    
% %         % *********************** generate mix, pre-process, extract features, and ground truth
% %     for cc = 1:length(nSp_vec)
% %         nSp     = nSp_vec(cc);
% %         
% %         for bb = 1:length(snr_vec)
% %             SNR     = snr_vec(bb);
% %             
% %             for aa = 1:length(mix_IDs)               
% %                 % ********** generate stimulus
% %                 mix_ID              = mix_IDs(aa);
% %                 tmpMix              = table2array(SpeakerMix_table(mix_ID,[2:(nSp+1),7:(nSp+6)]));
% %                 azmVec              = zeros(nSp,1);
% %                 for zz = 1:nSp
% %                     azmVec(zz,1) = stimuli.azimuths(aa,zz);
% %                 end
% %                 [stereo_mix,sig_components] = noisy_mix(nSp, azmVec, SNR, params.prep_params, tmpMix,'testing');
% %                 
% %                 % ********** pre-process mixture and extract features
% %                 [pre_processed]     = pre_process(stereo_mix, params.prep_params);
% % 
% %                 feature_types{1,1}  = 'AZIMUTH';                feature_types{1,2}  = 'PERIODICITY';                
% %                 [feature_mats]      = feature_ex(pre_processed, feature_types, params.feature_params, params.prep_params);
% %                 summed_power_mat    = feature_mats.power_mat;
% %                 
% %                 % ********** pre-process individual components and get ground truth (requires separate signal components)
% %                 [ground_truth]      = pre_process_ground_truth(sig_components, params.prep_params);
% %                 
% %                 % ********** do the local clustering/ glimpse estimation
% %                 [clustering_out]    = local_clustering(feature_mats, params);
% %                 
% %                 % ******** evaluate the output of local clustering as compared to the available ground_truth
% %                 for xx = 1:num_bands                                        % within filterband transitions
% %                     SV_vec = clustering_out.SV_mat(2:2:end,xx*2-1,:);
% %                     c_vec  = ground_truth.ideal_contours(2:2:end,xx*2-1);
% %                     contr  = ground_truth.ideal_contrasts(2:2:end,xx*2-1);
% %                     all_inb_SV_vecs{xx} = [all_inb_SV_vecs{xx}; SV_vec]; % collect data over all mix IDs
% %                     all_inb_c_vecs{xx}  = [all_inb_c_vecs{xx};c_vec];
% %                     all_inb_contrasts{xx}= [all_inb_contrasts{xx};contr];
% %                 end
% %                 for yy = 1:num_bands-1                                      % across filterband transitions
% %                     SV_vec = clustering_out.SV_mat(1:2:end,yy*2,:);
% %                     c_vec  = ground_truth.ideal_contours(1:2:end,yy*2);
% %                     contr  = ground_truth.ideal_contrasts(1:2:end,yy*2);
% %                     all_ac_SV_vecs{yy}  = [all_ac_SV_vecs{yy}; SV_vec];
% %                     all_ac_c_vecs{yy}   = [all_ac_c_vecs{yy}; c_vec];
% %                     all_ac_contrasts{yy}= [all_ac_contrasts{yy}; contr];
% %                 end
% %                 
% %                 % evaluate the binary labels/contours
% %                 [Cout_struct]               = evaluate_contours(ground_truth.ideal_contours,clustering_out.contours);
% %                 filebased_hit_fa(bb,aa)     = Cout_struct.hfar;
% %                 filebased_per_corr(bb,aa)   = Cout_struct.per_corr;               
% %                 all_contours_bacc           = [all_contours_bacc ; Cout_struct.bacc];
% %                 filebased_acc               = [filebased_acc ; Cout_struct.acc];                       
% %                 
% %                 % evaluate the soft labels/contrasts
% %                 [out_struct]                = evaluate_contrasts(ground_truth.ideal_contrasts,clustering_out.contrast_map);
% %                 filebased_ROC{bb,aa}        = out_struct.ROC_data;
% %                 filebased_vaf(bb,aa)        = out_struct.vaf;
% %                 filebased_mse(bb,aa)        = out_struct.mse;
% %                 filebased_mae(bb,aa)        = out_struct.mae;
% %                 filebased_rmse(bb,aa)       = out_struct.rmse;
% %                 
% %                 % evaluate the glimpse_maps
% %                 glimpse_opts{1} = ground_truth.ideal_source_map;
% %                 glimpse_opts{2} = 1; % target source for labeled accuracy
% %                 [out_GL]        = evaluate_glimpses(ground_truth.ideal_glimpses,clustering_out.glimpse_map,{'MEAN_JACCARD','W_MEAN_JACCARD','LABELED_ACC'},glimpse_opts);
% %                 
% %                 filebased_jaccard(bb,aa)        = out_GL.mean_jaccard;
% %                 filebased_labeled_acc(bb,aa)    = out_GL.labeled_acc;
% %                 filebased_nEstGlimpses(bb,aa)   = out_GL.nEstGlimpses;
% %                 filebased_nTrueGlimpses(bb,aa)  = out_GL.nTrueGlimpses;
% %                 filebased_nPixels(bb,aa)        = out_GL.nPixels;
% %                 
% %                 
% %                 
% %                 if strcmp(evaluated_stage,'local_and_global_grouping')==1   % do global grouping WITH and WITHOUT prior local clustering (second version is the baseline system)
% %                     params.global_params.est_num_sources    = nSp+1;        % assuming the true number of sources (number of speakers + background noise)
% %                                        
% %                     %% baseline system 1 (no glimpses, individual T-F units are assigned to azm peaks)
% %                     params_baseline                             = params;
% %                     params_baseline.global_params.feature_types = {'peak_azm_val','direct'};
% %                     params_baseline.global_params.group_by      = 'absolute_glimpse_features';
% %                     params_baseline.global_params.grouping_method = 'histogram1D';
% %                     
% %                     mini_glimpses                               = reshape([1:numel(summed_power_mat(:,:,1))],[size(summed_power_mat,1),size(summed_power_mat,2)]); % reshape a vector with integers from 1 to number of t-f-units
% %                     baseline_in.glimpse_map                     = mini_glimpses+1;      % in this glimpse map, each t-f-unit is its own glimpse
% % 
% %                     [baseline_out]                              = global_grouping(baseline_in.glimpse_map, feature_mats, params_baseline);
% %                     [baseline_source_map, ~, ~, baseline_eval]  = match_comp_masks(baseline_out.est_source_map, ground_truth.ideal_source_map);
% %                     [baseline_rec_target, ~]                    = sum_reconstructed_target(ground_truth.ideal_source_map, baseline_source_map, feature_mats.summed_power_mat);
% %                   
% %                     %% baseline system 2 (glimpses are assigned based on ground truth T-F unit labels) i.e. best possible glimpse assignment
% %                     [gt_glimpse_assignment, glimpseLabels]      = ground_truth_global( clustering_out.glimpse_map, ground_truth, params );
% %                     [~, ~, ~, git_global_out]                   = match_comp_masks( gt_glimpse_assignment, ground_truth.ideal_source_map);
% %                     
% %                     %% glimplse assignment based on feature estimates 
% %                     [global_out]                                = global_grouping(clustering_out.glimpse_map, feature_mats, params);
% %                     
% %                     % averaging over all sources
% %                     [est_source_map,est_matched_IBMs,~,res_glob]= match_comp_masks(global_out.est_source_map, ground_truth.ideal_source_map);
% %                     [rec_target, ~]                             = sum_reconstructed_target(ground_truth.ideal_source_map, est_source_map, feature_mats.summed_power_mat);
% %                     
% %                     % considering only a single target source
% %                     [glimpse_IBM1_bacc]                         = evaluate_IBM(ground_truth.IBMs(:,:,1),est_matched_IBMs(:,:,1));
% %                     
% %                     %% ROC analysis of glimpse similarity + clustering are used to achieve glimpse assignment                    
% %                     if strcmp(params.global_params.grouping_method,'linkage') == 1 || strcmp(params.global_params.grouping_method,'k_means') || strcmp(params.global_params.grouping_method,'k_means_weight') || strcmp(params.global_params.grouping_method,'k_meansSARI')
% %                         [glimpse_sizes]         = average_glimpse_features(feature_mats,clustering_out.glimpse_map,{'glimpse_size'},params);
% %                         big_glimpse_IDs         = find(glimpse_sizes{1,1}>=params.global_params.min_glimpse_size); % consider only glimpses that are above a minimum size
% %                         
% %                         num_glimpses            = numel(big_glimpse_IDs);
% %                         glimpsePairSims         = global_out.glimpse_pair_sim;   % similarity values of glimpse pairs
% %                         glimpsePairLabels       = zeros(((num_glimpses*(num_glimpses-1))/2),1);  % based on ground truth assignment of glimpses
% %                         counter                 = 1;
% %                         for ee = 1:(num_glimpses-1)
% %                             for dd = (ee+1):num_glimpses
% %                                 if glimpseLabels(big_glimpse_IDs(ee)) == glimpseLabels(big_glimpse_IDs(dd))
% %                                     glimpsePairLabels(counter,1) = 1;
% %                                 elseif glimpseLabels(big_glimpse_IDs(ee)) ~= glimpseLabels(big_glimpse_IDs(dd))
% %                                     glimpsePairLabels(counter,1) = 0;
% %                                 end
% %                                 counter = counter + 1;
% %                             end
% %                         end
% %                         allJointPairs       = glimpsePairSims(find(glimpsePairLabels));         % class with connected pairs
% %                         allSepPairs         = glimpsePairSims(find(glimpsePairLabels*(-1)+1));  % class with separated pairs
% %                         allJointPairs       = allJointPairs(1:100:end);       % these are way too many datapoints, so only take every 100 th point
% %                         allSepPairs         = allSepPairs(1:100:end);
% %                         
% %                         s_data              = unique(sort([allSepPairs; allJointPairs]));          % Sorted data points
% %                         s_data(isnan(s_data)) = [];                                 % Delete NaN values
% %                         d_data              = diff(s_data);                         % Difference between consecutive points
% %                         if(isempty(d_data))
% %                             ROCglobal.param.AROC = 0.5;
% %                         else
% %                             ROCglobal           = roc_curve(allSepPairs,allJointPairs, 0, 0);
% %                         end
% %                     else
% %                         ROCglobal.param.AROC = 0.5;   % no ROC analysis for setting 'histogram1D' as it doesn't use glimpse similarities                        
% %                     end
% %                     
% %                     %% 
% %                     filebased_target_rec(bb,aa)         = rec_target;           % reconstructed target energy (when applying the estimated binary masks to the mixed signal)
% %                     filebased_pcorr(bb,aa)              = res_glob.per_corr;
% %                     filebased_IBM_bacc(bb,aa)           = glimpse_IBM1_bacc;  % accuracy for an IBM
% %                     filebased_IBM_pcorr(bb,aa)          = res_glob.per_corr(1);  % percent correct for the one
% %                     filebased_global_bestpcorr(bb,aa)   = git_global_out.per_corr;  % percent correct when applying the best possible labeling on est. glimpses
% %                     
% %                     filebased_baseline_IBM_pcorr(bb,aa) = baseline_eval.per_corr;
% %                     filebased_baseline_target_rec(bb,aa)= baseline_rec_target; %
% %                     filebased_global_ROC{bb,aa}         = ROCglobal;
% %                     
% %                 end
% %                 
% %                 
% %                 clear stereo_mix decomp_seg_sigs feature_mats pre_processed features_all ground_truth_data sig_components ground_truth
% %                 
% %             end
% %         end
% %     end

    evaluation_data.model_name                  = clustering_out.model_name;
    evaluation_data.azimuths                    = stimuli.azimuths;
    
    % metrics as summed over all datapoints in all files
    evaluation_data.overall_inb_SV_vecs{1}      = all_inb_SV_vecs;                 % all estimated within band similarity values
    evaluation_data.overall_inb_c_vecs{1}       = all_inb_c_vecs;                  % all estimated within band contours
    evaluation_data.overall_inb_contrasts{1}    = all_inb_contrasts;               % all estimated within band contrasts
    evaluation_data.overall_ac_SV_vecs{1}       = all_ac_SV_vecs;                  % all estimated across band similarity values
    evaluation_data.overall_ac_c_vecs{1}        = all_ac_c_vecs;                   % all estimated across band contours
    evaluation_data.overall_ac_contrasts{1}     = all_ac_contrasts;                % all estimated across band contrasts
    
    % metrics as summed only over datapoints in individual separate files    
    evaluation_data.filebased_bacc              = all_contours_bacc;              % balanced accuracy of estimated contours
    evaluation_data.filebased_acc               = filebased_acc;                  % accuracy of estimated contours
    evaluation_data.filebased_rmse              = filebased_rmse;                 % root mean squared error of estimated contrasts
    evaluation_data.filebased_vaf               = filebased_vaf;                  % variance accounted for of estimated contrasts
    evaluation_data.filebased_ROC               = filebased_ROC;                  % receiver operating characteristics estimated contrasts&ideal contours
    evaluation_data.filebased_mae               = filebased_mae;                  % mean absolute error
    evaluation_data.filebased_mse               = filebased_mse;                  % mean squared error
    evaluation_data.filebased_hfa_rate          = filebased_hit_fa;               % hit minus false alarm rate
    evaluation_data.filebased_per_corr          = filebased_per_corr;             % percent correct of estimated contours
    
    evaluation_data.filebased_jaccard           = filebased_jaccard;              % jaccard index of estimated/ideal glimpse maps
    evaluation_data.filebased_labeled_acc       = filebased_labeled_acc;          % accuracy of an IBM with best possible labeling of the glimpses
    evaluation_data.filebased_nEstGlimpses      = filebased_nEstGlimpses(bb,aa);  % number of estimated glimpses
    evaluation_data.filebased_nTrueGlimpses     = filebased_nTrueGlimpses(bb,aa); % number of true glimpses
    evaluation_data.filebased_nPixels           = filebased_nPixels(bb,aa);       % number of TF units
    
    if strcmp(evaluated_stage,'local_and_global_grouping')==1
        
        evaluation_data.filebased_IBM_bacc                  = filebased_IBM_bacc; % balanced accuracy for each estimated IBM
        evaluation_data.filebased_IBM_pcorr                 = filebased_IBM_pcorr; % percent correct for each estimated IBM
        evaluation_data.filebased_global_bestpcorr          = filebased_global_bestpcorr;
        evaluation_data.filebased_target_rec                = filebased_target_rec; % reconstructed target energy (when applying the estimated binary masks to the mixed signal)
        evaluation_data.filebased_global_ROC                = filebased_global_ROC; % receiver operating characteristics estimated glimpse similarity & ideal glimpse pair label
        evaluation_data.filebased_global_pcorr              = filebased_pcorr;      % overall percent correctly labeled T-F units
        
        evaluation_data.filebased_baseline_IBM_bacc         = filebased_baseline_IBM_bacc;  % balanced accuracy for baseline IBM (no glimpse formation, just spatial assigment of T-F units)
        evaluation_data.filebased_baseline_IBM_pcorr        = filebased_baseline_IBM_pcorr; % percent correct for baseline IBM (no glimpse formation, just spatial assigment of T-F units)
        evaluation_data.filebased_baseline_target_rec       = filebased_baseline_target_rec;% reconstructed target energy for baseline IBM (no glimpse formation, just spatial assigment of T-F units)
    end

%% if the model doesn't exist, skip this condition, set all results to zero
elseif isfile(['local_grouping',filesep,'contrast_estimation',filesep,'trained_models',filesep,model_name,'.mat']) == false && strcmp(params.local_params.model_confs.model_type,'none') ~= 1
    disp('local model does not exist, condition is skipped')
    evaluation_data                         = [];
    evaluation_data.model_name              = model_name;
    evaluation_data.azimuths                = stimuli.azimuths;
    evaluation_data.nSp                     = stimuli.nSp_vec;        

    evaluation_data.overall_inb_SV_vecs{1}  = 0;                evaluation_data.overall_inb_c_vecs{1}   = 0;
    evaluation_data.overall_inb_contrasts{1}= 0;                evaluation_data.overall_ac_SV_vecs{1}   = 0;
    evaluation_data.overall_ac_c_vecs{1}    = 0;                evaluation_data.overall_ac_contrasts{1} = 0;
    
    evaluation_data.filebased_bacc          = 0;                evaluation_data.filebased_acc           = 0;
    evaluation_data.filebased_rmse          = 0;                evaluation_data.filebased_vaf           = 0;
    evaluation_data.filebased_mae           = 0;                evaluation_data.filebased_mse           = 0;
    evaluation_data.filebased_hfa_rate      = 0;                evaluation_data.filebased_per_corr      = 0; 
                    
    evaluation_data.filebased_ROC{1,1}.param.AROC = 0;
    evaluation_data.filebased_jaccard           = 0;
    evaluation_data.filebased_labeled_acc       = 0;
    evaluation_data.filebased_nEstGlimpses      = 0;
    evaluation_data.filebased_nTrueGlimpses     = 0;
    evaluation_data.filebased_nPixels           = 0;

    if strcmp(evaluated_stage,'local_and_global_grouping')==1

        evaluation_data.filebased_IBM_bacc                  = 0;
        evaluation_data.filebased_IBM_pcorr                 = 0; %
        evaluation_data.filebased_global_bestpcorr          = 0;
        evaluation_data.filebased_target_rec                = 0;
        evaluation_data.filebased_global_ROC                = 0;
        evaluation_data.filebased_global_pcorr              = 0;
        
        evaluation_data.filebased_baseline_IBM_bacc         = 0;
        evaluation_data.filebased_baseline_IBM_pcorr        = 0;
        evaluation_data.filebased_baseline_target_rec       = 0;
        
        
    end    

end

end