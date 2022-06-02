% process_individual_sentence_mix
% 
% 
% 1) Settings 
% 2) generate stimulus 
% 3) formation of glimpses, binary mask estimation
% 4a) if ground truth (i.e. separate speaker signals) is not available, 
%     plot features and estimated contrasts and binary masks
% 4b) if ground truth is available, compare the estimated masks and contrasts 
%     to ground representations


clear 
close all

plot_ground_truth = input('Are separate signals available? y/n  : ','s');  % if this is available, the estimated contrasts can be compared, if not, just estimate contrasts


%% 1) SETTINGS
% model settings for contrast estimation
[eval_settings, var_names] = load_model_settings('model_settings'); % txt file with model specific parameter sets (see evaluation>specs_eval_conditions)
eval_settings              = eval_settings(1,:);       
% eval_settings(1,:) -> best performing model in paper
% eval_settings(2,:) -> contrast estimation based on single relative relative feature (here: correlation of azimuth probability) without prior training


num_tested_models       = size(eval_settings,1);       % number of models that are supposed to be testet
new_settings            = cell(2,size(eval_settings,2));
new_settings(1,:)       = var_names;
new_settings(2,:)       = eval_settings;

% define independent parameters
p_prep          = def_prep_params;
p_feat          = def_feature_params;
p_local         = def_local_params;
p_global        = def_global_params; 

% initialize and update parameter structs
params          = update_parameters_new(new_settings, p_prep, p_feat, p_local, p_global);

%% 2) stimulus 
numSpeakers = 2;    % up to 4 speakers
azmVec      = [-90; -60; 10; 0; 10; 90]; % speaker positions in azimuth plane
SNR         = 5;  % SNR single speaker to diffuse background noise
params.global_params.est_num_sources = numSpeakers+1;

room.name   = 'anechoic'; 
htc{1}      = false;    % option to add harmonic tone complex = false
[noisy_mix, sig_components]   = noisy_mix(numSpeakers,azmVec,SNR, params.prep_params,[],'example',room,htc);



%% 3) glimpse extraction, binary mask estimation
[main_out]          = spatial_glimpses_main(noisy_mix,params);

% ******************** baseline system (without previous glimpse formation) ******************** 
both_ears_power_mat = mean(main_out.pre_processed.power_mat,3);     % power mat of both ears (average over both ears)
int_vec             = 1:numel(both_ears_power_mat);    % vector with integers from 1 to total number of t-f-units
mini_glimpses       = reshape(int_vec,[size(both_ears_power_mat,1),size(both_ears_power_mat,2)]);
[global_out_nogl]   = global_grouping(mini_glimpses, main_out.feature_mats, params); % global out without glimpse formation


%% 4) plots
both_ears_power_mat         = mean(main_out.pre_processed.power_mat,3);     % power mat averaged over both ears (win_length X num_windows X filter_channels)
t_end                       = length(noisy_mix)*1/params.prep_params.fs;
time_vec                    = linspace(0,t_end);
freq_vec                    = linspace(1,params.prep_params.fs/2);
bands                       = size(both_ears_power_mat,2);

%% 4a) if separate signals are not available
if strcmp(plot_ground_truth,'n')==1
    figure(1)
    nSp         = params.global_params.est_num_sources;
    common_min  = min(min(10*log10(both_ears_power_mat)));
    common_max  = max(max(10*log10(both_ears_power_mat)));
    
    subplot(3,2,1) % cochleagram of mixture
    imagesc(time_vec, freq_vec, 10*log10(both_ears_power_mat).',[common_min common_max])
    colormap(gca,'parula')
    set(gca,'YDir','normal');   colorbar
    xlabel('time [ms]');        ylabel('freq. [Hz]')
    title('speaker mix, average over both ears')
    
    subplot(3,2,2)
    imagesc(time_vec,1:bands,main_out.clustering_out.contrast_map.')
    set(gca,'YDir','normal');   title(['estimated contrasts'],'FontSize',15)
    xlabel('Time (s)','FontSize',14); ylabel('filter band','FontSize',14)
    colormap(gca,flipud(gray))
    
    subplot(3,2,3)
    bar(main_out.global_out.feature_rates(:,1),main_out.global_out.feature_rates(:,2))
    title('estimated source locations','FontSize',15)
    
    subplot(3,2,4)
    Cmethod = params.local_params.contour_method;
    imagesc(time_vec,1:bands,main_out.clustering_out.contours.')
    set(gca,'YDir','normal');  title(['estimated contours using ', Cmethod],'FontSize',15)
    xlabel('Time (s)','FontSize',14);         ylabel('filter band','FontSize',14)
    colormap(gca,flipud(gray))
    
    subplot(3,2,6)
    imagesc(time_vec,1:bands,main_out.global_out.est_source_map.',[1 nSp])
    set(gca,'YDir','normal'); title('estimated source map','FontSize',15)
    xlabel('Time (s)','FontSize',14);         ylabel('filter band','FontSize',14)
    colormap(gca,'pink')
    
    
%% 4b) if separate signals are available   
elseif strcmp(plot_ground_truth,'y')==1   % evaluations and plots that are possible, if separate signal components are available
    [ground_truth]              = pre_process_ground_truth(sig_components, params.prep_params);
    
    % ideal labeling of estimated glimpses
    params.global_params.est_num_sources    = max(ground_truth.ideal_source_map(:)); % infer the true number of sources
    [gt_source_map, glimpseLabels]          = ground_truth_global( main_out.clustering_out.glimpse_map, ground_truth, params );
    [gt_source_map, ~, ~, ~]                = match_comp_masks( gt_source_map, ground_truth.ideal_source_map);
    
    [main_out.global_out.est_source_map, main_out.global_out.matched_est_IBMs, ~, evalWithGlimpse]         = match_comp_masks( main_out.global_out.est_source_map, ground_truth.ideal_source_map);
    [nogl_source_map, ~, ~, ~]              = match_comp_masks( global_out_nogl.est_source_map, ground_truth.ideal_source_map);
    
    [contrast_eval]                         = evaluate_contrasts(ground_truth.ideal_contrasts,main_out.clustering_out.contrast_map);
    
    %% ************* cochleagrams of single sources and mix
    figure(1)
    nSp         = params.global_params.est_num_sources;
    common_min  = min(min(10*log10(both_ears_power_mat)));
    common_max  = max(max(10*log10(both_ears_power_mat)));
    
    
    subplot((ceil(nSp/2)+1),2,1) % cochleagram of mixture
    imagesc(time_vec, freq_vec, 10*log10(both_ears_power_mat).',[common_min common_max])
    set(gca,'YDir','normal');   colorbar
    xlabel('time [ms]');        ylabel('freq. [Hz]')
    title('speaker mix, average over both ears')
    
    subplot((ceil(nSp/2)+1),2,2) % mixture with ideal glimpse contours
    small_mat       = 10*log10(both_ears_power_mat);
    large_plot_mat  = zeros(size(small_mat,1)*4,size(small_mat,2)*4);
    for aa = 1:4:size(large_plot_mat,1)
        for bb = 1:4:size(large_plot_mat,2)
            large_plot_mat(aa:(aa+3),bb:(bb+3))= repmat(small_mat((aa+3)/4,(bb+3)/4),4,4);
        end
    end
    
    temp_contours   = get_contours(ground_truth.ideal_source_map);
    large_contours  = zeros(size(large_plot_mat));
    for aa = 2:2:(size(large_contours,1)-1)
        for bb = 2:2:(size(large_contours,2)-1)
            large_contours((aa:aa+1),(bb:bb+1))=repmat((temp_contours(aa/2,bb/2)),2,2);
        end
    end
     
    
    im1 = imagesc(large_plot_mat.',[common_min common_max]);
    hold on
    colormap(gca,'parula')
    [xx,yy] = find(large_contours);
    plot(xx,yy,'ks','MarkerSize',2,'MarkerFaceColor', 'k')
    set(gca,'YDir','normal');   colorbar
    xlabel('time [ms]');        ylabel('freq. [Hz]')
    title('speaker mix and ideal glimpse contours')
    hold off
    
    for aa = 1:nSp      % cochleagrams of individual speakers
        subplot((ceil(nSp/2)+1),2,(aa+2))
        summed_sp = sum(ground_truth.power_mats{1,aa},3);
        
        imagesc(time_vec, freq_vec, 10*log10(summed_sp).',[common_min common_max])
        set(gca,'YDir','normal');   colorbar
        xlabel('time [ms]');        ylabel('freq. [Hz]')
        title(['speaker ',num2str(aa),', average over both ears'])
    end

    
    %%
    
    figure(2)
    common_scale = max(ground_truth.ideal_source_map(:));
    
    subplot(5,2,1)
    imagesc(time_vec,1:bands,(ground_truth.ideal_contrasts).')
    set(gca,'YDir','normal');   title('ideal contrasts','FontSize',15)
    xlabel('Time (s)','FontSize',14);ylabel('filter band','FontSize',14)
    
    colormap(flipud(gray))
    subplot(5,2,2)
    imagesc(time_vec,1:bands,main_out.clustering_out.contrast_map.')
    set(gca,'YDir','normal');   title(['estimated contrasts'],'FontSize',15)
    xlabel('Time (s)','FontSize',14); ylabel('filter band','FontSize',14)
    
    subplot(5,2,3)
    imagesc(time_vec,1:bands,ground_truth.ideal_contours.')
    set(gca,'YDir','normal');   title('ideal contours','FontSize',15)
    xlabel('Time (s)','FontSize',14);         ylabel('filter band','FontSize',14)
    
    subplot(5,2,4)
    Cmethod = params.local_params.contour_method;
    imagesc(time_vec,1:bands,main_out.clustering_out.contours.')
    set(gca,'YDir','normal');  title(['estimated contours using ', Cmethod],'FontSize',15)
    xlabel('Time (s)','FontSize',14);         ylabel('filter band','FontSize',14)
    
    subplot(5,2,5)
    imagesc(time_vec,1:bands,ground_truth.ideal_source_map.',[1 common_scale])
    set(gca,'YDir','normal');   title('ideal source map','FontSize',15)
    xlabel('Time (s)','FontSize',14);         ylabel('filter band','FontSize',14)
    colormap(gca,'pink')
    
    subplot(5,2,6)
    imagesc(time_vec,1:bands,main_out.global_out.est_source_map.',[1 common_scale])
    set(gca,'YDir','normal'); title('estimated source map','FontSize',15)
    xlabel('Time (s)','FontSize',14);         ylabel('filter band','FontSize',14)
    colormap(gca,'pink')
    
    subplot(5,2,7)
    imagesc(time_vec,1:bands,nogl_source_map.',[1 common_scale])
    set(gca,'YDir','normal'); title('estimated source map without glimpse formation','FontSize',15)
    xlabel('Time (s)','FontSize',14);         ylabel('filter band','FontSize',14)
    colormap(gca,'pink')
    
    subplot(5,2,8)
    bar(main_out.global_out.feature_rates(:,1),main_out.global_out.feature_rates(:,2))
    title('estimated source locations','FontSize',15)    
    
    subplot(5,2,9)
    imagesc(time_vec, 1:bands, gt_source_map.',[1 common_scale])
    set(gca,'YDir','normal')
    colormap(gca,'pink')
    title('estimated glimses with ground truth labeling','FontSize',15)
    xlabel('Time (s)','FontSize',14);         ylabel('filter band','FontSize',14)
    
end




%%  feature maps
figure(3)
azm_vec = -90:5:90;
num_frames = size(main_out.feature_mats.complete_NAC_mat,2);
tau_vec     = params.feature_params.minLags:params.feature_params.maxLags;
f0_vec      = params.prep_params.fs./tau_vec;

subplot(3,1,1)  % illustrate pitch salience
imagesc(time_vec, 1:bands, (squeeze(max(main_out.feature_mats.complete_NAC_mat,[],1))).',[0.5 1.1])
set(gca,'YDir','normal');   colorbar
xlabel('time [ms]');        ylabel('filterband')
title('max autocorrelation value per t-f-unit (pitch salience)')

subplot(3,1,2)  % illustrate azimuth estimates
[~,azm_I] = max(main_out.feature_mats.log_az_probs,[],1);
imagesc(azm_vec(squeeze(azm_I)).',[-90 90])
set(gca,'YDir','normal');   colorbar
xlabel('time [ms]');        ylabel('filterband')
title('azimuth with max probability per t-f-unit')

subplot(3,1,3)  % illustrate pitch estimates
f0_I = zeros(num_frames,bands);
for aa= 1:num_frames
    for bb = 1:bands
        [~,b,~,~] = findLocalPeaks(main_out.feature_mats.complete_NAC_mat(:,aa,bb),1);
        if isempty(b)==1
            [~,b,~,~] = findLocalPeaks(main_out.feature_mats.complete_NAC_mat(:,aa,bb),2);
        end
        f0_I(aa,bb) = b;
    end
end
imagesc(f0_vec(squeeze(f0_I)).',[100 300])
set(gca,'YDir','normal');   colorbar
xlabel('time [ms]');        ylabel('filterband')
title('f_0 with max autocorrelation per t-f-unit')







