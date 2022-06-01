function [out_struct] = evaluate_contrasts(tru_contrast_map,est_contrast_map)
%EVALUATE_CONTRASTS compare estimated contrast map and ground truth
%contrast map, different metrics and ROC analysis
%   IN
%   tru_contrast_map : ground truth. dims(2*nFrames-1  X 2*nBands-1)
%   est_contrast_map : (2*nFrames-1  X 2*nBands-1)
%   OUT 
%   - mse: mean squared error
%   - rmse: root mean square error
%   - vaf: variance accounted for
%   - AROC: area under ROC curve
%   - ...

num_bands       = ceil(size(tru_contrast_map,2)/2);
nFrames         = ceil(size(tru_contrast_map,1)/2);
tru_contours    = (tru_contrast_map >= 0.5);

all_Contrast_vecs       = [];
all_c_vecs              = [];
est_contrasts_inBand    = zeros(nFrames-1, num_bands);
tru_contrasts_inBand    = zeros(nFrames-1, num_bands);
est_contrasts_acBand    = zeros(nFrames, num_bands-1);
tru_contrasts_acBand    = zeros(nFrames, num_bands-1);

for bb = 1:num_bands            % within filterbands
    inb_Contrast_vec            = est_contrast_map(2:2:end, bb*2-1);
    all_Contrast_vecs           = [all_Contrast_vecs; inb_Contrast_vec(:)];
    inb_c_vec                   = tru_contours(2:2:end,bb*2-1);   % before: ground_truth.ideal_contours
    all_c_vecs                  = [all_c_vecs; inb_c_vec(:)];
    est_contrasts_inBand(:,bb)  = est_contrast_map(2:2:end,bb*2-1);
    tru_contrasts_inBand(:,bb)  = tru_contrast_map(2:2:end,bb*2-1);
end
for gg = 1:num_bands-1          % across filterbands
    ac_Contrast_vec             = est_contrast_map(1:2:end, gg*2);
    all_Contrast_vecs           = [all_Contrast_vecs; ac_Contrast_vec(:)];
    ac_c_vec                    = tru_contours(1:2:end,gg*2);
    all_c_vecs                  = [all_c_vecs; ac_c_vec(:)];
    est_contrasts_acBand(:,gg)  = est_contrast_map(1:2:end, gg*2);
    tru_contrasts_acBand(:,gg)  = tru_contrast_map(1:2:end, gg*2);
end

all_contour_contrasts   = all_Contrast_vecs(find(all_c_vecs));
all_joint_contrasts     = all_Contrast_vecs(find(all_c_vecs*(-1)+1));
est_contrasts           = [est_contrasts_inBand(:) ; est_contrasts_acBand(:)];
tru_contrasts           = [tru_contrasts_inBand(:) ; tru_contrasts_acBand(:)];

%  metrics for evaluating continuous labels
mse    = mean(abs(tru_contrasts-est_contrasts).^2);         % mean squared error
mae    = mean(abs(tru_contrasts-est_contrasts));            % mean absolute error, mae, weighs all differences equally
rmse   = sqrt(mean(abs(tru_contrasts-est_contrasts).^2));   % (the rmse penalizes higher differences more than the mae)
vaf    = ((sum((tru_contrasts-mean(tru_contrasts)).*(est_contrasts-mean(est_contrasts))))/sqrt(sum((tru_contrasts-mean(tru_contrasts)).^2)*sum((est_contrasts-mean(est_contrasts)).^2)))^2;

[high_contrast_inds,~]  = find(tru_contrasts>=0.5);
tru_high_contrasts      = tru_contrasts(high_contrast_inds);
est_high_contrasts      = est_contrasts(high_contrast_inds);
rmset                   = sqrt(mean((tru_high_contrasts-est_high_contrasts).^2)); % rmse between the contrasts that are above 0.5

% receiver operating characteristics
ROC_data                = roc_curve(all_joint_contrasts, all_contour_contrasts, 0, 0);

% ***********
out_struct.mse      = mse;
out_struct.mae      = mae;
out_struct.rmse     = rmse;
out_struct.rmset    = rmset;
out_struct.vaf      = vaf;
out_struct.AROC     = ROC_data.param.AROC;
out_struct.ROC_data = ROC_data;

end

