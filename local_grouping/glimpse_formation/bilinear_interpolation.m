function [interp_contrast_map] = bilinear_interpolation(contrast_map, weighting_option, interp_mode)
%DO A BILINEAR INTERPOLATION FOR AN EMPTY ELEMENT IN A 2 DIMENSIONAL MATRIX
% INPUTS:  - x1, x2 :           the two neighbouring values in x-direction
%          - y1, y2 :           the two neighbouring values in y-direction
% OUTPUTS: - interpolated_val:  the result of the interpolation
%
% it's assumed that the interpolated value has the same distance
% to all 4 input values (so they are all weighted equally)
% 
% [interpolated_val] = bilinear_interp(x1,x2,y1,y2)
% interpolated_val = (x1 +x2 +y1 +y2 )/4;
% 
% interpolate 'vertices_only', 'vertices_and_TFunits'



%% interpolate the vertices
if islogical(contrast_map)~=1               % interpolate conintuous values
    if strcmp(weighting_option,'none') == 1 % alle Richtungen gleichberechtigt
        for mm = 2:2:size(contrast_map,1)-1         % time frames
            for oo = 2:2:size(contrast_map,2)-1 % filterbands
                test_sum = contrast_map(mm-1,oo)+contrast_map(mm+1,oo)+contrast_map(mm,oo-1)+contrast_map(mm,oo+1);
                contrast_map(mm,oo) = test_sum/4;
            end
        end
        
    elseif strcmp(weighting_option,'favor_high_contrasts') == 1 || strcmp(weighting_option,'high_low') == 1  || strcmp(weighting_option,'high_high') == 1
        for mm = 2:2:size(contrast_map,1)-1         % time frames
            for oo = 2:2:size(contrast_map,2)-1 % filterbands
                test_sum1 = contrast_map(mm-1,oo)+contrast_map(mm+1,oo);
                test_sum2 = contrast_map(mm,oo-1)+contrast_map(mm,oo+1);
                if test_sum1 >= test_sum2
                    contrast_map(mm,oo) = test_sum1/2;  % use the mean of the lower SV values (this is more likely a contour)
                elseif test_sum2 > test_sum1
                    contrast_map(mm,oo) = test_sum2/2;
                end
            end
        end
        
    elseif strcmp(weighting_option,'favor_low_contrasts') == 1   % rather connect = use the lower contrast
        for mm = 2:2:size(contrast_map,1)-1         % time frames
            for oo = 2:2:size(contrast_map,2)-1 % filterbands
                test_sum1 = contrast_map(mm-1,oo)+contrast_map(mm+1,oo);
                test_sum2 = contrast_map(mm,oo-1)+contrast_map(mm,oo+1);
                if test_sum1 >= test_sum2
                    contrast_map(mm,oo) = test_sum2/2;  % use the mean of the lower SV values (this is more likely a contour)
                elseif test_sum2 > test_sum1
                    contrast_map(mm,oo) = test_sum1/2;
                end
            end
        end
 
    else
        error('define a weighting option')
    end
    
elseif islogical(contrast_map)==1
    % alle Richtungen gleichberechtigt, bzw. default 'none', or favor high contrasts for interpolating the vertices and low contrasts for T-F units
    if strcmp(weighting_option,'none') == 1 || strcmp(weighting_option,'favor_high_contrasts')==1 || strcmp(weighting_option,'high_low') == 1  || strcmp(weighting_option,'high_high') == 1
        
        for mm = 2:2:size(contrast_map,1)-1         % time frames
            for oo = 2:2:size(contrast_map,2)-1 % filterbands
                test_sum = contrast_map(mm-1,oo)+contrast_map(mm+1,oo)+contrast_map(mm,oo-1)+contrast_map(mm,oo+1);
                if test_sum >= 2
                    contrast_map(mm,oo) = 1;  % use the mean of the lower SV values (this is more likely a contour)
                elseif test_sum < 2
                    contrast_map(mm,oo) = 0;
                end
            end
        end
        
    elseif strcmp(weighting_option,'favor_low_contrasts') == 1   % rather connect = use the lower contrast
        for mm = 2:2:size(contrast_map,1)-1         % time frames
            for oo = 2:2:size(contrast_map,2)-1 % filterbands
                test_sum = contrast_map(mm-1,oo)+contrast_map(mm+1,oo)+contrast_map(mm,oo-1)+contrast_map(mm,oo+1);
                if test_sum3 >= 3
                    contrast_map(mm,oo) = 1;  % use the mean of the lower SV values (this is more likely a contour)
                elseif test_sum < 3
                    contrast_map(mm,oo) = 0;
                end
            end
        end
        
    else
        error('define a weighting option')
    end
    
end

halp_point = 1;

%% optionally interpolage the placeholders for t-f-units as well

if strcmp(interp_mode,'vertices_only') == 1  % do not interpolate the placeholders for the t-f-units (used when regionsgrow is applied later)
    interp_contrast_map = contrast_map;

elseif strcmp(interp_mode,'vertices_and_TFunits') == 1  % interpolate the t-f-units as well
    
    if strcmp(weighting_option,'high_low') == 1
        % favor low contrasts when interpolating the t-f units
        for mm = 3:2:size(contrast_map,1)-1         % time frames
            for oo = 3:2:size(contrast_map,2)-1 % filterbands
                test_sum1 = contrast_map(mm-1,oo)+contrast_map(mm+1,oo);
                test_sum2 = contrast_map(mm,oo-1)+contrast_map(mm,oo+1);
                if test_sum1 >= test_sum2
                    contrast_map(mm,oo) = test_sum2/2;  % use the mean of the lower SV values (this is more likely a contour)
                elseif test_sum2 > test_sum1
                    contrast_map(mm,oo) = test_sum1/2;
                end
            end
        end
        
    elseif strcmp(weighting_option,'high_high') == 1
        % favor low contrasts when interpolating the t-f units
        for mm = 3:2:size(contrast_map,1)-1         % time frames
            for oo = 3:2:size(contrast_map,2)-1 % filterbands
                test_sum1 = contrast_map(mm-1,oo)+contrast_map(mm+1,oo);
                test_sum2 = contrast_map(mm,oo-1)+contrast_map(mm,oo+1);
                if test_sum1 >= test_sum2
                    contrast_map(mm,oo) = test_sum1/2;  % use the mean of the lower SV values (this is more likely a contour)
                elseif test_sum2 > test_sum1
                    contrast_map(mm,oo) = test_sum2/2;
                end
            end
        end
        
    else
        % the default setting is to not favor any direction
        for mm = 3:2:size(contrast_map,1)-1     % time frames
            for oo = 3:2:size(contrast_map,2)-1 % filterbands
                test_sum = contrast_map(mm-1,oo)+contrast_map(mm+1,oo)+contrast_map(mm,oo-1)+contrast_map(mm,oo+1);
                contrast_map(mm,oo) = test_sum/4;                
%                 test_sum = [contrast_map(mm-1,oo);contrast_map(mm+1,oo);contrast_map(mm,oo-1);contrast_map(mm,oo+1)];
%                 contrast_map(mm,oo) = median(test_sum);  % using the median instead of mean, smoothes the edges of large enough glimpses
            end
        end
        
    end
    
    
    % for first/last filterband/frame, don't favor any direction
    for aa = 3:2:(size(contrast_map,2)-2)
        contrast_map(1,aa)  = (contrast_map(1,aa-1)+contrast_map(1,aa+1)+contrast_map(2,aa))/3;
        contrast_map(end,aa)= (contrast_map(end,aa-1)+contrast_map(end,aa+1)+contrast_map(end-1,aa))/3;
    end
    
    for aa = 3:2:(size(contrast_map,1)-2)
        contrast_map(aa,1)  = (contrast_map(aa-1,1)+contrast_map(aa+1,1)+contrast_map(aa,2))/3;
        contrast_map(aa,end)= (contrast_map(aa-1,end)+contrast_map(aa+1,end)+contrast_map(aa,end-1))/3;
    end
    
    contrast_map(1,1)       = (contrast_map(1,2) +contrast_map(2,1))/2;
    contrast_map(end,1)     = (contrast_map(end-1,1) +contrast_map(end,2))/2;
    contrast_map(1,end)     = (contrast_map(1,end-1) +contrast_map(2,end))/2;
    contrast_map(end,end)   = (contrast_map(end,end-1) +contrast_map(end-1,end))/2;
    interp_contrast_map     = contrast_map;
    
else
    error('define interpolation mode')
end


end