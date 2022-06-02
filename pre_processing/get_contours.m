function [contour_matrix] = get_contours(source_map)
%DERIVE CONTOURS BETWEEN GLIMPSES (LOCAL CLUSTERS)
%CONTOUR_MATRIX IS SET TO 1, IF TWO ADJACENT TIME-FREQUENCY UNITS BELONG TO
%DIFFERENT GLIMPSES AND ELSE 0
% 
% INPUTS:   - source_map: matrix (DIMS: num_win X nFilts) in which each
%             t-f-unit is labeled with the speaker ID (or signal component)
%             it belongs to. The labels are integers 1:num_components
% 
% OUTPUTS:  - contour_matrix: matrix (DIMS: num_win*2-1) X (nFilts*2-1)
%             the Input source_map is contoured by 1s in this
%             matrix. There are placeholders for each t-f-units, that are 0
contour_matrix = zeros(size(source_map,1)*2-1,size(source_map,2)*2-1);

% within filterband transitions
for ii = 1:size(source_map,1)-1    
    for jj = 1:size(source_map,2)       
        if source_map(ii,jj)==source_map(ii+1,jj)
            contour_matrix(ii*2,jj*2-1)=0;
        else 
            contour_matrix(ii*2,jj*2-1)=1;            
        end    
    end
end

% across filterband transitions
for kk = 1:size(source_map,1)    
    for ll = 1:size(source_map,2)-1               
        if source_map(kk,ll)==source_map(kk,ll+1)
            contour_matrix(kk*2-1,ll*2)=0;
        else 
            contour_matrix(kk*2-1,ll*2)=1;
        end        
    end    
end

% complete the contour matrix by adding ones or zeros in the corners (gaps between tiles and known contours)
% (elements that are not a direct transition between two t-f-units) the bin is set to 1 if 2 or more neighbours are also set to 1
for mm = 2:2:size(contour_matrix,1)-1
    for oo = 2:2:size(contour_matrix,2)-1 
        test_sum = contour_matrix(mm-1,oo)+contour_matrix(mm+1,oo)+contour_matrix(mm,oo-1)+contour_matrix(mm,oo+1);
        if test_sum >=2
            contour_matrix(mm,oo)=1;     
        elseif test_sum <=1
            contour_matrix(mm,oo)=0;
        end
    end    
end

end


