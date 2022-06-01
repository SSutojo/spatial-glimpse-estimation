function Rxy = pearsons_corr(vec1, vec2)
% CALCULATE PEARSONS CORRELATION COEFFICIENT BETWEEN TWO DATASETS OF SAME SIZE
% (vectors of same length)
% Pearson's correlation coefficient
% INPUTS:   vec1 - vector with dataset 1
%           vec2 - vector with dataset 2 (same length as vec1)
% OUTPUTS:  Rxy - correlation coefficient between vec1 and vec2 (scalar)



% make sure both inputs are column vectors
vec1 = vec1(:);
vec2 = vec2(:);

if mean(vec1 == vec2)==1        % this is important in case some t-f-units are set to silent or to zero, in that case
    Rxy = 1;                    % the calculation of the correlation as below would cause an error
                                % in this the correlation is 1 because vec1 and vec2 are the same (most likely they are both set to zero)    
elseif mean(vec1 == vec2)~=1 && sum(vec1-repmat(mean(vec1),length(vec1),1))==0 || mean(vec1 == vec2)~=1 && sum(vec2-repmat(mean(vec2),length(vec2),1))==0
    Rxy = 0;                    % if only one of the two t-f-units is set to 0 the calculation of the Rxy would also cause an error
                                % it can be assumed that the correlation is very low    
elseif mean(vec1==vec2)~=1 && sum(vec1-repmat(mean(vec1),length(vec1),1))~=0 && sum(vec2-repmat(mean(vec2),length(vec2),1))~=0    
                                % if vec1 and vec2 are different (at least some elements are different) and neither of the t-f-units is 
                                % set to zero. Pearson correlation coeff. can be calculated
    R_numerator = sum((vec1-mean(vec1)).*(vec2-mean(vec2)));
    R_denominator = sqrt(sum((vec1-mean(vec1)).^2))*sqrt(sum((vec2-mean(vec2)).^2));
%     *****alternative Formulierung 
%     R_numerator = sum(vec1.*vec2)-length(vec1)*mean(vec1)*mean(vec2);
%     R_denominator = sqrt(sum(vec1.^2)-length(vec1)*mean(vec1)^2)*sqrt(sum(vec2.^2)-length(vec2)*mean(vec2)^2);    
    Rxy = R_numerator/R_denominator;    
else
    error('unknown case for correlation coefficient')
end  

end