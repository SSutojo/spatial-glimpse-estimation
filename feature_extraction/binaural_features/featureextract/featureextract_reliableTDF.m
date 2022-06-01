function [ out ] = featureextract_reliableTDF(FEXT, TDFs, rHat)
%FEATUREEXTRACT_RELIABLETDF generate a mask that identifies reliable time
%domain features

out = rHat > FEXT.minXCreliable;
if FEXT.removeEndpoints
    out = out & (TDFs(:,:,1)<(FEXT.maxITDsec-.5/FEXT.fs));
    out = out & (TDFs(:,:,1)>-(FEXT.maxITDsec-.5/FEXT.fs));
    if ndims(TDFs)>2
        out = out & (TDFs(:,:,2)<FEXT.maxTDOAsec);
        out = out & (TDFs(:,:,2)>-FEXT.maxTDOAsec);
        %out = out & (TDFs(:,:,3)<FEXT.maxTDOAsec);
        %out = out & (TDFs(:,:,3)>-FEXT.maxTDOAsec);
    end        
end

end

