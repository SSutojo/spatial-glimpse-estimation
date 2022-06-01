load logtest
data=gmmsamp(mix,1e4);

types=cell(1,4);
[types{:}]=deal('spherical','diag','full','ppca');

options=foptions;
options(1)=0;  % 1 to see progress in fits

for i1=1:length(types)
  types{i1}
  mix1=gmm(mix.nin,mix.ncentres,types{i1});
  mix1=gmminit(mix1,data,options);
  mixl=mix1;
  mix1=gmmem(mix1,data,options);
  mixl=loggmmem(mixl,data,options);
  err=[ mean(abs(mix1.priors-mixl.priors))...
	mean(abs(mix1.centres(:) - mixl.centres(:)))...
	mean(abs(mix1.covars(:) - mixl.covars(:)))];
  if i1==4 %ppca
    err=[err mean(abs( mix1.U(:) - mixl.U(:)))...
	 mean(abs( mix1.lambda - mixl.lambda))];
  end % if i1
  err
end % for i1


  