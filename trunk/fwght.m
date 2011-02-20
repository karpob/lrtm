function wght=fwght(layers,masterindex)

% function fwght
%
% This function checks to see that double contributions are included in the weighting functions
%


sl=size(layers,2);
for n=1:sl
   for k=1:max(masterindex(:,n))
      ch=masterindex(:,n)==masterindex(k,n);	% are there double contributions?
      %matchindextemp=masterindex(:).*ch  % Zero out other layer indices
      matchlayers=ch.*layers(:,n); % matched layer values
      wght(k,n)=sum(matchlayers);						% add those layers together 
   end
end

% Beamcoupling factor
return
N=size(wght,2);

num=sum((wght.*bw),2);
den=sum(N*1);

wf=num/den;

sP=size(P,1);
sa=zeros(sP-size(wf,1),1);

wfwa=[wfw;sa];
wfwb=[wfwb;sb];
wfwc=[wfwc;sc];
wfwd=[wfwd;sd];

weightingfactor=(wfwa+wfwb+wfwc+wfwd)/4;
