function wght=wghtcalc(alpha,hght,costheta);% Revised 10/17/00
%The indexing is: i=1 is the TOP of the atmosphere and its tau must be set to zero
%so that the first REAL layer has a tau above it of zero.
i=size(alpha,1);


tauatm=alpha.*hght./costheta;
trans=1-exp(-tauatm);
for i=1:i;
   tauabove(i)=sum(tauatm(1:i));
   wght(i+1)=trans(i+1).*exp(-tauabove(i));
end
