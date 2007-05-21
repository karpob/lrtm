function ph3decay=ph3photolysis(P,xPH3)


sp=size(P,1);
Ps=P(2:sp);					% strip the 'fake' zero from P
sph3=size(xPH3,1);
xPH3s=xPH3(2:sph3);		% strip the 'fake' zero from P

% Start deep, stop high

Pdecay_high=0.15;			% Pressure of highest altitude given by Orton 2000 
Pdecay_deep=0.645;		% Pressure (bars) where ph3 begins to photolyse
pph3_high=4.1e-7;			% Orton
pph3_deep=12e-6;			% Orton

Pstartindex=sum(Pdecay_deep>=Ps);		% Gives index (P) of smallest P at or before Pdecay
Pstopindex=sum(Pdecay_high>=Ps);

logPstart=log10(Ps(Pstartindex));
logPstop=
