function decay=ph3decay(P,xPH3)

ortonPt=[0.15,0.645]';
%ortonxPH3=[4.1e-7,1.2e-5]';
ortonxPH3=[4.1e-7,max(xPH3)]';
y=log(ortonxPH3);
x=log(ortonPt);

slope=(y(2)-y(1))/(x(2)-x(1));
b=y(1)-slope.*x(1);
% Now have the equation for the line


sp=size(P,1);
Ps=P(2:sp);					% strip the 'fake' zero from P
sph3=size(xPH3,1);
xPH3s=xPH3(2:sph3);		% strip the 'fake' zero from P

% Start deep, stop high

Pdecay_high=ortonPt(1);			% Pressure of highest altitude given by Orton 2000 
Pdecay_deep=ortonPt(2);			% Pressure (bars) where ph3 begins to decay
pph3_high=ortonxPH3(1);			% Orton
pph3_deep=ortonxPH3(2);			% Orton

Pstartindex=sum(Pdecay_deep>=Ps);		% Gives index (P) of smallest P at or before Pdecay
Pstopindex=sum(Pdecay_high>=Ps);

%logPstart=log10(Ps(Pstartindex));
%logPstop=log10(Ps(Pstopindex));


% want to linear interp, but in the loglog of P and mixing ratio
lnP=log(Ps);
lnxPH3=log(xPH3s);

%%%%%%%%%%%%%%%%%%%%%%%%EXAMPLE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tab =
%    1950    150.697
%    1960    179.323
%    1970    203.212
%    1980    226.505
%    1990    249.633

%then the population in 1975, obtained by table lookup within the matrix tab, is
% p = interp1(tab(:,1),tab(:,2),1975)
% p =
%     214.8585
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

%testlnPt=linspace(min(lnPt),max(lnPt),10)
%testPt=exp(testlnPt);
ln_interp_xPH3=slope.*lnP(1:Pstartindex)+b;
interp_xPH3=exp(ln_interp_xPH3);
decay=xPH3s;
decay(1:Pstartindex)=interp_xPH3;
decay=[0;decay];			% reinsert the zero pressure condition
%figure

%semilogy(ortonxPH3,ortonPt,'h')
%hold on
%semilogy(decay,P,'r.')

%figure
%loglog(ortonxPH3,ortonPt,'h')
%hold on
%loglog(decay,P,'r.')

