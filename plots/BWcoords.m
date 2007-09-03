function [Xbeam, Ybeam, Zbeam]=BWcoords(Raydirection,Rayorigin,deltaradian)

X=Raydirection(1);
Y=Raydirection(2);
Z=Raydirection(3);


% Generate FOUR alternate 'lasers' off Beam center

% init
Xbeam=zeros(4,1);
Ybeam=Xbeam;
Zbeam=Ybeam;
%

% Convert to spherical  (Reminder Matlab defines:THETA is the angle in the xy plane counterclockwise from the
% positive x axis.  PHI is the elevation of the direction vector from the xy plane 

[TH,PHI,R]=cart2sph(X,Y,Z);
TH_low=TH-deltaradian;
TH_hi=TH+deltaradian;
PHI_low=PHI-deltaradian;
PHI_hi=PHI+deltaradian;
[Xbeam(1) Ybeam(1) Zbeam(1)]=sph2cart(TH_low,PHI,R);
[Xbeam(2) Ybeam(2) Zbeam(2)]=sph2cart(TH_hi,PHI,R);
[Xbeam(3) Ybeam(3) Zbeam(3)]=sph2cart(TH,PHI_low,R);
[Xbeam(4) Ybeam(4) Zbeam(4)]=sph2cart(TH,PHI_hi,R);
