function [b1,b2,b3,b1w,b2w,b3w]=beamsample(N_ring_one,BWHM)

% For generating beam samples with 2pi symm in theta, with beamwidth of phi
% Nphi=3 !!! You're stuck with this, unless you change dimensions of Tb to
% include phi ie., get rid of Ta,Tb,Tc,Td convention in maintamone
Nphi=3; %Don't touch this unless you want a wrong answer.
dphi_degree=BWHM/Nphi;
phi_degree=cumsum(dphi_degree.*ones(Nphi,1))
n=phi_degree./dphi_degree % radius multiple from first ring
Ntheta=N_ring_one.*(2*(n)-1); % number of samples for each ring

% center phi on z-axis (around pi/2)

phi_degree_norm=90+phi_degree;
phi=phi_degree_norm*pi/180; % convert to radian

lasttheta1=2*pi-2*pi/Ntheta(1);
theta1=linspace(0,lasttheta1,Ntheta(1));
lasttheta2=2*pi-2*pi/Ntheta(2);
theta2=linspace(0,lasttheta2,Ntheta(2));
lasttheta3=2*pi-2*pi/Ntheta(3);
theta3=linspace(0,lasttheta3,Ntheta(3));
phi1=phi(1);
phi2=phi(2);
phi3=phi(3);

r=1;
[R1,TH1,PHI1]=meshgrid(r,theta1,phi1);
[R2,TH2,PHI2]=meshgrid(r,theta2,phi2);
[R3,TH3,PHI3]=meshgrid(r,theta3,phi3);

[x1,y1,z1]=sph2cart(TH1,PHI1,R1);
[x2,y2,z2]=sph2cart(TH2,PHI2,R2);
[x3,y3,z3]=sph2cart(TH3,PHI3,R3);

x1=x1(:);
y1=y1(:);
z1=z1(:);

x2=x2(:);
y2=y2(:);
z2=z2(:);

x3=x3(:);
y3=y3(:);
z3=z3(:);

b1=[x1';y1';z1'];
b2=[x2';y2';z2'];
b3=[x3';y3';z3'];

% NOW DO BEAMWEIGHTING
deltabeam=dphi_degree/2	;	% since phi is sampled around by 2pi, the deltaphi (from axis is1/2)
delta(1)=deltabeam;
delta(2)=deltabeam*2;
delta(3)=deltabeam*3;
wt=exp(-2.76*(delta./BWHM).^2);
b1w=wt(1);
b2w=wt(2);
b3w=wt(3);

