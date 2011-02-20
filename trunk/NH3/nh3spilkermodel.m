function [alphanh3]=nh3spilkermodel(f,T,P,H2mr,Hemr,NH3mr) % Spilker model

%input f as a vector in GHz, T in Kelvin, P in Bars, H2mr, Hemr and NH3mr in mole
%fraction, so for a 1% mixture of NH3, NH3mr = 0.01, etc.  f is the only
%variable that can be a vector.  Opacity (alphanh3) is returned as a vector matching f
%in dB/km.

PH2=H2mr*P/1.01325; %convert to atm
PHe=Hemr*P/1.01325;
PNH3=NH3mr*P/1.01325;

load nh3lincatPK16 % Load line catalog from Poynter and Kakar (1975)
gammaNH3o=Width*760/1000;
fo=Fo/1000;
x=2.98.*J.*(J+1);
y=1.09.*(K.^2);
z=4.8./T;
d=exp(z.*(y-x));

l=(2.*J+1).*(K.^2);
m=J.*(J+1);
n=(l./m);
q=PNH3./(T.^(3.5));

A =1230.*n.*(fo.^2).*S.*q.*d;

To=300;
xi=2/3;
xi2=1;
Tdiv=To/T;

r=8.79*exp(-T/83);
v=2.122*exp(-T/116.8);
p=(exp(9.024-(T/20.3)))-0.9918+PH2;
w=(v/(p^r));
gnu1=2.34*(1-w);

%gnu1=2.34*(1-((2.122*exp(-T/116.8))/((exp(9.024-(T/20.3))-0.9918+PH2)^r)));
gnu2=0.46+(T/3000);
gnu3=0.74;
gH2=gnu1*PH2;
gHe=gnu2*PHe;
gNH3=gnu3*PNH3*gammaNH3o;
gamma=((gH2+gHe)*((Tdiv)^xi))+gNH3*(300/T)^(xi2);

znu1=5.7465-7.7644*gnu1+(9.1931*(gnu1^2))-(5.6816*(gnu1^3))+(1.2307*(gnu1^4));
znu2=0.28-(T/1750);
znu3=0.50;
zH2=znu1*PH2;
zHe=znu2*PHe;
zNH3=znu3*PNH3*gammaNH3o;
zeta=((zH2+zHe)*((Tdiv)^xi))+zNH3*(300/T)^(xi2);


n=size(f,2);  %returns the number of columns in f
m=size(fo,1); %returns the number of rows in fo
% f1 f2 f3 f4 ....fn  n times where n is the number of frequency steps
% f1 f2 f3 f4 ....fn				in the observation range                            
% ...
% f1 f2 f3 f4 ....fn 
% m times where m is the number of spectral lines
zetasize=size(fo,1);
pst=-0.45*PNH3*ones(zetasize,1);	

nones=ones(1,n);
mones=ones(m,1); 
f_matrix=mones*f;
fo_matrix=fo*nones;
dnu_matrix=gamma*nones;    
ce_matrix=zeta*nones;
pst_matrix=pst*nones;

Aa=(2)*((f_matrix./fo_matrix).^2);			

Bb=(dnu_matrix-ce_matrix).*(f_matrix.^2);
Cc=dnu_matrix+ce_matrix;
Dd=((fo_matrix+pst_matrix).^2) + (dnu_matrix.^2)-(ce_matrix.^2);
Ee=f_matrix.^2;
Jj=(fo_matrix+pst_matrix).^2;
Gg=dnu_matrix.^2;
Hh=ce_matrix.^2;
Ii=4*(f_matrix.^2).*(dnu_matrix.^2);
Ff=Aa.*(Bb+(Cc.*Dd))./(((Ee-Jj-Gg+Hh).^2)+Ii);						% m x n matrix

Fbr=Ff;

n=size(f,2);
nones=ones(1,n);
Tb=180;
if T<Tb
C=(Tb^2/70600-0.337)+(1/110.4-Tb/35300)*T;
else    
C=-0.337+(T/110.4)-((T.^2)/70600);
end
A_matrix = A*nones;
br_alpha_matrix=A_matrix.*Fbr.*C;
alpha_opdep=sum(br_alpha_matrix,1);
alphanh3=alpha_opdep*434294.5;