function [gammaNH3,psiNH3] = bgparametersnh3(To,T,PH2,PHe,PNH3)
% To is reference Temperature-usually 300 k
% T is actual temperature
% % dP is the partial pressure of each molecule contribution
% B&G

% Calculates the linewidth for a gas mixture of mixtures (i) and produces
% a delta_nu line broadening in GHz for each line (j)
% delta_nu will be a vector of length (j)
% delta_nu(j)=sum{i}(delta_nu(i,j)*dP(i)*(To/T)^-(epsilson)
% where epsilon is (m+1)/2(m-1)   m is related to force law 
% m=3 dipole, m=infi is hardball

% CHANGE: ORDER MUST BE H2,He,NH3

% LOAD in CORRECT PRESSURE BROADENING FILES	

% B&G Ammonia
% Redefine reference temperature since Poynter Kakar  use 295
%To=295;
To=300;

pow1=2/3;
Tdiv=To/T;


load nh3lin gNH3o;				% NH3 line catalog  (only has REAL gNH3o for (1:200))
										%							^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


%r=8.79*exp(-T/83);
%gH2=2.34*(1-((2.122*exp(-T/116.8))/((exp(9.024-(T/20.3))-0.9918+PH2)^r)));
%gHe=0.46+T/3000;
%gNH3=0.74*gNH3o;

gH2=2.318*PH2;
gHe=0.79*PHe;
gNH3=0.75*PNH3*gNH3o;

gammaNH3=((gH2+gHe)*((Tdiv)^pow1))+gNH3*(Tdiv);


%psiH2=5.746-7.76*gH2+9.193*gH2^2-5.682*gH2^3+1.231*gH2^4;
%psiHe=0.28-T/1750;
%psiNH3_only=0.50*gNH3o;


psiH2=1.92*PH2;
psiHe=0.3*PHe;
psiNH3_only=0.49*PNH3*gNH3o;
psiNH3=((psiH2+psiHe)*((Tdiv)^pow1))+psiNH3_only*Tdiv;

%psiNH3=(psiH2*PH2+psiHe*PHe)*(Tdiv)^pow1+psiNH3_only*PNH3*Tdiv;

