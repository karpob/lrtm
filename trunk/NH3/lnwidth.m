function gammaNH3 = lnwidth(To,T,PH2,PHe,PNH3)
% To is reference Temperature-usually 300 k
% T is actual temperature
% Mol's is for extending functionality to use more than just H2,He,PH3
% dP is the partial pressure of each molecule contribution

% Calculates the linewidth for a gas mixture of mixtures (i) and produces
% a delta_nu line broadening in GHz for each line (j)
% delta_nu will be a vector of length (j)
% delta_nu(j)=sum{i}(delta_nu(i,j)*dP(i)*(To/T)^-(epsilson)
% where epsilon is (m+1)/2(m-1)   m is related to force law 
% m=3 dipole, m=infi is hardball

% CHANGE: ORDER MUST BE H2,He,NH3

% LOAD in CORRECT PRESSURE BROADENING FILES	

% Spilker Ammonia
expo1=2/3;
Tdiv=To/T;


load nh3lin gNH3o;				% NH3 line catalog  (only has REAL gNH3o for (1:200))
										%							^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


r=8.79*exp(-T/83);
gH2=2.34*(1-((2.122*exp(-T/116.8))/((exp(9.024-(T/20.3))-0.9918+PH2)^r)));
gHe=5.746-7.76*gH2+9.193*gH2^2-5.682*gH2^3+1.231*gH2^4;
gNH3=0.74*gNH3o;

gammaNH3=((gH2*PH2+gHe*PHe)*(Tdiv)^expo1)+gNH3*PNH3*(Tdiv);
