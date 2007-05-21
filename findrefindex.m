function refindex=findrefindex(T,P_H2,P_He)
% NOT ITERATIVE -FINDS REFRACTIVE INDEX PROFILE
% P must be in bars, T in kelvin


Nr.H2=124.4.*P_H2.*(293./T);
Nr.He=35.83.*P_He.*(293./T);
Nr.tot=Nr.He+Nr.H2;
n=(Nr.tot./10^6)+1;
refindex=[n;1];
% last digit just a dummy for the raypath loop
