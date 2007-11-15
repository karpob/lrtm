function refindex=findrefindex(T,P_H2,P_He)
%
% function findrefindex
%
% Calculate the refractive index of each layer. Note refractive index only accounts for H2 and He.
% 
%             VARIABLE DEFINITIONS:
%
%                        -->  INPUT:
%
%                             ->T: temperature of the layer in K
%                             ->P_H2: Partial pressure of Hydrogen (H2) in bars
%                             ->P_He: Partial pressure of Helium in bars
%
%                        <-- OUTPUT:
%
%                              <--Refractive index of the layer


%JPH NOT ITERATIVE -FINDS REFRACTIVE INDEX PROFILE
%JPH P must be in bars, T in kelvin


Nr.H2=124.4.*P_H2.*(293./T);
Nr.He=35.83.*P_He.*(293./T);
Nr.tot=Nr.He+Nr.H2;
n=(Nr.tot./10^6)+1;
refindex=[n;1];
% last digit just a dummy for the raypath loop
