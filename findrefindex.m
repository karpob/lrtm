function refindex=findrefindex(T,P_H2,P_He,reference_select)
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
%                             ->reference_select: Select reference source
%                                               0=No refraction ie. n=1.0D0
%                                               1=Original DeBoer/Hoffman He/H2
%                                               2=Karpowicz He/H2
%                                               3=Karpowicz Everything (He,H2,NH3, etc..)
%
%                        <-- OUTPUT:
%
%                              <--Refractive index of the layer


%JPH NOT ITERATIVE -FINDS REFRACTIVE INDEX PROFILE
%JPH P must be in bars, T in kelvin

if (reference_select==0)
	Nr.H2=0.0.*P_H2.*(293./T);
        Nr.He=0.0.*P_He.*(293./T);
end

if (reference_select==1)
	Nr.H2=124.43.*P_H2.*(293./T);
	Nr.He=35.83.*P_He.*(293./T);
end

if (reference_select==2)
	Nr.H2=136.*(P_H2./1.01325).*(273./T);
	Nr.He=35.*(P_He./1.01325).*(273./T);
end
if (reference_select==3)
        Nr.H2=136.*(P_H2./1.01325).*(273./T);
        Nr.He=35.*(P_He./1.01325).*(273./T);
%      Add some other cool/hot gases at some point
end

Nr.tot=Nr.He+Nr.H2;
n=(Nr.tot./10^6)+1;
refindex=[n;1];
% last digit just a dummy for the raypath loop
