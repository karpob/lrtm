function refindex=findrefindex(T,P_H2,P_He,P_CH4,P_H2O,reference_select)
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
%                             ->P_He: Partial pressure of Helium (He) in bars
%                             ->P_CH4: Partial pressure of Methane (CH4) in bars 
%                             -P_H2O; Partial pressure of Water vapor (H2O) in bars
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
	Nr.H2=0.0.*P_H2.*(293./T); Should give a n=1.0 (N=0)
        Nr.He=0.0.*P_He.*(293./T);
end

if (reference_select==1)
	Nr.H2=124.4.*P_H2.*(293./T); % Note Hoffman truncated/rounded to the nearest tenth
	Nr.He=35.83.*P_He.*(293./T); % Both Equations from DeBoer's thesis
end

if (reference_select==2)
	Nr.H2=136.*(P_H2./1.01325).*(273./T); % From Essen
	Nr.He=35.*(P_He./1.01325).*(273./T); % From Essen
end
if (reference_select==3)
        Nr.H2=136.*(P_H2./1.01325).*(273./T); % From Essen
        Nr.He=35.*(P_He./1.01325).*(273./T); % From Essen
        Nr.CH4=440.*(P_CH4./1.01325).*(273./T); % From Spilker's Thesis 1990
        Nr.H2O=3.73e5.*(P_H2O.*1000).*(1./(T.*T)); % From Thayer, 1974
%      Add some other cool/hot gases at some point
end
if(reference_select<3)
    Nr.tot=Nr.He+Nr.H2;
else
    Nr.tot=Nr.He+Nr.H2+Nr.CH4+Nr.H2O;
end
n=(Nr.tot./10^6)+1;
refindex=[n;1];
% last digit just a dummy for the raypath loop where pressure is zero. ie.
% a place for cold space P=0 T=2.7 K
