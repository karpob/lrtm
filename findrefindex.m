function refindex=findrefindex(T,P_H2,P_He,P_CH4,P_H2O,D_SOL,f,reference_select)
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
%                             ->P_H2O: Partial pressure of Water vapor (H2O) in bars
%                             ->D_sol: Cloud density of the solution cloud H20-NH3 in g/cm^3
%                             ->f: Frequency in GHz
%                             ->reference_select: Select reference source
%                                               0=No refraction ie. n=1.0D0
%                                               1=Original DeBoer/Hoffman He/H2
%                                               2=Karpowicz He/H2
%                                               3=Karpowicz (H2,He,CH4,H2O,No Clouds) 
%                                               4=Karpowicz (H2,He,CH4,H2O, Solution Cloud)
%
%                        <-- OUTPUT:
%
%                              <--Refractive index of the layer


%JPH NOT ITERATIVE -FINDS REFRACTIVE INDEX PROFILE
%JPH P must be in bars, T in kelvin

if (reference_select==0)
	Nr.H2=0.0.*P_H2.*(293./T); %Should give a n=1.0 (N=0)
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
        cd refractivity
        Nr.H2O=refractivity_water_vapor_rueger(P_H2O,T); % 
        cd ..
end
if (reference_select==4)
        Nr.H2=136.*(P_H2./1.01325).*(273./T); % From Essen
        Nr.He=35.*(P_He./1.01325).*(273./T); % From Essen
        Nr.CH4=440.*(P_CH4./1.01325).*(273./T); % From Spilker's Thesis 1990
        cd refractivity
        Nr.H2O=refractivity_water_vapor_rueger(P_H2O,T); % From Rueger,2002
        Nr.Solution_Cloud=-1*refractivity_solution_cloud_93(T,D_SOL,1,f); % Minus 1 necessary because Liebe is interested in Delay, leaves out the - in Re(-K)
        cd ..
end
if(reference_select<3)
    Nr.tot=Nr.He+Nr.H2;
end
if(reference_select==3)
    Nr.tot=Nr.He+Nr.H2+Nr.CH4+Nr.H2O;
end
if(reference_select==4)
    Nr.tot=Nr.He+Nr.H2+Nr.CH4+Nr.H2O+Nr.Solution_Cloud;
end


n=(Nr.tot./10^6)+1;
refindex=n;
