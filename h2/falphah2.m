function alphah2=falphah2(f,T,PH2,PHe,PCH4)

% function alphah2=falphah2(f,T,PH2,PHe,PCH4)
% ((1.545*10^-6)./(lambda.^2)).*PH2[PH2....

if PH2==0
   alphah2=0;
   return
end


atmperbar=0.987;
% Change the pressures from bars to atm for this
PH2=PH2*atmperbar;
PHe=PHe*atmperbar;
PCH4=PCH4*atmperbar;


f=f.*1e9;			% Hz
lambda=(3e10)./f;		% Change into frequency (cm)
lambdasq=lambda.^2;
A=((3.557e-11)./lambdasq).*PH2;		% First part
B=PH2*((273/T)^3.12)+(1.382*PHe*((273/T)^2.24))+(9.322*PCH4*((273/T)^3.34));
alphah2=A.*B;
