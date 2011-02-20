function alphah2=falphah2_goodman_by_joiner(f,T,PH2,PHe)


if PH2==0
   alphah2=0;
   return
end


atmperbar=0.987;
% Change the pressures from bars to atm for this
PH2=PH2*atmperbar;
PHe=PHe*atmperbar;


c=2.99792458e10;
f=f.*1e9;			% Hz
lambda=c./f;		% Change into wavelength (cm)
lambdasq=1./lambda.^2;



coef=lambdasq*PH2*(4e-11);

A=PH2*(273/T)^(2.8);		% First part
B=PHe*1.7*(273/T)^(2.61);
alphah2=coef*(A+B);
