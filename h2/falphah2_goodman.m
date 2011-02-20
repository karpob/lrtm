function alphah2=falphah2_goodman(f,T,PH2,PHe)


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
nu=f./(c);		% Change into wavenumber (cm^-1)
nusq=nu.^2;

N_He=PHe*(273.15/T)*(2.687e19);
N_H2=PH2*(273.15/T)*(2.687e19);

coef=N_H2*nusq*(1/c)*1e-38;

A=0.377*N_H2*((T/100)^(-0.8));		% First part
B=0.535*N_He*((T/100)^(-0.61));
alphah2=coef*(A+B);
