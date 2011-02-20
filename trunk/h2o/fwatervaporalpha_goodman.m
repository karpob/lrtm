function alpha=fwatervaporalpha_goodman(f,P,T,XH2O,XH2,XHe)


nu=((f*1e9)/3e8)/100;
PH2O=XH2O*P*750.06;
P=P*750.06;

%stuff from equation B13
front_part=PH2O*((273/T)^(13/3))*nu^2;

delta_nu_1=0.1*((P/760)*(273/T)^(2/3))*(0.810*XH2+0.35*XHe);

frac_minus_074=(delta_nu_1)/((nu-0.74)^2 +delta_nu_1^2);

frac_plus_074=(delta_nu_1)/((nu+0.74)^2 +delta_nu_1^2);

middle_part=(1.073e-8)*(frac_minus_074+frac_plus_074);

end_part=(17.20e-8)*delta_nu_1;

alpha=front_part*(middle_part+end_part);

%alpha=(1.75e-8)*PH2O*((273/T)^5)*(P/760)*nu^2;