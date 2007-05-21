function alpha=goodman_model(f,P,T,XH2O)
OpticaldepthstodB=434294.5;
nu=((f*1e9)/3e8)/100;

front_part=XH2O*P*((273/T)^(13/3))*nu^2;

frac_minus_074=(delta_nu_1)/((nu-0.74)^2 +delta_nu_1^2);
frac_plus_074=(delta_nu_1)/((nu+0.74)^2 +delta_nu_1^2);

middle_part=(1.073e-8)*(frac_minus_074+frac_plus_074);

end_part=(17.20e-8)*delta_nu_1;

alpha=OpticaldepthstodB*front_part*middle_part+end_part;
