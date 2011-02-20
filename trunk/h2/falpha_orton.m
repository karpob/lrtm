function alphah2=falpha_orton(f,T,P)
if(P<1e-6)
    alphah2=0;
else
GHztoHz=1e9;
f=f*GHztoHz;
c=2.99792458e10; %cm/sec
nu=f./c;

Lo=2.68719e19; %Loschmidt number molecules/cm^3 at stp

[tmp1,tmp2]=h2h2_mod(T,nu,nu,1);
spec=tmp1(1);
freq=tmp2(1);

Pascal_per_bar=100000;
RH2=4124.18; %Nm/Kg/K
cm3_per_m3=1e6;
rho=(Pascal_per_bar*P)/(RH2*T);
rho=rho/cm3_per_m3;
grams_per_kg=1000;
rho=rho*grams_per_kg;
grams_per_mole=2*1.007822;
molecules_per_mole=6.0221415e23;
rho=rho/grams_per_mole;
rho=rho*molecules_per_mole;
amagat=rho/Lo;
alphah2=spec(1)*amagat*amagat;
end;

