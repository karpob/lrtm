function Z=Z_H2_atreya(rho,T)
R=8.314472/100;
a=0.2453;
b=0.0266;
Vm=1/rho;
Z=1/(1-(b/Vm))-a/(R*T*Vm);