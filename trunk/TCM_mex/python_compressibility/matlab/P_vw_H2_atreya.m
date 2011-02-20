function Pvw=Pvw_H2_atreya(P,T)
R=8.314472*10;
a=0.2453;
b=0.0266;
rho=P/(R*T);
Pvw=(R*T)/((1/rho)-b)-a*rho^2;
