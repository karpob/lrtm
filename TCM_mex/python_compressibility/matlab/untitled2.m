% test integration
R=8.314472/2.01594;
To=298.15
Po=0.101325
cp_normal_h2_ideal(0.000001)
%syms T
quad(@cp_normal_h2_ideal,0.0001,To)
%int(R*cp_normal_h2_idel(T)/T,T,0,To)