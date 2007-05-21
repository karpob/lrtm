X=internormal(:,1);
Y=internormal(:,2);
Z=internormal(:,3);

[TH,PHI,R] = CART2SPH(X,Y,Z);
R=R./(max(max(R)));
polar(PHI,R,'r')