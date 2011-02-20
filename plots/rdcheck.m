X=rd(:,1);
Y=rd(:,2);
Z=rd(:,3);

[TH,PHI,R] = CART2SPH(X,Y,Z);

polar(PHI,R)