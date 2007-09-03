X=intercept(1:461,1);
Y=intercept(1:461,2);
Z=intercept(1:461,3);

[TH,PHI,R] = CART2SPH(X,Y,Z);

polar(PHI,R,'r')