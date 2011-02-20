load alphadata
f=alphadata(:,1)';
xH2=alphadata(1,2);
xHe=alphadata(1,3);
xPH3=alphadata(1,4);
P=alphadata(1,5);
T=alphadata(1,6);

vvw_alpha=alphashape(f,xH2,xHe,xPH3,P,T)