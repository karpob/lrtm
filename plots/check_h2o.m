%
clear
P=6.89;
T=373;
XH2=0.8467;
XHe=0.1333;
XH2O=0.02;
OpticaldepthstodB=434294.5;
P_H2O=P*XH2O;

P_H2=XH2*P;
P_He=XHe*P;

cd h2o
f=0.5:0.01:30;
for k=1:length(f)
    alpha_watervapor_deboer_corrected(k)=fwatervaporalpha_deboer_corrected(f(k),T,P_H2,P_He,P_H2O);
    alpha_watervapor_deboer_original(k)=fwatervaporalpha_deboer_original(f(k),T,P_H2,P_He,P_H2O);
    alpha_watervapor_goodman(k)=fwatervaporalpha_goodman(f(k),P,T,XH2O,XH2,XHe);
end
figure(1)
semilogy(f,alpha_watervapor_deboer_corrected*OpticaldepthstodB)
hold on;
semilogy(f,alpha_watervapor_deboer_original*OpticaldepthstodB,'k')
semilogy(f,alpha_watervapor_goodman*OpticaldepthstodB,'r')
hold off;
