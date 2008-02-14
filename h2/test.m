clear all;
OpticaldepthstodB=434294.5;
f=0.1:1:500;
P=100; %pressure 1 bar
T=293;
PH2=P;
PHe=0.13*P;
PCH4=0.02*P;
for i=1:length(f)
    alphah2_g_by_j(i)=falphah2_goodman_by_joiner(f(i),T,PH2,PHe);
    alphah2_g(i)=falphah2_goodman(f(i),T,PH2,PHe);
    alphah2_j(i)=falphah2(f(i),T,PH2,PHe,PCH4);
    alphah2_b(i)=falpha_borysow(f(i),T,P);
    alphah2_o(i)=falpha_orton(f(i),T,P);
    alphah2_o_mod(i)=falpha_orton_mod(f(i),T,P);
    alphah2_b_pph2(i)=falpha_borysow(f(i),T,PH2);
    alphah2_o_pph2(i)=falpha_orton(f(i),T,PH2);
end
save('1_500GHz_293K_100bar.mat')
%plot(f,OpticaldepthstodB*alphah2_o,'k')
%hold on;
%plot(f,OpticaldepthstodB*alphah2_o_mod,'r')
%hold off
%save('1_30GHz.mat')
%plot(f,OpticaldepthstodB*spec*amagat^2)
%hold on;
%plot(f,OpticaldepthstodB*spec_mod*amagat^2,'k')


