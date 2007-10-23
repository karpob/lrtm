clear all;
f=1:0.1:30;
P=100; %pressure 1 bar
T=293;
PH2=0.85*P;
PHe=0.13*P;
PCH4=0.02*P;
for i=1:length(f)
    alphah2_g_by_j(i)=falphah2_goodman_by_joiner(f(i),T,PH2,PHe);
    alphah2_g(i)=falphah2_goodman(f(i),T,PH2,PHe);
    alphah2_j(i)=falphah2(f(i),T,PH2,PHe,PCH4);
    alphah2_b(i)=falpha_borysow(f(i),T,P);
    alphah2_o(i)=falpha_orton(f(i),T,P);
    alphah2_b_pph2(i)=falpha_borysow(f(i),T,PH2);
    alphah2_o_pph2(i)=falpha_orton(f(i),T,PH2);
end
size(alphah2_o)
size(alphah2_j)

save('1_30GHz.mat')
%plot(f,OpticaldepthstodB*spec*amagat^2)
%hold on;
%plot(f,OpticaldepthstodB*spec_mod*amagat^2,'k')


