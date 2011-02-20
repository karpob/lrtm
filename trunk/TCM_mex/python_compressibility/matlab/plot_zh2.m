%plot Z
clear all
P=[0.1:0.1:100];
T=[293;375;450;500;525];
R=8.314472/100;
for i=1:length(P)
    for j=1:length(T)
    Z(i,j)=Z_H2(P(i),T(j));
    rho(i,j)=P(i)/(R*T(j)*Z(i,j));
    Za(i,j)=Z_H2_atreya(rho(i,j),T(j));
    err(i,j)=100*(Za(i,j)-Z(i,j))/Z(i,j);
    end
end
hold off;
figure(1)
plot(P,Z)
xlabel('Actual Pressure (bars)')
ylabel('Compressibility (Z=\frac{P_{real}}{P_{ideal})')
figure(2)
plot(P,Za)
xlabel('Actual Pressure (bars)')
ylabel('Compressibility (Z=\frac{P_{real}}{P_{ideal})')
figure(3)
plot(P,err)
xlabel('Actual Pressure (bars)')
ylabel('Percent Error (\frac{Z_{Lemmon}-Z_{Atreya}}{Z_{Lemmon}} \times 100)')
% %Test points from Lemmon, Huber, and Leachman

T_test=[200;300;400;500;200];
P_test=[10;100;500;2000;2000];

for i=1:length(P_test)
    Z_test(i)=Z_H2(P_test(i),T_test(i));
end

rho_test=P_test./(R.*T_test.*transpose(Z_test));