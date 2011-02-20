h2=0.92338;
he=0.0563;
nh3=0.02032;
T=186;
f=[32:0.1:41]';
P=1.964;
quickbg
quickspil
[alphanh3,alphaspilker]*434294.5
alpha=[alphanh3]*434294.5
x=[32.7810 35.6773 37.8499 40.0236];
y=[405.3069 260.8589 220.6590 175.0842];
e=[3.6725 1.5115 1.8937 1.1408];
plot(f,alpha);
hold on;
errorbar(x,y,e,'k--');
hold off;
xlabel('Frequency(GHz)');
ylabel('Opacity(dB/km)');
title('P=1.964 bars NH3, T=186 K')