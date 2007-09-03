% Need the mixing ratio profiles for 0.83<P<1.26

load SXkappa
load SXkappa5
load SXkappa10


XLindal=[0.00217246
0.00449081
6.59E-03
9.28E-03
1.12E-02];

XPLindal=[0.83176
0.87096
0.91201
0.95499
1];
plot(XLindal,XPLindal,'k')
Xerr=0.001.*ones(size(XLindal));
hold on
xX=[XLindal Xerr];
yX=[XPLindal zeros(size(XPLindal))];
eplot(xX,yX,'k')
hold on


SLindal=[1.41E-03
2.95E-03
3.66E-03
5.25E-03];

SPLindal=[1.09648
1.14815
1.20226
1.25893];

plot(SLindal,SPLindal,'k')
hold on
Serr=0.001.*ones(size(SLindal));
xS=[SLindal Serr];
yS=[SPLindal zeros(size(SPLindal))];
eplot(xS,yS,'k')
hold on
axis ij
hold on


plot(Skappa,P,'k')
hold on
plot(Xkappa,P,'k')

plot(Skappa5,P5,':k')
hold on
plot(Xkappa5,P5,':k')

plot(Skappa10,P10,'--k')
hold on
plot(Xkappa10,P10,'--k')

axis([0 0.015 0.65 1.4])

xlabel('\fontsize{16}Absorption (dB/km)')
ylabel('\fontsize{16}Pressure (bar)')