% Need the mixing ratio profiles for 0.83<P<1.26
donew=1
if donew==1
pres=cd
cd(pres)
load tcme.out
%load c:\models\neptunedave\tcm.out
%tcme=tcm;
%cd c:\models\rtm\imaging
recordlength=size(tcme,1);
k=(0:recordlength-1);
%
%P	T	dr XH2	XHe	XH2S	XNH3	XH2O	XCH4	XPH3	clouds	DNH4SH	DH2S	DNH3	DH2O	DCH4	DPH3	DSOL	g	mu	refr_w/o	refr_w/


% Extract Prameters, flip around
% Add a zero for space
P=[0;(tcme(recordlength-k,1))];
T=[2.7;(tcme(recordlength-k,2))];
% third is the dR vector -handled elsewhere since need two
major=1e5.*[(tcme(recordlength-k,3))];
%minor=1e5.*[(tcmp(recordlength-k,3))];
xH2=[0;(tcme(recordlength-k,4))];
xHe=[0;(tcme(recordlength-k,5))];
xH2S=[0;(tcme(recordlength-k,6))];
xNH3=[0;(tcme(recordlength-k,7))];
xH2O=[0;(tcme(recordlength-k,8))];
xCH4=[0;(tcme(recordlength-k,9))];
xPH3=[0;(tcme(recordlength-k,10))];
% cause ph3 decay (photolyse)
%decay=ph3decay(P,xPH3);
decay=ph3decay(P,xPH3);
xPH3=decay;

% Convert to partial pressures
P_H2=P.*xH2;
P_He=P.*xHe;
P_H2S=P.*xH2S;
P_NH3=P.*xNH3;
P_H2O=P.*xH2O;
P_CH4=P.*xCH4;
P_PH3=P.*xPH3;

f=2.3;
Skappa=findkappanosave(f,T,P_H2,P_He,P_NH3,P_H2O,P_CH4,P_PH3,P_H2S);
f=8.4;
Xkappa=findkappanosave(f,T,P_H2,P_He,P_NH3,P_H2O,P_CH4,P_PH3,P_H2S);

% to dB/km
todb=434294.5;
Skappa=todb.*Skappa;
Xkappa=todb.*Xkappa;
else
end






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

%plot(XLindal,XPLindal,'b')
%Xerrhigh=XLindal+0.001;
%plot(Xerrhigh,XPLindal,'bh')
%Xerrlow=XLindal-0.001;
%plot(Xerrlow,XPLindal,'bo')
plot(XLindal,XPLindal,'b')
Xerr=0.001.*ones(size(XLindal));
%errorbar(XLindal,XPLindal,Xerr,'y')
hold on
xX=[XLindal Xerr];
yX=[XPLindal zeros(size(XPLindal))];
eplot(xX,yX,'b')
hold on

SLindal=[1.41E-03
2.95E-03
3.66E-03
5.25E-03];

SPLindal=[1.09648
1.14815
1.20226
1.25893];

plot(SLindal,SPLindal,'r')
hold on
Serr=0.001.*ones(size(SLindal));


xS=[SLindal Serr];
yS=[SPLindal zeros(size(SPLindal))];
eplot(xS,yS,'r')
hold on
axis ij
hold on


plot(Skappa,P,'m')
hold on
plot(Xkappa,P,'g')

