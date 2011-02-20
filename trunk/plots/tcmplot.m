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
cnh4sh=[0;(tcme(recordlength-k,12))];
ch2s=[0;(tcme(recordlength-k,13))];
cnh3=[0;(tcme(recordlength-k,14))];
ch2o=[0;(tcme(recordlength-k,15))];
cch4=[0;(tcme(recordlength-k,16))];
cph3=[0;(tcme(recordlength-k,17))];
caqu=[0;(tcme(recordlength-k,18))];


figure
plotyy(


loglog(cnh4sh,P,'k')
hold on
loglog(cnh4sh,P,'b')
hold on
%figure
%loglog(ch2s,P)
%figure
loglog(cnh3,P,'.-y')
%figure
loglog(ch2o,P,'r.')
%figure
loglog(caqu,P,'r')



axis ij
axis([1e-9 1e-2 0.01 100])