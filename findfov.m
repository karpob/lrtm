function [x_min,x_max,y_min,y_max,z_min,z_max]=findfov(ao,bo,co)

% Find max extent of x fov given current X0,Y0,Z0 and fixed Yd,Zd
load craft
% Setup
X0=Rayorigin(1);
Y0=Rayorigin(2);
Z0=Rayorigin(3);

Xd_notunit=Raydirection(1);		% Not normalized, need to
Yd_notunit=Raydirection(2);
Zd_notunit=Raydirection(3);


Xc=0;						% ASSUMING Centered
Yc=0;
Zc=0;

Sr=1;						% For ellipse always equal to unity

Raydirection_length=vectorlength(Raydirection);

Xd=Xd_notunit/Raydirection_length;
Yd=Yd_notunit/Raydirection_length;
Zd=Zd_notunit/Raydirection_length;
% Have the numbers, now find the ranges of view

options(1)=0;
options(2)=1e-12;
options(3)=1e-12;

maxextent=max([ao bo co])
xx=fmin('fobjfunXd',0,maxextent,options,Yd,Zd,X0,Y0,Z0,ao,bo,co);
x_max=xx;
xx=fmin('fobjfunXd',-maxextent,0,options,Yd,Zd,X0,Y0,Z0,ao,bo,co);
x_min=xx;

yy=fmin('fobjfunYd',0,maxextent,options,Xd,Zd,X0,Y0,Z0,ao,bo,co);
y_max=yy;
yy=fmin('fobjfunYd',-maxextent,0,options,Xd,Zd,X0,Y0,Z0,ao,bo,co);
y_min=yy;

zz=fmin('fobjfunZd',0,maxextent,options,Xd,Yd,X0,Y0,Z0,ao,bo,co);
z_max=zz;
zz=fmin('fobjfunZd',-maxextent,0,options,Xd,Yd,X0,Y0,Z0,ao,bo,co);
z_min=zz;


