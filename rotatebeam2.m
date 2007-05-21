function [Vrx,Vry,Vrz]=rotatebeam2(Vx,Vy,Vz,alpha)

V=[Vx;Vy;Vz];


Rx11=1;
Rx12=0;
Rx13=0;
Rx21=0;
Rx22=cos(alpha);
Rx23=sin(alpha);
Rx31=0;
Rx32=-sin(alpha);
Rx33=cos(alpha);

Rx=[Rx11 Rx12 Rx13;Rx21 Rx22 Rx23;Rx31 Rx32 Rx33];
Vrx=Rx*V;

Ry11=cos(alpha);
Ry12=0;
Ry13=sin(alpha);
Ry21=0;
Ry22=1;
Ry23=0;
Ry31=-sin(alpha);
Ry32=0;
Ry33=cos(alpha);

Ry=[Ry11 Ry12 Ry13;Ry21 Ry22 Ry23;Ry31 Ry32 Ry33];
Vry=Ry*V;

Rz11=cos(alpha);
Rz12=sin(alpha);
Rz13=0;
Rz21=-sin(alpha);
Rz22=cos(alpha);
Rz23=0;
Rz31=0;
Rz32=0;
Rz33=1;

Rz=[Rz11 Rz12 Rz13;Rz21 Rz22 Rz23;Rz31 Rz32 Rz33];
Vrz=Rz*V;
