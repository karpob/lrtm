function Vr=rotatebeam(A,V)


Xa=A(1);
Ya=A(2);
Za=A(3);

Xv=V(1);
Yv=V(2);
Zv=V(3);

gYaZa=sqrt(Ya^2+Za^2);

m11=gYaZa;
m12=-Xa*Ya/gYaZa;
m13=-Xa*Za/gYaZa;
m21=0;
m22=-Za/gYaZa;
m23=-Ya/gYaZa;
m31=Xa;
m32=Ya;
m33=Za;

Rzppa=[m11 m12 m13;m21 m22 m23;m31 m32 m33];
Vr=Rzppa*V;
