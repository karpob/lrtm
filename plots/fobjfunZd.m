function objfunZd=fobjfunZd(zz_notunit,Xd_notunit,Yd_notunit,Xo,Yo,Zo,ao,bo,co)

% find zero of this to find range of looks in x,y,z
a=1/(ao^2);
b=1/(bo^2);
c=1/(co^2);
Sr=1;			% for ellipse ==1

Rd=[Xd_notunit Yd_notunit zz_notunit];
Raydirection_length=vectorlength(Rd);

Xd=Xd_notunit/Raydirection_length;
Yd=Yd_notunit/Raydirection_length;
zz=zz_notunit/Raydirection_length;


A=a*(Xd^2)+b*(Yd^2)+c*(zz^2);					% Should equal '1' already (good error check)
B=2*(a*Xd*Xo+b*Yd*Yo+c*zz*Zo);
C=a*(Xo^2)+b*(Yo^2)+c*(Zo^2)-Sr^2;

squarerootterm=B^2-4*A*C;
objfunZd=abs(squarerootterm);