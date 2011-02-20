function objfunXd=fobjfunXd(xx_notunit,Yd_notunit,Zd_notunit,Xo,Yo,Zo,ao,bo,co)


% find zero of this to find range of looks in x,y,z
a=1/(ao^2);
b=1/(bo^2);
c=1/(co^2);
Sr=1;			% for ellipse ==1

Rd=[xx_notunit Yd_notunit Zd_notunit];
Raydirection_length=vectorlength(Rd);

xx=xx_notunit/Raydirection_length;
Yd=Yd_notunit/Raydirection_length;
Zd=Zd_notunit/Raydirection_length;


A=a*(xx^2)+b*(Yd^2)+c*(Zd^2);					% Should equal '1' already (good error check)
B=2*(a*xx*Xo+b*Yd*Yo+c*Zd*Zo);
C=a*(Xo^2)+b*(Yo^2)+c*(Zo^2)-Sr^2;

xx_notunit=xx
squarerootterm=B^2-4*A*C;
objfunXd=abs(squarerootterm);