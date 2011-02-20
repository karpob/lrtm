function objfunYd=fobjfunYd(yy_notunit,Xd_notunit,Zd_notunit,Xo,Yo,Zo,ao,bo,co)

% find zero of this to find range of looks in x,y,z
a=1/(ao^2);
b=1/(bo^2);
c=1/(co^2);
Sr=1;			% for ellipse ==1

Rd=[Xd_notunit yy_notunit Zd_notunit];
Raydirection_length=vectorlength(Rd);

Xd=Xd_notunit/Raydirection_length;
yy=yy_notunit/Raydirection_length;
Zd=Zd_notunit/Raydirection_length;



A=a*(Xd^2)+b*(yy^2)+c*(Zd^2);					% Should equal '1' already (good error check)
B=2*(a*Xd*Xo+b*yy*Yo+c*Zd*Zo);
C=a*(Xo^2)+b*(Yo^2)+c*(Zo^2)-Sr^2;

squarerootterm=B^2-4*A*C;
objfunYd=abs(squarerootterm);