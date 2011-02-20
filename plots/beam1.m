r=[0,1];
th=linspace(0,2*pi,15);
phi=linspace(pi/4,(3/4)*pi,15);
[R,TH,PHI]=meshgrid(r,th,phi);
[x,y,z]=sph2cart(TH,PHI,R);

x=x(:);
y=y(:);
z=z(:);

a=zeros(size(x));

quiver3(a,a,a,x,y,z,0,'.')
