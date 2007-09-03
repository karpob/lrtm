%look=[0 1 0];
[THlook,PHIlook,Rlook]=cart2sph(0,0,1);
deltatheta=10*pi/180;

r=1;
%th=linspace(THlook-deltatheta,THlook+deltatheta,15);
th=linspace(0,2*pi,15);

phi=linspace(PHIlook-deltatheta,PHIlook+deltatheta,15);
%phi=linspace(-pi,pi,15);

[R,TH,PHI]=meshgrid(r,th,phi);
[x,y,z]=sph2cart(TH,PHI,R);

x=x(:);
y=y(:);
z=z(:);

a=zeros(size(x));

h=quiver3(a,a,a,x,y,z,0,'.')
%rotate(h,[0 1 0],60,[0 0 0])
axis([-1 1 -1 1 0 1])