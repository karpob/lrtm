function  [nb,Zr]=rotrange(Raydirection,cb) 
x=[1 0 0];
z=[0 0 1];
%c=[-.27639 .85065 -.44721];  % Align z-axis with this unit vector
vl=vectorlength(Raydirection);
c=Raydirection./vl;
rvect_e=cross(z,c);
rvect_e=rvect_e/norm(rvect_e);
rang_e=acos(dot(z,c));
kx=rvect_e(1);
ky=rvect_e(2);
kz=rvect_e(3);
v=norm([ky kz]);
Rx=[1 0 0;
   0 kz/v ky/v;
   0 -ky/v kz/v];
Ry=[v 0 kx;
   0 1 0;
   -kx 0 v];
Rz=[cos(rang_e) -sin(rang_e) 0;
   sin(rang_e) cos(rang_e) 0;
   0 0 1];

% To transform any vector to align with the unit vector c
% Multiply that vector with this matrix
%
%  new_vector=Rx*Ry*Rz*inv(Ry)*inv(Rx)*old_vector
Zr=(Rx*Ry*Rz*inv(Ry)*inv(Rx)*z');
nb=(Rx*Ry*Rz*inv(Ry)*inv(Rx)*cb');

%   p=(Rx*Ry*Rz*inv(Ry)*inv(Rx)*p');
%   pp=[0 0 0;p'];
%   plot3(pp(:,1),pp(:,2),pp(:,3),'color','g')

