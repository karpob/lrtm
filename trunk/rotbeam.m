function  [Vr1,Zr]=rotbeam(Raydirection,b1) 
x=[1 0 0];
z=[0 0 1];

%c=[-.27639 .85065 -.44721];  % Align z-axis with this unit vector

%vl=vectorlength(Raydirection);
c=Raydirection;

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
%
%
% Example of what I think you're trying to do

%figure
%x=[1 0 0];
%y=[0 1 0];
%z=[0 0 1];
%x1=x;
%y1=y;
%z1=z;
%pp=[-x1(1) 0 0;x1(1) 0 0];
%plot3(pp(:,1),pp(:,2),pp(:,3))
%hold on
%pp=[0 -y1(2) 0;0 y1(2) 0];
%plot3(pp(:,1),pp(:,2),pp(:,3))
%pp=[0 0 -z1(3);0 0 z1(3)];
%plot3(pp(:,1),pp(:,2),pp(:,3))
%xlabel('x')
%ylabel('y')
%zlabel('z')
%grid on

%th=[-180:10:180]*pi/180;
%r=.2;
%pp=[];
%for j=1:length(th)
%   x1=r*sin(th(j));
%   y1=r*cos(th(j));
%   z1=r;
%   p=[x1 y1 z1];
%   pp=[0 0 0;p];
%   plot3(pp(:,1),pp(:,2),pp(:,3),'color','r')
%end
%pp=[];
%pp=[0 0 0;z(1) z(2) z(3)]
%plot3(pp(:,1),pp(:,2),pp(:,3),'color','r')

%pp=[];
%   x1=r*sin(th(j));
%   y1=r*cos(th(j));
%   z1=r;
%   p=[x1 y1 z1];
Zr=(Rx*Ry*Rz*inv(Ry)*inv(Rx)*z');
Vr1=(Rx*Ry*Rz*inv(Ry)*inv(Rx)*b1);
%Vr2=(Rx*Ry*Rz*inv(Ry)*inv(Rx)*b2);
%Vr3=(Rx*Ry*Rz*inv(Ry)*inv(Rx)*b3);

%   p=(Rx*Ry*Rz*inv(Ry)*inv(Rx)*p');
%   pp=[0 0 0;p'];
%   plot3(pp(:,1),pp(:,2),pp(:,3),'color','g')

