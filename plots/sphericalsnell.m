function [theta2,transmitted]=sphericalsnell(eta1,eta2,R1,R2,internormal,Raydirection)

% function [theta2,transmitted]=sphericalsnell(eta1,eta2,internormal,Raydirection)
% Finds theta2, which is the angle of the refracted ray, and the direction
% of that refracted ray, which is termed transmitted.
% The input parameters are: eta1 and eta2 which are the indices of refaction
% for the incident layer (1) and the transmitted layer (2),
% theta1, the incident angle,is the angle between the internormal (which is
% the normal to the surface) and the Raydirection.
% R1 is the radius of the first layer, R2 the second (R1-dz)
% This is done using Heckberts method (Ray Tracing, Glassner ed.)
% With the alteration of snell's law to spherical
% Spherical Snell's:  eta1R1sin(theta1)=eta2R2sin(theta2)


% The incidence angle is found from 

% To find the angle the ray makes with the normal use the dot product...
% theta1=acos(dot(a,b)/(vectorlength(a)*vectorlength(b)))
% where 'a' is the normal (internormal) and 'b' is the incident ray
a=internormal;
b=Raydirection;
theta1=acos(dot(a,b)/(vectorlength(a)*vectorlength(b)))

% To put things in Heckbert notation
I=Raydirection;
N=internormal;

% In Heckbert etarelative=sin(theta2)/sin(theta1)=(eta1/eta2)
% Here, with Spherical Snell's eta relative=(sin(theta2))/(sin(theta1))=(R1*eta1)/(R2*eta2)
etas=(R1*eta1)/(R2*eta2);
c1=cos(theta1);

% c2=cos(theta2) may be found from known quantities =sqrt(1-(etas^2)*(1-c1^2))
c2=sqrt(1-(etas^2)*(1-c1^2));

theta2=acos(c2);

% Trasmitted is found
transmitted=etas.*(I)+(etas*c1-c2).*N;


