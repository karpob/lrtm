function ellipses=findellipseradiusvector(ao,bo,co,major,minor)

% function ellipses=findellipseradiusvector(ao,bo,co,major,minor)
% from master ellipse (largest extent of the planet)
% finds the ellipsoidal shells for the elliptical-shell model
% ao,bo,co are the largest extent of each axis (x,y,z)
% From the TCM.out (which gives dR for each P relative to Ro(P=0.5)
% After running model (TCM) twice with Ro=major axis, minor axis) gives
% Vector of ellipsoid vertices, from that figures out the data points in between


a=ao+major;				% semi-minor axis is x
b=bo+major;				% semi-minor axis is also y
c=co+minor;				% semi-major axis is z

ellipses.a=[a;min(a)];
ellipses.b=[b;min(b)];
ellipses.c=[c;min(c)];
