function [beama,beamb,beamc]=generatebeam(Raydirection,Rayorigin)


% Get beam spread pattern
[b1,b2,b3]=beamsample;

% Rotate each beam to parallel the Raydirection


beama=rotatebeam(Raydirection,b1);
beamb=rotatebeam(Raydirection,b2);
beamc=rotatebeam(Raydirection,b3);





