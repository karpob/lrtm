function [Tatm,wlayers]=ftam(T,tau,tau_a,masterindex)

% function Tatm=ftam(T,tau,tau_a)
% T is the vector from ?top to bottom of thermal temperatures
% tau is the vector of integrated absorptions from top to bottom ?
% tau_a is the vector of absorptions for each layer top to bottom
% Units are kelvin and opdepths
% The first layer is the ray from the craft to the surface
% its tau, tau_a, T, are zero
% So T(1) should always be zero

loss=exp(-tau);
emit=T(masterindex).*(1-exp(-tau_a));
wlayers=loss.*(1-exp(-tau_a));
layers=loss.*emit;

Tatm=sum(layers);




