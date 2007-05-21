function [tau_a,tau]=ftau(kappas,s,masterindex)

% function tau=ftau(kappa,s)
% kappa should be total absorption for each layer (a priori)
% and 's' is the profile of pathlengths in each layer
% kappa (opdepth/cm), s should be in cm
% the tau is t(b,c)=sum[(a=b to c) tau(a)
% where tau(a(0:i))= integral[(0 to i)kappa(s)*ds]
% So find tau(a)'s which is the tau for each layer
% then find the tau(b,c)'s which  is the integrated from b to c
% Must be column vector

% tau(a)'s
% 

tau_a=kappas(masterindex).*s;		% now units of opdepth

% tau(b,c)'s
% I think this works cumsum?  For column vector use cumsum(X,1)

tau1=cumsum(tau_a,1);

st=size(tau1,1);
tau=[0;tau1(1:st-1)];
