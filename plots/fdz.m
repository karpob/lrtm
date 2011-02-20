function dz=fdz(P)

% Gets dz from daves tcm

dz = (1e3)*log10(P + 8.0);
% The 1e3 is for km

