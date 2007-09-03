load tcm.out
recordlength=size(tcm,1);
k=(0:recordlength-1);

% Extract Prameters, flip around
% Add a zero for space
drvector=[(tcm(recordlength-k,3));0];

