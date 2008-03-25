fu=fopen('BEAM3_V01.PAT','r');
[x,y]=fread(fu,'double');
fclose(fu);
data=reshape(x,1200,400);

x=[-2:0.01:1.99];
y=[-6:0.01:5.99];
xx=(pi/180)*x;
yy=(pi/180)*y;
save('beam_pattern.mat')
