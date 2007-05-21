%Inputs
clear;
dz=5;
XHe_i=0.07;
XH2S_i=410e-6;
XNH3_i=520.0e-6;
XH2O_i=1.0e-6;
XCH4_i=4.2e-3;
XPH3_i=0;
XCO=0.0;
funtimes={'Pressure' 'Temperature' 'XH_2' 'XHe' 'XH_2S' 'XNH_3' 'XH_2O' 'XCH_4' 'XPH_3' 'clouds' 'DNH_4SH' 'DH_2S' 'DNH_3' 'DH_2O' 'DCH_4' 'DPH_3' 'DSOL' 'g' 'mu' 'refr_w/o' 'refr_ w/'};
P_temp=70;
T_temp=489.5;
g0_i=900.0;
R0_i=6.03e9;
P0_i=0.5;
T_targ_i=146.2;
P_targ_i=1.29848;
P_term_i=0.06;
use_lindal='Y';
n_lindal_pts=27;
[a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v]=TCM(dz,XHe_i,XH2S_i,XNH3_i,XH2O_i,XCH4_i,XPH3_i,XCO,P_temp,T_temp,g0_i,R0_i,P0_i,T_targ_i,P_targ_i,P_term_i,use_lindal,n_lindal_pts,0,0,0,0);

previous_version=load('../TCM_saturn_test/modelfulldir/TCM.out');
current_version=[a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v];

[m,n]=size(previous_version)
for i=1:n
    for j=1:m
        fun(j,i)=current_version(j,i)-previous_version(j,i);
    end
end

for i=1:21
    figure(i)
    plot(previous_version(1:m,i),'r');
    hold on;
    plot(current_version(1:m,i),'g');
    hold on;
    title(funtimes(i));
    hold off;
end
figure(22)
plot(current_version(1:m,1),current_version(1:m,22));