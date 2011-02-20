%Inputs
clear;
dz=5;
XHe_i=0.11;
XH2S_i=2.9e-5;
XNH3_i=2.5e-4;
XH2O_i=1.5e-4;
XCH4_i=6.0e-4;
XPH3_i=0;
XCO=0.0;
funtimes={'Pressure' 'Temperature' 'z' 'XH_2' 'XHe' 'XH_2S' 'XNH_3' 'XH_2O' 'XCH_4' 'XPH_3' 'clouds' 'DNH_4SH' 'DH_2S' 'DNH_3' 'DH_2O' 'DCH_4' 'DPH_3' 'DSOL' 'g' 'mu' 'refr_w/o' 'refr_ w/'};
P_temp=500;
T_temp=1081;
g0_i=2417;
R0e_i=1.0e12;
P0_i=1;
T_targ_i=165;
P_targ_i=1.0;
P_term_i=0.141;
use_lindal=1;
n_lindal_pts=21;

SuperSatSelf_H2S=0;
SuperSatSelf_NH3=0;
SuperSatSelf_PH3=0;
SuperSatSelf_H2O=0;
supersatNH3=0;
supersatH2S=0;
AutoStep_constant=0;

cd ..
adams_data_directory='adams_model_3x';
adams_suffix='j08';
[adams,tcmp]=get_adams_data(1,adams_data_directory,adams_suffix);
cd TCM_mex
[a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v]=TCM(dz,XHe_i,XH2S_i,XNH3_i,XH2O_i,XCH4_i,XPH3_i,XCO,P_temp,T_temp,g0_i,R0e_i,P0_i,T_targ_i,P_targ_i,P_term_i,use_lindal,n_lindal_pts,SuperSatSelf_H2S,SuperSatSelf_NH3,SuperSatSelf_PH3,SuperSatSelf_H2O,supersatNH3,supersatH2S,AutoStep_constant);

%previous_version=load('../TCM_saturn_test/modelfulldir/TCM.out');
current_version=[a,b,v,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u];

[m,n]=size(current_version)
% for i=1:n
%     for j=1:m
%         fun(j,i)=current_version(j,i)-previous_version(j,i);
%     end
% end


for i=1:21
    figure(i)
    axes('YTicklabel',{'0.1','1','10','100','1000'},...
      'Ytick',[0.1,1,10,100,1000],...
      'Yscale','log',...
      'YminorTick','on',...
      'Ydir','reverse','FontSize',20)
  box('on')
  hold('on')

  %  plot(previous_version(1:m,i),'r');
  %  hold on;
    plot(current_version(1:m,i),current_version(1:m,1),'k');
    hold on;
    plot(adams(1:length(adams),i),adams(1:length(adams),1),'b');
    title(funtimes(i));
    %hold off;
end
