function [tcme,tcmp,Re,gravity]= get_orton_data(oblateness_factor,Planet)

%oblateness factor==ratio of polar to equatorial radius.

%Planet=String with name of the desired planet.

if(strcmp(Planet,'Jupiter'))
	name='Orton_compare/Jupiter_actual.txt'
        Xh2=0.864;
        Xhe=0.134;
        gravity=2400.0; %(cm/sec^2)
        Re=70000.0;  %km
        Re=Re*1e5; %convert to cm        
elseif(strcmp(Planet,'Saturn'))
        name='Orton_compare/Saturn_actual.txt'
        Xh2=0.881;
        Xhe=0.119;
        gravity=990.0; % cm/sec**2
        Re=60268.0; %km
        Re=Re*1e5; % convert to cm

elseif(strcmp(Planet,'Uranus'))
	name='Orton_compare/Uranus_actual.txt'
        Xh2=0.827;
        Xhe=0.15;
        gravity=913.0;% cm/sec**2
        Re=26228.0;
        Re=Re*1e5; %convert to cm
end
     
Alt_T_P_profile=load(name);
[me,l]=size(Alt_T_P_profile);
XH2=ones(me,1)*Xh2;
XHe=ones(me,1)*Xhe;
level=Alt_T_P_profile(:,1);
Z=Alt_T_P_profile(:,2);
P=Alt_T_P_profile(:,3);
T=Alt_T_P_profile(:,4);

tcme(1:me,1:22)=[P(1:me),T(1:me),Z(1:me),XH2(1:me),XHe(1:me),zeros(me,1),zeros(me,1),zeros(me,1),zeros(me,1),zeros(me,1),...
                 zeros(me,1),zeros(me,1),zeros(me,1),zeros(me,1),zeros(me,1),zeros(me,1),zeros(me,1),zeros(me,1),...
                 zeros(me,1),zeros(me,1),zeros(me,1),zeros(me,1)];

tcmp(1:me,1:22)=[tcme(1:me,1:2),oblateness_factor.*tcme(1:me,3),tcme(1:me,4:22)];
