function [DNH4SH,DH2S,DNH3,DH2O,DCH4,DPH3,DSOL]=filter_clouds(clouds,DNH4SH,DH2S,DNH3,DH2O,DCH4,DPH3,DSOL)
% function filter_clouds 
%
% This function is used by DeBoer_TCM.m to ensure that clouds exist above DeBoer's threshold values defined by his long int 
% "clouds" which states the cloud type, phase at each altitude.
%
%     Variable Definitions:
%               
%      ---->INPUT:         clouds: Long int taken from DeBoer's C program which gives cloud types and phase present at each
%                                  pressure level.
%     <---->INPUT/OUTPUT   DNH4SH: cloud density of the ammonium hydrosulfide cloud (g/cm^3)
%                          DH2S: cloud density of the hydrogen sulfide cloud (g/cm^3)
%                          DNH3: cloud density of the ammonia ice cloud (g/cm^3)
%                          DH2O: cloud density of the water ice cloud (g/cm^3)
%                          DCH4: cloud density of the methane cloud (g/cm^3)
%                          DPH3: cloud density of the phosphine cloud (g/cm^3)
%                          DSOL: cloud density of the ammonia-water solution cloud (g/cm^3)
%
%            
%
%
%Filter values for clouds, make sure that the "clouds" variable is
%consistent with cloud bulk density

SOL_filter='000001'; %makes sure cloud is a liquid cloud 
H2Oice_filter='000002'; % makes sure cloud is a H2O ice cloud
NH4SH_filter='000020'; %makes sure there is a NH4SH ice cloud
NH3_filter='000200'; %makes sure there is a NH3 ice cloud
H2S_filter='002000'; %makes sure there is a H2S ice cloud
CH4_filter='020000'; %makes sure there is a CH4 ice cloud
PH3_filter='200000'; %makes sure there is a PH3 ice cloud
me=length(DNH4SH);
%
%
for i=1:me
    cval=num2str(clouds(i)); %get cloud phase value can convert to string
%   Do kludge to pad the "cloud vector" with zeros    
    if(length(cval)==1)
        cval1=strcat('00000',cval);
    end
    if(length(cval)==2)
        cval1=strcat('0000',cval);
    end
    if(length(cval)==3)
        cval1==strcat('000',cval);
    end
    if(length(cval)==4)
        cval1==strcat('00',cval);
    end
    if(length(cval)==5)
        cval1=strcat('0',cval);
    end
% end padding kludge    

%  Now we check to make sure each cloud really exists, if it doesn't
%  set the cloud density to zero. Hey, its ugly, but it works...I think
    if(cval1(6)~=SOL_filter(6))
        DSOL(i)=0;
    end
    if(cval1(5)~=NH4SH_filter(5))
        DNH4SH(i)=0;
    end

   % if(cval1(4)~=NH3_filter(4))
   %     DNH3(i)=0;
   % end
    
    if(cval1(3)~=H2S_filter(3))
        DH2S(i)=0;
    end
    
    if(cval1(2)~=CH4_filter(2))
        DCH4(i)=0;
    end
    if(cval1(1)~=PH3_filter(1))
        DPH3(i)=0;
    end
end
