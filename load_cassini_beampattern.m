function [Beamz,Beam_weightz,beam_sum]=load_cassini_beampattern(path)

% function load_cassini_beampattern
% 
%   This function load the Cassini beam pattern stored in 'path'
%
%    VARIABLE DEFINITIONS:
%
%               ->  INPUT:
%                       ->path: path to cassini antenna data file
%
%               <-  OUTPUT: 
% 
%                       <-Beamz: Beam samples 2D array with dimensions
%                                X,Y,Z, ray number (kk)
%                       <-Beam_weightz: Beam weights sampled according to a gaussian beam shape
%                                       1D array ordered according to ray number
%                       <-beam_weighted_ave: weighted average of gaussian beam used for nomalization with
%                                            antenna beam and brightness temperature of boresight ray           
%        'BEAM3_V01.PAT'
      
fu=fopen(path,'r');
[data_prime,length_of_data]=fread(fu,'double');
fclose(fu);
data=reshape(data_prime,1200,400);
Beam_weight=data;
theta_degree=[-2:0.01:1.99];
phi_degree=[-6:0.01:5.99];


theta=theta_degree*(pi/180);
phi=phi_degree*(pi/180);


Beam_X=[];
Beam_Y=[];
Beam_Z=[];
beam_weighted_ave=0;

[TH,PHI]=meshgrid(theta,phi);
[mm,nn]=size(TH);
R=ones(mm,nn);
[Beam_X,Beam_Y,Beam_Z]=sph2cart(TH,PHI,R);

% reshape is used to redimension arrays and to replace the following:
%
%kk=1;
%for i_phi=1:length(phi_degree)
%    for j_theta=1:length(theta)
%        Beamz(:,kk)=[Beam_X(j_theta,i_phi);Beam_Y(j_theta,i_phi);Beam_Z(j_theta,i_phi)];
%        Beam_weightz(kk)=Beam_weight(i_phi,j_theta);
%        kk=kk+1;
%    end    
%end



[m,n]=size(Beam_X);
Beamz_X=reshape(Beam_X,1,m*n);
Beamz_Y=reshape(Beam_Y,1,m*n);
Beamz_Z=reshape(Beam_Z,1,m*n);

Beamz=[Beamz_X;Beamz_Y;Beamz_Z];

Beam_weightz=reshape(Beam_weight,1,m*n);
if(Beam_weightz(:)<0)
    Beam_weightz(:)=0;
end

beam_sum=sum(Beam_weightz);
