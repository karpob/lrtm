function [Beamz,Beam_weightz,beam_weighted_ave]=load_cassini_beampattern(path)

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
r=1;

Beam_X=[];
Beam_Y=[];
Beam_Z=[];
beam_weighted_ave=0;
[Beam_X,Beam_Y,Beam_Z]=meshgrid(r,phi,theta)
[R,TH,PHI]=meshgrid(r,theta,phi);
[Beam_X,Beam_Y,Beam_Z]=sph2cart(TH,PHI,R);
kk=1;
for i_phi=1:length(phi_degree)
    for j_theta=1:lenth(theta)
        Beamz(:,kk)=[Beam_X(i_phi,j_theta);Beam_Y(i_phi,j_theta);Beam_Z(i_phi,j_theta)];
        Beam_weightz(kk)=Beam_weight(i_phi,j_theta);
        kk=kk+1;
    end    
end
beam_weighted_ave=1;
