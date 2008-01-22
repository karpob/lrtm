function [Beamz,Beam_weightz,beam_weighted_ave]=beamsample(Nphi_rings,N_ring_one,BWHM)

% function beamsample
% 
%   This function samples a gaussian beam using Nphi_rings in the phi direction and equally spaced rays in theta. The number of
%   samples in the first ray is supplied, and the other 2 rings are equally spaced according to Nrings_k=N_ring_one*(2*k-1)
%   where k is either the second or third ring. This is called by maintamone.m. Note that this is a major revision from 
%   original Hoffman code where Nphi_rings was hardcoded with a value of 3 (ie. He Ran findraypath for three rings in phi, and averaged
%   them in maintamone.)
%
%
%    VARIABLE DEFINITIONS:
%
%               ->  INPUT:
%                       ->Nphi_rings: Number of rings in the "phi" direction
%                       ->N_ring_one: number of samples in the first ring.
%                       ->BWHM: Beamwidth half-maximum or the 3dB beamwidth of the gaussian antenna beam.
%
%               <-  OUTPUT: 
% 
%                       <-Beamz: Beam samples 2D array with dimensions
%                                X,Y,Z, ray number (kk)
%                       <-Beam_weightz: Beam weights sampled according to a gaussian beam shape
%                                       1D array ordered according to ray number
%                       <-beam_weighted_ave: weighted average of gaussian beam used for nomalization with
%                                            antenna beam and brightness temperature of boresight ray                 


dphi_degree=BWHM/Nphi_rings;
phi_degree=cumsum(dphi_degree.*ones(Nphi_rings,1))
n=phi_degree./dphi_degree % radius multiple from first ring
Ntheta=N_ring_one.*(2*(n)-1); % number of samples for each ring

% center phi on z-axis (around pi/2)

phi_degree_norm=90+phi_degree;
phi=phi_degree_norm*pi/180; % convert to radian
r=1;
delta_beam=dphi_degree/2; % since phi is sampled around by 2pi the delta_phi (from axis is 1/2) AKA halfspace
Beam_X=[];
Beam_Y=[];
Beam_Z=[];
Theta=zeros(Nphi_rings,max(Ntheta));
beam_weighted_ave=0;

for i_ring=1:Nphi_rings
     last_theta_in_ring=2*pi-2*pi/Ntheta(i_ring);
     t=linspace(0,last_theta_in_ring,Ntheta(i_ring));
     Theta(i_ring,1:length(t))=t';
     [R(i_ring,1:length(t)),TH(i_ring,1:length(t)),PHI(i_ring,1:length(t))]=meshgrid(r,Theta(i_ring,1:length(t)),phi(i_ring));
     [x(i_ring,1:length(t)),y(i_ring,1:length(t)),z(i_ring,1:length(t))]=sph2cart(TH(i_ring,1:length(t)),PHI(i_ring,1:length(t)),R(i_ring,1:length(t)));
     Beam_X(i_ring,1:length(t))=x(i_ring,1:length(t));
     Beam_Y(i_ring,1:length(t))=y(i_ring,1:length(t));
     Beam_Z(i_ring,1:length(t))=z(i_ring,1:length(t));
     delta(i_ring)=delta_beam*i_ring;
     Beam_weight(i_ring)=(1/(N_ring_one*(2*i_ring-1)))*exp(-2.76*(delta(i_ring)/BWHM)^2);
     beam_weighted_ave=((N_ring_one*(2*i_ring-1)))*Beam_weight(i_ring)+beam_weighted_ave;
     theta_length(i_ring)=length(t);
end
kk=1;
for i_ring=1:Nphi_rings
    for j_theta=1:theta_length(i_ring)
        Beamz(:,kk)=[Beam_X(i_ring,j_theta);Beam_Y(i_ring,j_theta);Beam_Z(i_ring,j_theta)];
        Beam_weightz(kk)=Beam_weight(i_ring);
        kk=kk+1;
    end    
end
