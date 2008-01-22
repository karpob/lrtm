function [intercept,internormal,d,t,masterindex,missflag]=findraypath(recordlength,refindex,P,ellipses,Rayorigin,Raydirection)

% function findraypath
%        Flow: maintamone <-findraypath
% 
%        This function calculates the path (ds', neat way of avoiding the use of mu) and direction for all rays.
%
%        Variable Definitions:
%                 
%                     ->INPUT:
%                             ->recordlength: number of layers
%                             ->refindex: refractive index of the layer
%                             ->P: layer(s) pressure(s) in bars
%                             ->ellipses: dimensions for elliptical shells
%                                  ellipses.a: dimensions for elliptical shells along X dimension
%                                  ellipses.b: dimensions for elliptical shells along Y dimension
%                                  ellipses.c: dimensions for elliptical shells along Z dimension
%                             ->Rayorigin: X Y Z position of the spacecraft
%                             ->Raydirection: X Y Z orientation of the spacecraft  
%
%                     <-OUTPUT:
%                              <-intercept: intercept values of each ray with each elliptical shell
%                              <-internormal: normal direction values for each ray/ellipse intersection
%                              <-d: length of each ray (ds) for each layer
%                              <-t: angle relative to normal for each bent ray (theta_2)
%                              <-masterindex: index which maps ray direction arrays?
%                              <-missflag: flag for when the ray misses the atmosphere



%JPH function [intercpt,internormal,d,t]=findraypath(recordlength,refindex,P,ellipses)
%JPH d is the raypath distance, t is the theta2,
%JPH ellipses is an object with ellipses.a/b/c for sizes of ellisphiod axes
%JPH Pmin/max are the min and max pressure stop values -tells when to terminate

global CRITICALFLAG

% Initiate Matrices
intercept=zeros((recordlength+1),3);
internormal=zeros((recordlength+1),3);
rd=zeros((recordlength),3);

% Does spacecraft see the planet intitially???

ellipse.a=ellipses.a(1);
ellipse.b=ellipses.b(1);
ellipse.c=ellipses.c(1);
CRITICALFLAG=0;
limbflag=0;
missflag=0;
%Raydirection
[A,B,C,limbflag]=rayellipseint(Rayorigin,Raydirection,ellipse);
if limbflag==1
   disp('Ray segment misses planet')
   missflag=1;
   intercept=[0 0 0];
   internormal=[0 0 0];
   d=0;
   t=0;
   masterindex=1;
   return
end

% I guess it does if limbflag not set to 1
% First layer-space layer  T=2.7K, P=0 etc...

% Have checked to make sure planet in view, it is, so now start
% Without the if statement this does not do limb sounding
% The if kicks out to do outward going part of the limb sounding



% HEY-NOW NEED TO PRECALC RADIUS OF PLANET SHELLS

masterindex=[];
for k=1:length(P)
   P_ray=P(k);
   ellipse.a=ellipses.a(k);
   ellipse.b=ellipses.b(k);
   ellipse.c=ellipses.c(k);

   [A,B,C,limbflag]=rayellipseint(Rayorigin,Raydirection,ellipse);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  If we get a limb flag meaning the ray is tangent to the ellipsoid otherwise skip this bit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if limbflag==1
      sprintf('Limb sounding? Ray has passed out of atmosphere at k= %.0f\n',k)
      m=k-1;									% decrement k back to pre-limb case
      kk=k;						% Value of k from outer loop
      while P_ray >= min(P)   
         % All stays same but sphere-to be tested is now larger one (from previous k)    
         ellipse.a=ellipses.a(m);
         ellipse.b=ellipses.b(m);
         ellipse.c=ellipses.c(m);
         [A,B,C,limbflag]=rayellipseint(Rayorigin,Raydirection,ellipse);
         intercept(kk,:)=A;					% k has still be incremented-keep as index?
         internormal(kk,:)=B;
         d(kk)=C;
         eta1=refindex(m);
         eta2=refindex(m-1);
         % Using regular snells  -taken care of sphere with ray-sphere intersection test
         [theta2,transmitted]=snells(eta1,eta2,internormal(k,:),Raydirection);
         if CRITICALFLAG==1
            d=[d,inf]';
            m=m-1;
            masterindex=[masterindex;m;m-1];
            disp('hit critical')
            return %this should kick back to maintam
         end
         
         t(kk)=theta2;		% saves value of theta2-mostly for debugging
         % New Rayorigin is current intercept
         Rayorigin=intercept(kk,:);
         % New Raydirection is 'transmitted'
         Raydirection=transmitted;
                  
         rd(kk,:)=Raydirection;		% saving ray directions
        
         Psave(kk)=P_ray; 
           
         m=m-1;    %step backward in m (or in pressure)
         masterindex=[masterindex;m]; 
         P_ray=P(m);			% Pressure profile-out of planet
      end
      % Done with limb sounding-strip and return
      d=d'; 
      return
   end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Back to normal (not 90 tangent to ellipsoid)
   intercept(k,:)=A;
   internormal(k,:)=B;
     
   d(k)=C;% Found intercept now find new direction
   eta1=refindex(k);
   if(k<length(P))
       eta2=refindex(k+1);
   else
       eta2=refindex(k);refindex(k);%refindex(k);
   end
   % Using regular snells  -taken care of sphere with ray-sphere intersection test
   [theta2,transmitted]=snells(eta1,eta2,internormal(k,:),Raydirection);
   t(k)=theta2;		% saves value of theta2-mostly for debugging
   % New Rayorigin is current intercept
   Rayorigin=intercept(k,:);
   % New Raydirection is 'transmitted'
   Raydirection=transmitted;
   rd(k,:)=Raydirection;		% saving ray directions
   Psave(k)=P_ray;				% Debug-saves pressure profile
   masterindex=[masterindex;k];
end

d=d';
% comes out with shallowest first (top->bottom)
