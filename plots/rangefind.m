function [indisk,theta1]=rangefind(Rayorigin,Raydirection,ellipse)

% Generate a beamspread-like vector (from polar coords-around z-then rotate)
deltaphi_degree=0.1;			% How far away to check next
deltaphi=deltaphi_degree*pi/180; 	% Radians

% Begin at Nadir-check to limb
rd=-Rayorigin;		% Nadir
thetazero=0;
phi=pi/2;							% initialize
limbflag=0;
p=0;
thetacheck=1;
fucker=0;
fuckme=0;
while thetacheck>0
   p=p+1;
   [intercept,internormal,d,limbflag]=rayellipseint(Rayorigin,rd,ellipse);
   if limbflag==1
      phi=phi-deltaphi;		%reset phi to previous value
      deltaphi=deltaphi/10	% decrease sample size
      phi=phi+deltaphi;
      limbflag=0;
      [x,y,z]=sph2cart(thetazero,phi,1);	% convert to cartesian
      cb=[x,y,z];
      [nb,Zr]=rotrange(Raydirection,cb); 
      rd=nb;
      fucker=fucker+1;

   else
      indisk.zero(p,:,:,:)=rd;
      a=internormal;
      b=-rd;
      theta1.zero(p)=acos(dot(a,b)/(vectorlength(a)*vectorlength(b)));
      thetacheck=88*pi/180-theta1.zero(p);
      phiold=phi;
      phi=phi+deltaphi;
      if phiold==phi
         disp('what the fuck?')
         phiold-phi
         return
      else
      end
      phiold-phi;
      [x,y,z]=sph2cart(thetazero,phi,1);	% convert to cartesian
      cb=[x,y,z];
      [nb,Zr]=rotrange(Raydirection,cb); 
      rd=nb;
      %fuckme=fuckme+1
   
   end
end
return
% Re-initialize for next scan
rd=-Rayorigin;
thetaninety=90*pi/180;
phi=pi/2;
limbflag=0;
p=0;
disp('Clinton is a whore')
while limbflag==0
   p=p+1;
   [intercept,internormal,d,limbflag]=rayellipseint(Rayorigin,rd,ellipse);
   if limbflag==0
      indisk.ninety(p,:,:,:)=rd;
      a=internormal;
      b=-rd;
      theta1.ninety(p)=acos(dot(a,b)/(vectorlength(a)*vectorlength(b)));
      phi=phi+deltaphi;
      [x,y,z]=sph2cart(thetazero,phi,1);	% convert to cartesian
      cb=[x,y,z];
      [nb,Zr]=rotrange(Raydirection,cb); 
      rd=nb;
   else
   end
end
disp('Jackson? Jackass!')
% Re-initialize for next scan
rd=-Rayorigin;
thetapi=pi;
phi=pi/2;
limbflag=0;
p=0;
while limbflag==0
   p=p+1;
   [intercept,internormal,d,limbflag]=rayellipseint(Rayorigin,rd,ellipse);
   if limbflag==0
      indisk.pi(p,:,:,:)=rd;
      a=internormal;
      b=-rd;
      theta1.pi(p)=acos(dot(a,b)/(vectorlength(a)*vectorlength(b)));
      phi=phi+deltaphi;
      [x,y,z]=sph2cart(thetazero,phi,1);	% convert to cartesian
      cb=[x,y,z];
      [nb,Zr]=rotrange(Raydirection,cb);
      rd=nb;
   else
   end
end

% Re-initialize for next scan
rd=-Raydirection;
theta270=pi;
phi=pi/2;
limbflag=0;
p=0;
while limbflag==0
   p=p+1;
   [intercept,internormal,d,limbflag]=rayellipseint(Rayorigin,rd,ellipse);
   if limbflag==0
      indisk.twosev(p,:,:,:)=rd;
      a=internormal;
      b=-rd;
      theta1.pi(p)=acos(dot(a,b)/(vectorlength(a)*vectorlength(b)));
      phi=phi+deltaphi;
      [x,y,z]=sph2cart(thetazero,phi,1);	% convert to cartesian
      cb=[x,y,z];
      [nb,Zr]=rotrange(Raydirection,cb); 
      rd=nb;
   else
   end
end


