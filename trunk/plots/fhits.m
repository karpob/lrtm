function  [hit_x_pos,hit_x_neg,hit_y_pos,hit_y_neg,hit_z_pos,hit_z_neg]=fhits(Rayorigin,ao,bo,co,theta_res)

% From knowledge of spacecraft position-find field of view FOV
theta_res_radians=theta_res*pi/180;
dx=atan(theta_res_radians)*6*ao;		% 1/2 size of planet ;rough -depends on looking
% dx is how big in cartesian-the spot beam is

% number of samples per spot beam (estimate)
N=2;
resolution=dx/N;

% Find nadir

m=vectorlength(Rayorigin)
xo=-Rayorigin(1)/m;
yo=-Rayorigin(2)/m;
zo=-Rayorigin(3)/m;

nadir=[xo yo zo];

% Where can I see
ellipse.a=ao;
ellipse.b=bo;
ellipse.c=co;

% x towards pos direction
limbflag=0;
k=-1;								% want to incremement k at begin-so k=-1, gives k=0
while limbflag==0
   k=k+1							% index increment only if 'while' satisfied
   xd=xo+resolution*k;
   Raydirection=[xd yo zo];
   [intercept,internormal,d,limbflag]=rayellipseint(Rayorigin,Raydirection,ellipse);
   hit_x_pos(k+1,:)=intercept;
end

% x towards neg direction
limbflag=0;
k=0;								% prev loop does k=0+1 already
while limbflag==0
   k=k+1							% index increment only if 'while' satisfied
   xd=xo-resolution*k;
   Raydirection=[xd yo zo];
   [intercept,internormal,d,limbflag]=rayellipseint(Rayorigin,Raydirection,ellipse);
   hit_x_neg(k+1)=intercept;
end

% y towards pos direction
limbflag=0;
k=-1;								% want to incremement k at begin-so k=-1, gives k=0
while limbflag==0
   k=k+1							% index increment only if 'while' satisfied
   yd=yo+resolution*k;
   Raydirection=[xo yd zo];
   [intercept,internormal,d,limbflag]=rayellipseint(Rayorigin,Raydirection,ellipse);
   hit_y_pos(k+1)=intercept;
end

% y towards neg direction
limbflag=0;
k=0;								% prev loop does k=0+1 already
while limbflag==0
   k=k+1							% index increment only if 'while' satisfied
   yd=yo-resolution*k;
   Raydirection=[xo yd zo];
   [intercept,internormal,d,limbflag]=rayellipseint(Rayorigin,Raydirection,ellipse);
   hit_y_neg(k+1)=intercept;
end


% z towards pos direction
limbflag=0;
k=-1;								% want to incremement k at begin-so k=-1, gives k=0
while limbflag==0
   k=k+1							% index increment only if 'while' satisfied
   zd=zo+resolution*k;
   Raydirection=[xo yo zd];
   [intercept,internormal,d,limbflag]=rayellipseint(Rayorigin,Raydirection,ellipse);
   hit_z_pos(k+1)=intercept;
end

% z towards neg direction
limbflag=0;
k=0;								% prev loop does k=0+1 already
while limbflag==0
   k=k+1							% index increment only if 'while' satisfied
   zd=zo-resolution*k;
   Raydirection=[xo yo zd];
   [intercept,internormal,d,limbflag]=rayellipseint(Rayorigin,Raydirection,ellipse);
   hit_z_neg(k+1)=intercept;
end

