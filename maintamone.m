function [Tbeam,jims_zenith,wfwa,refindex,...
          intercepts_boresight,intercepts_b,...
          intercepts_c,intercepts_d]= maintamone(Spherecenter,Sphereradius,Raydirection,Rayorigin,...
                                                 tcme,tcmp,ao,bo,co,f,no_ph3,...
                                          select_h2h2_model,select_ammonia_model,select_water_model,include_clouds,N_ring_one,BWHM,refractivity_source)
global CRITICALFLAG
CRITICALFLAG=0				% If is '1' then critical refraction reached for that ray
USEBEAM=1;

% Want to flip order from low atltitude to high over to top-down
recordlength1=size(tcme,1);
recordlength2=size(tcmp,1);
if(recordlength1>recordlength2)
    recordlength=recordlength2;
end
if(recordlength2>recordlength1)
    recordlength=recordlength1;
end
if(recordlength1==recordlength2)
    recordlength=recordlength1;
end
k=(0:recordlength-1);
%
% Extract Prameters, flip around
% Add a zero for space
P=[0;(tcme(recordlength-k,1))];
T=[2.7;(tcme(recordlength-k,2))];
% third is the dR vector -handled elsewhere since need two
major=1e5.*[(tcme(recordlength-k,3))];
minor=1e5.*[(tcmp(recordlength-k,3))];
xH2=[0;(tcme(recordlength-k,4))];
xHe=[0;(tcme(recordlength-k,5))];
xH2S=[0;(tcme(recordlength-k,6))];
xNH3=[0;(tcme(recordlength-k,7))];
xH2O=[0;(tcme(recordlength-k,8))];
xCH4=[0;(tcme(recordlength-k,9))];
xPH3=[0;(tcme(recordlength-k,10))];
% Get cloud Densities from TCME
DNH4SH=[0;(tcme(recordlength-k,12))];
DH2S=[0;(tcme(recordlength-k,13))];
DNH3=[0;(tcme(recordlength-k,14))];
DH2O=[0;(tcme(recordlength-k,15))];
DSOL=[0;(tcme(recordlength-k,18))];

% % cause ph3 decay (photolyse)
  if (no_ph3>0)
      look_out_for_decaying_ph3=1
      xPH3(1)=xPH3(1);
      xPH3(recordlength);
      decay=ph3decay(P,xPH3);
      xPH3=decay;
  end

smallestdr=min(length(major),length(minor));
minor=minor(1:smallestdr);		% unify lengths
major=major(1:smallestdr);    % unify lengths

% Convert to partial pressures
P_H2=P.*xH2;
P_He=P.*xHe;
P_H2S=P.*xH2S;
P_NH3=P.*xNH3;
P_H2O=P.*xH2O;
P_CH4=P.*xCH4;
P_PH3=P.*xPH3;

sph3=size(P_PH3);


% Find the refractive index profile
refindex=findrefindex(T,P_H2,P_He,P_CH4,P_H2O,DSOL,f,refractivity_source);

% Find the elliptical shells used by 'findraypath' (major/minor) are vectors of how radius changes with each index
ellipses=findellipseradiusvector(ao,bo,co,major,minor);

% IMAGING (mini)


% DOES BORESIGHT MISS-IF SO=> END
[intercept,internormal,d,t,masterindexa,missflag]=findraypath(recordlength,refindex,P,ellipses,Rayorigin,Raydirection);

if missflag==1
   disp('Boresight Misses Planet')
   Tbeam=2.7;
   weightingfactor=0;
   return				% This will help with full imaging (if this file becomes a function
   % that gets passed various raydirection(boresights)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Calculate Jim's zenith
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% what Jim Hoffman calls zenith is the angle between normal to the surface,
% and observer. Not what you'd think when he talks about nadir, and
% zenith. Quite confusing if you ask me.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
jims_zenith=acos(dot(-internormal(1,:),Raydirection))*(180/pi)

%Calculate absorption coeff as a function of TP 
kappa=findkappa(f,T,P,P_H2,P_He,P_NH3,P_H2O,P_CH4,P_PH3,P_H2S,...
                xH2,xHe,xNH3,xH2O,DNH4SH,DH2S,DNH3,DH2O,DSOL,...
                select_h2h2_model,select_ammonia_model,select_water_model,include_clouds);
done_with_abs_coeff=0

%  Check Boresight
if CRITICALFLAG==1;
   disp('Critical Refraction')
   disa=d;
   bs=1;
   [tau_a,tau]=ftau(kappa,d,masterindexa);	% Comes out with deepest first
   [Tatma,wlayersa]=ftam(T,tau,tau_a,masterindexa);
   sura=intercept(1,:);
   intercepts_boresight=intercept;
   CRITICALFLAG=0;
else
   [tau_a,tau]=ftau(kappa,d,masterindexa);	% Comes out with deepest first
   disa=d;
   bs=1;
   [Tatma,wlayersa]=ftam(T,tau,tau_a,masterindexa);
   sura=intercept(1,:);
   intercepts_boresight=intercept;
end   

% BORESIGHT DIDN'T MISS
% Now its worth calculating kappas to get taus

% d is the pathlength, t is theta, masterindex are the indices of the raypaths
% Need to gets tau's
% tau_a is the tau of that layer, tau is the cumulative summations of those layers
% Now do the beamspread and rotate beampattern to lookvector

[b1,b2,b3,b1w,b2w,b3w]=beamsample(N_ring_one,BWHM);

[Vr1,Vr2,Vr3,Zr]=rotbeam(Raydirection,b1,b2,b3);

% initialize wlayersa,b,c,d
wlayersb=zeros(length(P),length(Vr1));
wlayersc=zeros(length(P),length(Vr2));
wlayersd=zeros(length(P),length(Vr3));
% masterindex for wght fnctns
windexb=zeros(length(P),length(Vr1));
windexc=zeros(length(P),length(Vr2));
windexd=zeros(length(P),length(Vr3));
missb=0;
missc=0;
missd=0;

VRONE=0
for p=1:length(Vr1);
   bs=bs+1;						% bs-beamspread- keeps track of beamspread samples
   Rd=[Vr1(:,p)]';
   [intercept,internormal,d,t,masterindexb,missflag]=findraypath(recordlength,refindex,P,ellipses,Rayorigin,Rd);
   windexb(1:size(masterindexb,1),p)=masterindexb;
   disb=d;
   if missflag==1
      disp('This beam-sample  (b)Misses Planet')
      Tatmb(p)=2.7;
      wlayersb(:,p)=0;
      missb=missb+1;
      surb(p,:)=[0 0 0];
   elseif CRITICALFLAG==1;
      disp('Critical Refraction')
      [tau_a,tau]=ftau(kappa,d,masterindexb);	% Comes out with deepest first
      [Tatmb(p),wtemp]=ftam(T,tau,tau_a,masterindexb);
      wlayersb(1:size(wtemp,1),p)=wtemp;
      surb(p,:)=intercept(1,:);
      CRITICALFLAG=0;
   else
      [tau_a,tau]=ftau(kappa,d,masterindexb);	% Comes out with deepest first
      [Tatmb(p),wlayersb(:,p)]=ftam(T,tau,tau_a,masterindexb);
      surb(p,:)=intercept(1,:);
      intercepts_b(p,:,:)=intercept;
      
   end
end
VRTWO=0
for p=1:length(Vr2)
   bs=bs+1;						% bs-beamspread- keeps track of beamspread samples
   Raydirection=[Vr2(:,p)]';
   [intercept,internormal,d,t,masterindexc,missflag]=findraypath(recordlength,refindex,P,ellipses,Rayorigin,Raydirection);
   windexc(1:size(masterindexc,1),p)=masterindexc;
   disc=d;
   if missflag==1
      disp('This beam-sample(c) Misses Planet')
      Tatmc(p)=2.7;
      wlayersc(:,p)=0;
      crap=12;
      surc(p,:)=[0 0 0];
      missc=missc+1;
   elseif CRITICALFLAG==1;
      disp('Critical Refraction')
      [tau_a,tau]=ftau(kappa,d,masterindexc);	% Comes out with deepest first
      [Tatmc(p),wtemp]=ftam(T,tau,tau_a,masterindexc);
      wlayersc(1:size(wtemp,1),p)=wtemp;
      surc(p,:)=intercept(1,:);
      CRITICALFLAG=0;
   else
      [tau_a,tau]=ftau(kappa,d,masterindexc);	% Comes out with deepest first
      [Tatmc(p),wlayersc(:,p)]=ftam(T,tau,tau_a,masterindexc);
      surc(p,:)=intercept(1,:);
      intercepts_c(p,:,:)=intercept;
   end
end
VRTHREE=0   
for p=1:length(Vr3)
   bs=bs+1;						% bs-beamspread- keeps track of beamspread samples
   Raydirection=[Vr3(:,p)]';
   [intercept,internormal,d,t,masterindexd,missflag]=findraypath(recordlength,refindex,P,ellipses,Rayorigin,Raydirection);
   windexd(1:size(masterindexd,1),p)=masterindexd;
   disd=d;
   if missflag==1
      disp('This beam-sample (d) Misses Planet')
      Raydirection
      Tatmd(p)=2.7;
      wlayersd(:,p)=0;
      surd(p,:)=[0 0 0];
      missd=missd+1;

   elseif CRITICALFLAG==1;
      disp('Critical Refraction')
      [tau_a,tau]=ftau(kappa,d,masterindexd);	% Comes out with deepest first
      [Tatmd(p),wtemp]=ftam(T,tau,tau_a,masterindexd);
      wlayersd(1:size(wtemp,1),p)=wtemp;
      surd(p,:)=intercept(1,:);
      CRITICALFLAG=0;
      
   else
      [tau_a,tau]=ftau(kappa,d,masterindexd);	% Comes out with deepest first
      [Tatmd(p),wlayersd(:,p)]=ftam(T,tau,tau_a,masterindexd);
      surd(p,:)=intercept(1,:);
      intercepts_d(p,:,:)=intercept;
   end
end

% Apply Beam weights (Beam coupling)
% T=sum(1,N)Tn*(exp(-2.76*(deltabw/3db)^2))
b1w
b2w
b3w
Na=length(Tatma);
Nb=length(Tatmb);
Nc=length(Tatmc);
Nd=length(Tatmd);

numa=sum(Tatma.*1);
dena=sum(Na*1);
numb=sum(Tatmb.*b1w);
denb=sum(Nb*b1w);
numc=sum(Tatmc.*b2w);
denc=sum(Nc*b2w);
numd=sum(Tatmd.*b3w);
dend=sum(Nd*b3w);

Twa=sum(Tatma)*1./Na;
Twb=sum(Tatmb)*b1w./Nb;
Twc=sum(Tatmc)*b2w./Nc;
Twd=sum(Tatmd)*b3w./Nd;
Tbeam=(Twa+Twb+Twc+Twd)/(1+b1w+b2w+b3w);

% Find weighting function
wfa=1.*fwght(wlayersa,masterindexa);
wfb=b1w.*fwght(wlayersb,windexb);
wfc=b2w.*fwght(wlayersc,windexc);
wfd=b3w.*fwght(wlayersd,windexd);

% To make Matlab happy-make all the weight vectors the same size by padding with zeros
% to the size of the Pressure vector
sP=size(P,1);
sa=zeros(sP-size(wfa,1),1);
sb=zeros(sP-size(wfb,1),1);
sc=zeros(sP-size(wfc,1),1);
sd=zeros(sP-size(wfd,1),1);

% wfwa is already only 1 column
wfwa=[wfa;sa];
% turn into columns
wfwb=[wfb;sb];
wfwc=[wfc;sc];
wfwd=[wfd;sd];


sumb=sum(wfwb,2)./(size(Vr1,2)-missb);
sumc=sum(wfwc,2)./(size(Vr2,2)-missc);
sumd=sum(wfwd,2)./(size(Vr3,2)-missd);


weightingfactor=(wfwa+sumb+sumc+sumd)/(1+b1w+b2w+b3w);


% find out how many 
% example of how to index the Pressure with the wf (using plot)
% plot(P(1:max(masterindexa)),wfa)



% need to pad the wfw's somehow to be same length-need to average the
% multiple beams together
% need to take of care the zeros that come from limb sound-screws wieght

