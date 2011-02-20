function alphanh3=devaraj_nh3_model_0(f,T,P,H2mr,Hemr,NH3mr)

load ammonia_lines
%Loads the Poynter Pickett (JPL) nh3 line catalog.  The workspace contains the arrays for Freq (fo) in GHz, intensity (Io) in inverse cm,
%lower state energy (Eo) in inverse cm and self broadening due to NH3 (gammaNH3o) in GHz/bar.
%Constants
GHztoinv_cm=1/29.9792458;	            %for converting GHz to inverse cm
OpticaldepthstodB=434294.5;				%converts from cm^-1 to dB/km
torrperatm=760;
bartoatm=0.987;
GHztoMHz=1000;
hc=19.858252418E-24;			%planks (J.s) light (cm/s)
k=1.38*10^-23;					%boltzmann's in J/K or N.m/K
No=6.02297e23;					%Avogadros Number [mole^-1]
R=8.31432e7;					%Rydberg's [erg/mole-K]
To=300;			                %Ref temp for P/P Catalogue
dynesperbar=1e6;				%dyne=bar/1e6;
coef=dynesperbar*No/R;          %See Appendix D: Using the Poyter-Pickett Catalogs

variable_hanley=[  1.2817    1.4867      0.5698    0.9336    0.3060   -0.3000    0.7902    0.2014    2.3513];
    variable_rot=[0.2000    1.9960    2.0000    0.2002     2.0000   -3.0000    3.0    2.0000    2.0000];
    variable_v2=[0.2000    2.0000 0.0001     0.9243  0.5000    0.2641    0.2148     2.0000    1.1000];

    PH2=P*H2mr; %partial pressures
    PHe=P*Hemr;
    PNH3=P*NH3mr;
    %calculate vector linewidth
    xi1=variable_hanley(1);
    xi2=2/3;
    xi3=1;
    xi12=variable_hanley(2);
    xi22=2/3;
    xi32=variable_hanley(3); 
    Tdiv=To/T;
    gnu1=variable_hanley(4); 
    gnu2=0.75; 
    gnu3=variable_hanley(5); 
    gH2=gnu1*PH2;
    gHe=gnu2*PHe;
    gNH3=gnu3*PNH3*gammaNH3o;
    gamma=((gH2)*((Tdiv)^(xi1))+(gHe)*((Tdiv)^(xi2))+gNH3*(295/T)^(xi3));

    delt=variable_hanley(6)*gamma; 
    znu1=variable_hanley(7);
    znu2=0.3;
    znu3=variable_hanley(8); 
    zH2=znu1*PH2;
    zHe=znu2*PHe;
    zNH3=znu3*PNH3*gammaNH3o;
    zeta=(zH2)*((Tdiv)^(xi12))+(zHe)*((Tdiv)^(xi22))+zNH3*(295/T)^(xi32);

    zetasize=size(fo,1);
    pst=delt;       							% answer in GHz
    %Coupling element, pressure shift and dnu or gamma are in GHz, need to convert brlineshape to inverse cm which is done below

    n=size(f,2);  %returns the number of columns in f
    m=size(fo,1); %returns the number of rows in fo
    % f1 f2 f3 f4 ....fn  n times where n is the number of frequency steps
    % f1 f2 f3 f4 ....fn				in the observation range
    % ...
    % f1 f2 f3 f4 ....fn
    % m times where m is the number of spectral lines

    nones=ones(1,n);
    mones=ones(m,1);
    f_matrix=mones*f;
    fo_matrix=fo*nones;

    % The 10^6 allows use of P(bar) for P(dynes/cm^2)

    eta=3/2;			% for symmetric top molecule
    expo=-(1/T-1/To)*Eo*hc/k;
    ST=Io.*exp(expo);	% S(T) =S(To)converted for temperature
    Con=variable_hanley(9);%(0.9349+PH2/T*0.5388-(PH2/T)^2*0.102);%(0.9349+PH2/T*0.535+(PH2/T)^2*0.0312);%(0.9404);
    alpha_noshape=Con*coef*(PNH3/To)*((To/T)^(eta+2)).*ST;%0.9387
    %Alpha Max Found

    %Ben Reuven lineshape calculated by the brlineshape function gives the answer in GHz
    %Here we change from GHz to inverse cm.
    dnu_matrix=gamma*nones;
    ce_matrix=zeta*nones;
    pst_matrix=pst*nones;
    Aa=(2/pi)*((f_matrix./fo_matrix).^2);
    Bb=(dnu_matrix-ce_matrix).*(f_matrix.^2);
    Cc=dnu_matrix+ce_matrix;
    Dd=((fo_matrix+pst_matrix).^2) + (dnu_matrix.^2)-(ce_matrix.^2);
    Ee=f_matrix.^2;
    Jj=(fo_matrix+pst_matrix).^2;
    Gg=dnu_matrix.^2;
    Hh=ce_matrix.^2;
    Ii=4*(f_matrix.^2).*(dnu_matrix.^2);
    Ff=Aa.*(Bb+Cc.*Dd)./(((Ee-Jj-Gg+Hh).^2)+Ii);

    Fbr=(1/GHztoinv_cm).*Ff;

    alpha_noshape_matrix=alpha_noshape*nones;
    
    br_alpha_matrix=alpha_noshape_matrix.*Fbr;



    % Computing the absorption contributed by the rotational lines
    ST_rot=Io_rot.*(exp((1/To-1/T)*Eo_rot*hc/k));
    
     %Factor GAMMA:
    xi1=variable_rot(1); xi2=2/3; xi3=1;
    %Factor Z:
    xi12=variable_rot(2); xi22=2/3; xi32=variable_hanley(3);
    %Factor temperature:
    Tdiv=To/T;
    %Factor nu:
    gnu1=variable_rot(4); gnu2=0.75; gnu3=variable_hanley(5);
    %Factor gamma:
    gH2=gnu1*PH2*gH2_rot; gHe=gnu2*PHe*gHe_rot;gNH3=gnu3*PNH3*gNH3_rot;
    %Total broadening
    gamma=((gH2)*((Tdiv)^(xi1))+(gHe)*((Tdiv)^(xi2))+gNH3*(295/T)^(xi3));
    %Factor delta
    delt=variable_rot(6)*gamma;
    %Factor zeta
    znu1=variable_rot(7); znu2=0.3; znu3=variable_rot(8);
    %Factor z
    zH2=znu1*PH2*gH2_rot; zHe=znu2*PHe*gHe_rot;zNH3=znu3*PNH3*gNH3_rot;
    %Factor zeta:
    zeta=(zH2)*((Tdiv)^(xi12))+(zHe)*((Tdiv)^(xi22))+zNH3*(295/T)^(xi32);

    zetasize=size(fo_rot,1);
    pst=delt;							% answer in GHz
    
    %Coupling element, pressure shift and dnu or gamma are in GHz,
    % need to convert brlineshape to inverse cm which is done below
    
    n=size(f,2);  %returns the number of columns in f
    m=size(fo_rot,1); %returns the number of rows in fo
    nones=ones(1,n);
    mones=ones(m,1);
    f_matrix=mones*f;
    fo_matrix=fo_rot*nones;
    % The 10^6 allows use of P(bar) for P(dynes/cm^2)
    eta=3/2;			% for symmetric top molecule
    expo=-(1/T-1/To)*Eo_rot*hc/k;
    ST=Io_rot.*exp(expo);	% S(T) =S(To)converted for temperature
    Con=variable_rot(9);
    alpha_max_rot=Con*coef*(PNH3/To)*((To/T)^(eta+2)).*ST;
    %Alpha Max Found

    %Ben Reuven lineshape calculated by the brlineshape function gives the
    % answer in GHz. Here we change from GHz to inverse cm.
    dnu_matrix=gamma*nones;
    ce_matrix=zeta*nones;
    pst_matrix=pst*nones;
    Aa=(2/pi)*((f_matrix./fo_matrix).^2);
    Bb=(dnu_matrix-ce_matrix).*(f_matrix.^2);
    Cc=dnu_matrix+ce_matrix;
    Dd=((fo_matrix+pst_matrix).^2) + (dnu_matrix.^2)-(ce_matrix.^2);
    Ee=f_matrix.^2;
    Jj=(fo_matrix+pst_matrix).^2;
    Gg=dnu_matrix.^2;
    Hh=ce_matrix.^2;
    Ii=4*(f_matrix.^2).*(dnu_matrix.^2);
    Ff=Aa.*(Bb+Cc.*Dd)./(((Ee-Jj-Gg+Hh).^2)+Ii);
    Fbr=(1/GHztoinv_cm).*Ff;
    alpha_rot_matrix=alpha_max_rot*nones;
    alpha_rot=alpha_rot_matrix.*Fbr;

    % Computing the absorption contributed by the v2 rotovibrational lines
    ST_v2=Io_v2.*(exp((1/To-1/T)*Eo_v2*hc/k));
    Io_v2(1)=Io_v2(1)*1.25;
 % Broadening parameters for the v2 vibrational lines
    gH2_v2=linspace(0.7069,0.7069,length(fo_v2))';
    gHe_v2=linspace(0.2000,0.2000,length(fo_v2))';
    gNH3_v2=linspace(13.2037, 13.2037,length(fo_v2))';
    
     %Factor GAMMA:
    xi1=variable_v2(1); xi2=2/3; xi3=1;
    %Factor Z:
    xi12=variable_v2(2); xi22=2/3; xi32=variable_v2(3);
    %Factor temperature:
    Tdiv=To/T;
    %Factor nu:
    gnu1=variable_v2(4); gnu2=0.75; gnu3=variable_v2(5);
    %Factor gamma:
    gH2=gnu1*PH2*gH2_v2; gHe=gnu2*PHe*gHe_v2;gNH3=gnu3*PNH3*gNH3_v2;
    %Total broadening
    gamma=((gH2)*((Tdiv)^(xi1))+(gHe)*((Tdiv)^(xi2))+gNH3*(295/T)^(xi3));
    %Factor delta
    delt=variable_v2(6)*gamma;
    %Factor zeta
    znu1=variable_v2(7); znu2=0.3; znu3=variable_v2(8);
    %Factor z
    zH2=znu1*PH2*gH2_v2; zHe=znu2*PHe*gHe_v2;zNH3=znu3*PNH3*gNH3_v2;
    %Factor zeta:
    zeta=(zH2)*((Tdiv)^(xi12))+(zHe)*((Tdiv)^(xi22))+zNH3*(295/T)^(xi32);

    zetasize=size(fo_v2,1);
    pst=delt;							% answer in GHz
    
    %Coupling element, pressure shift and dnu or gamma are in GHz,
    % need to convert brlineshape to inverse cm which is done below
    
    n=size(f,2);  %returns the number of columns in f
    m=size(fo_v2,1); %returns the number of rows in fo
    nones=ones(1,n);
    mones=ones(m,1);
    f_matrix=mones*f;
    fo_matrix=fo_v2*nones;
    % The 10^6 allows use of P(bar) for P(dynes/cm^2)
    eta=3/2;			% for symmetric top molecule
    expo=-(1/T-1/To)*Eo_v2*hc/k;
    ST=Io_v2.*exp(expo);	% S(T) =S(To)converted for temperature
    Con=variable_v2(9);
    alpha_max_v2=Con*coef*(PNH3/To)*((To/T)^(eta+2)).*ST;
    %Alpha Max Found

    %Ben Reuven lineshape calculated by the brlineshape function gives the
    % answer in GHz. Here we change from GHz to inverse cm.
    dnu_matrix=gamma*nones;
    ce_matrix=zeta*nones;
    pst_matrix=pst*nones;
    Aa=(2/pi)*((f_matrix./fo_matrix).^2);
    Bb=(dnu_matrix-ce_matrix).*(f_matrix.^2);
    Cc=dnu_matrix+ce_matrix;
    Dd=((fo_matrix+pst_matrix).^2) + (dnu_matrix.^2)-(ce_matrix.^2);
    Ee=f_matrix.^2;
    Jj=(fo_matrix+pst_matrix).^2;
    Gg=dnu_matrix.^2;
    Hh=ce_matrix.^2;
    Ii=4*(f_matrix.^2).*(dnu_matrix.^2);
    Ff=Aa.*(Bb+Cc.*Dd)./(((Ee-Jj-Gg+Hh).^2)+Ii);
    Fbr=(1/GHztoinv_cm).*Ff;
    alpha_v2_matrix=alpha_max_v2*nones;
    alpha_v2=alpha_v2_matrix.*Fbr;
     
    alpha_opdep=sum(br_alpha_matrix,1)+sum(alpha_rot,1)+sum(alpha_v2,1);
    alphanh3=alpha_opdep*434294.5;



