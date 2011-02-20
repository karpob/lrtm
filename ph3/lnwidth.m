function delta_nu = lnwidth(To,T,Mol_1,dP_1,Mol_2,dP_2,Mol_3,dP_3)
% To is reference Temperature-usually 300 k
% T is actual temperature
% Mol's is for extending functionality to use more than just H2,He,PH3
% dP is the partial pressure of each molecule contribution

% Calculates the linewidth for a gas mixture of mixtures (i) and produces
% a delta_nu line broadening in GHz for each line (j)
% delta_nu will be a vector of length (j)
% delta_nu(j)=sum{i}(delta_nu(i,j)*dP(i)*(To/T)^-(epsilson)
% where epsilon is (m+1)/2(m-1)   m is related to force law 
% m=3 dipole, m=infi is hardball


% LOAD in CORRECT PRESSURE BROADENING FILES	
if(strcmp(Mol_1,'H2'))
   load h2pbroad
   M1=h2pbroad;
   m(1)=7;										% force law exponent
elseif(strcmp(Mol_1,'He'))
   load hepbroad
   M1=hepbroad;
   m(1)=7;										% force law exponent
elseif(strcmp(Mol_1,'PH3'))
   load ph3pbroad
   M1=ph3pbroad;
   m(1)=3;										% force law exponent
else
   error('Unknown Molecule label.');
end


if(strcmp(Mol_2,'H2'))
   load h2pbroad
   M2=h2pbroad;
	m(2)=7;										% force law exponent   
elseif(strcmp(Mol_2,'He'))
   load hepbroad
   M2=hepbroad;
   m(2)=7;										% force law exponent
elseif(strcmp(Mol_2,'PH3'))
   load ph3pbroad
   M2=ph3pbroad; 
   m(2)=3;										% force law exponent
else
   error('Unknown Molecule label.');
end

if(strcmp(Mol_3,'H2'))
   load h2pbroad
   M3=h2pbroad;
   m(3)=7;										% force law exponent
elseif(strcmp(Mol_3,'He'))
   load hepbroad
   M3=hepbroad;
   m(3)=7;										% force law exponent
elseif(strcmp(Mol_3,'PH3'))
   load ph3pbroad
   M3=ph3pbroad;
   %m(3)=3;	% force law exponent
   m(3)=3;
else
   error('Unknown Molecule label.');
end

% LOADING COMPLETE
%a1=[15 17 19 21 23 25 26 28 30 32 33 34 35 36 37 38 39 40 41];
%a2=[2:14 16 18 20 22 24 27 29 31];
%a1=a1-1;
%a2=a2-1;

%a1=[1 4 7 9 11 14 16 18 20 22 24 25 27 29 31:40];
%a2=[2 3 5 6 8 10 12 13 15 17 19 21 23 26 28 30]; 


%M1(a1)=M1(a1).*A(3);
%M1(a2)=M1(a2).*A(4);
%M1(41:105)=M1(41:105).*A(6);

%M2(a1)=M2(a1).*A(3);
%M2(a2)=M2(a2).*A(4);
%M2(41:105)=M2(41:105).*A(6);

%M3(a1)=M3(a1).*A(5);
%M3(a2)=M3(a2).*A(6);
%M3(41:105)=M3(41:105).*A(9);

i=1:3;											% Three gas mixture
epsilon(i)=(m(i)+1)./(2*(m(i)-1));
dnu1=M1.*dP_1.*(To/T).^(epsilon(1));
dnu2=M2.*dP_2.*(To/T).^(epsilon(2));
dnu3=M3.*dP_3.*(To/T).^(epsilon(3));
delta_nu=dnu1+dnu2+dnu3;					% GHz
