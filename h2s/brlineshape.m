function F = brlineshape(f,fo,df,ce,pst)

% This function determines the Ben-Reuven (BR) lineshape over
% a range of freqencies which are user defined as is the frequency step
% the resulting vector is passed to this subroutine as f.
% fo is a vector containing the center frequencies (GHz) of the spectral 
% lines of the given molecule (ie PH3) as supplied by the "Submillimeter,
% Millimeter, and Microwave Spectral Catalog" aka Poynter-Pickett Line Catalog
% but in GHz, not in MHz as provided by the catalog.
% Actually converted to inverse cm
% df is elsewhere called d_nu in the main program and is the line half-width at
% half-maximum or simply linewidth
% ce is the coupling element  which is calc by function
% pst is the pressure shift term which is calc by function
% (BR) Becomes (VVW) if ce=pst=0
% (BR) Becomes Gross or kinetic if df=ce and pst=0


% NOTE: All frequencies are or have been converted to inverse cm!!!

% fbr(f,fo,df,ce,pst)=(2/pi)*(f/fo)^2 * ( (df-ce)f^2 + (df+ce)*[(fo+pst)^2 + df^2 - ce^2])
%                                        -------------------------------------------------
%                                        [f^2 - (fo+pst)^2 - df^2 + ce^2]^2 + 4f^2df^2

% F =                          A      * (    B        +    C    *[           D          ])
%                                        ------------------------------------------------
%                                        [E  -     J       - G    +  H  ]^2 + I

n=size(f,1);
m=size(fo,2);
% f(1) f(1) f(1) f(1) ....f(n) n times where n is the number of frequency steps
% f(2) f(2) f(2) f(2)				in the observation range                            
% ...
% f(n) f(n) f(n) f(n)   n times where n is the number of frequency steps
% m times where m is the number of spectral lines

nones=ones(n,1);
mones=ones(1,m);
f_matrix=f*mones;
fo_matrix=nones*fo;
df_matrix=nones*df;
ce_matrix=nones*ce;
pst_matrix=nones*pst;


A=(2/pi)*(f_matrix./fo_matrix).^2;			% f(1)*fo(1)  f(2)*fo(1) f(3)*fo(1)
                                          % f(1)*fo(2)  f(2)*fo(2) f(3)*fo(2)
B=(df_matrix-ce_matrix).*f_matrix.^2;
C=df_matrix+ce_matrix;
D=((fo_matrix+pst_matrix).^2) + (df_matrix.^2)-(ce_matrix.^2);
E=f_matrix.^2;
J=(fo_matrix+pst_matrix).^2;
G=df_matrix.^2;
H=ce_matrix.^2;
I=4*(f_matrix.^2).*(df_matrix.^2);
F=df_matrix.*A.*(B+C.*D)./(((E-J-G+H).^2)+I);						% m x n matrix

% Note: the extra df factor comes in because of the expression for alpha has
% an additional df factor which cancels with the alpha_max, but is added here
% for computational reasons.  result is actually df*F
