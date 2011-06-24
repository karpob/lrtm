#from alpha_normal_h2 import alpha_normal_h2
from pylab import shape
from pylab import zeros
from enthalpy_methane import enthalpy_methane
#from helmholtz_functions.alpha_normal_h2 import alpha_normal_h2
	#/* Temperature, Pressure, Density, Int. Energy, Enthalpy, Entropy, Cv, Cp, Cp0, Helmholtz */
	#/* (C), (MPa), (kg/m^3), (kJ/kg), (kJ/kg), (kJ/kg-K), (kJ/kg-K), (kJ/kg-K), (kJ/kg-K), (kJ/kg) */
	
#	Temperature	Pressure	Density	Int. Energy	Enthalpy	Entropy	Cv	Cp	Cp0	Helmholtz
									

[m,n]=shape(reference_data)
absolute=zeros(m)
#normalized=absolute
ideal=zeros(m)
residual=zeros(m)
h=zeros(m)
Z=zeros(m)
s=zeros(m)
cp=zeros(m)
cv=zeros(m)
speed_o_sound=zeros(m)
R=8.31451
R_methane=R/16.0428
#print 'Z'+'    '+'Enthalpy'+'    '+'Entropy'+'    '+'Cv'+'    '+'Cp'+'    '+'Helmholtz\n' 
#print 1e3*reference_data[i][1]/(reference_data[i][2]*R_methane*(reference_data[i][0]+273.15))-Z[i]
for i in range(0,len(reference_data)):
        [Z[i],h[i],s[i],absolute[i],cv[i],cp[i],speed_o_sound[i],ideal[i],residual[i]]=enthalpy_methane(float(reference_data[i][0])+273.15,float(reference_data[i][2]))
#        print reference_data[i][0],'\t' ,1e3*reference_data[i][1]/(reference_data[i][2]*R_methane*(reference_data[i][0]+273.15))-Z[i],'\t',reference_data[i][4]-h[i],'\t',s[i]-reference_data[i][5],'\t', cv[i]-reference_data[i][6],'\t',cp[i]-reference_data[i][7],'\t',absolute[i]-reference_data[i][9]
#print '\n'        
#print 'Z\t\t\t','Enthalpy\t\t','Entropy\t\t','Cv\t\t\t','Cp\t\t\t','Helmholtz\n'         
#for i in range(0,m):
#        print Z[i],'\t',h[i],'\t',s[i],'\t', cv[i],'\t',cp[i],'\t',absolute[i],speed_o_sound[i],ideal[i],residual[i]