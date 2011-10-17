from pylab import shape
from pylab import zeros
from enthalpy_h2o import enthalpy_h2o
## T   rho   P  cv   w   s
reference_data=[[300.0, 0.9965560e3,0.992418352e-1,0.413018112e1,0.150151914e4,0.393062643],
 [300.0,0.1005308e4,0.200022515e2,0.406798347e1,0.153492501e4,0.387405401],
 [300.0,0.1188202e4,0.700004704e3,0.346135580e1,0.244357992e4,0.132609616],
 [500.0,0.4350000,0.999679423e-1,0.150817541e1,0.548314253e3,0.794488271e1],
 [500.0,0.4532000e1,0.999938125,0.166991025e1,0.535739001e3,0.682502725e1],
 [500.0,0.8380250e3,0.100003858e2,0.322106219e1,0.127128441e4,0.256690919e1],
 [500.0,0.1084564e4,0.700000405e3,0.307437693e1,0.241200877e4,0.203237509e1],
 [647.00000000000000000,358.0000,0.220384756e2,0.618315728e1,0.252145078e3,0.432092307e1],   
 [900.0,0.2410000,0.100062559,0.175890657e1,0.724027147e3,0.916653194e1],
 [900.0,0.5261500e2,0.200000690e2,0.193510526e1,0.698445674e3,0.659070225e1],
 [900.0,0.8707690e3,0.700000006e3,0.266422350e1,0.201933608e4,0.417223802e1]]
[m,n]=shape(reference_data) 
R_h2o = 0.46151805  #kJ kg^-1 K^-1
absolute=zeros(m)
h=zeros(m)
Z=zeros(m)
s=zeros(m)
cp=zeros(m)
cv=zeros(m)
speed_o_sound=zeros(m)

#print 'Z'+'    '+'Enthalpy'+'    '+'Entropy'+'    '+'Cv'+'    '+'Cp'+'    '+'Helmholtz\n' 
for i in range(0,m):
        [Z[i],h[i],s[i],absolute[i],cv[i],cp[i],speed_o_sound[i]]=enthalpy_h2o(float(reference_data[i][0]),float(reference_data[i][1]))
        #print Z[i]#reference_data[i][0],reference_data[i][1],reference_data[i][2],reference_data[i][3],reference_data[i][4],reference_data[i][5]#cv[i]#,h[i],s[i],absolute[i],cv[i],cp[i],speed_o_sound[i]
        
#        print (1e3*reference_data[i][2]/(R_h2o*reference_data[i][1]*reference_data[i][0]))-(Z[i]),cv[i]-reference_data[i][3],'\t', speed_o_sound[i]-reference_data[i][4],'\t',s[i]-reference_data[i][5]
        
