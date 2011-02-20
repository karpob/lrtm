import matplotlib
matplotlib.use('PDF')
from scipy.optimize import leastsq
from residual_enthalpy_data import residual_enthalpy_data
from residual_pXT_data import residual_pXT_data
from calc_enthalpy import calc_enthalpy
from residual_enthalpy_and_PVT import residual_enthalpy_and_PVT
from Calc_B_from_B12 import Calc_B_from_B12
from pylab import *
from numpy import *

MPa_to_kPa=1000.0
Bars_to_kPa=100.0
#p=asarray([1,1,1,1,
#           1,1,1,1,
#           1,1,1,1,
#           1,1,1,1])
#[  7.63484836e-01,   4.46348033e-01,   1.15824901e+00,   1.25218138e+00,
#             6.87207340e-02,   3.46965438e-03,  -2.10640173e-02,  -5.98956557e-02,
#             7.12992049e+00,   6.78165187e+00,   4.45732550e+00,   6.23065727e+00,
#             8.60867222e+00,   1.05133196e+01,   9.99989620e-01,   8.79722224e+00]

p=asarray([ -3.00486044e+01,   6.78078409e-01 , -7.78933663e+01,   7.19303666e-01,
           -1.30262397e+00 ,  3.58122469e-02 ,  3.34024729e+00 , -9.06332331e+00,
             8.14449972e+00,   8.26908517e+00,  -4.29107942e-01,   8.44396579e+00,
             1.01431484e+00,   7.16209966e-01,   2.21352138e+00,   1.47909401e+00,
             1.73021814e-01,  -9.36093815e-02,   1.05574602e+00,   7.28203524e-01])





#p=asarray(  [ -6.84724158e+01,   2.76510561e+00,  -1.72902015e+02,   3.36805346e+00,
#               8.43730166e-02,   1.20304163e-02,   4.85353759e+00,  -9.45732780e+00,
#               2.96892622e+01,   5.66963126e+00,  -4.72763978e-01,   5.68600592e+00,
#               1.01325950e+00,   8.75427966e-01,   2.25904893e+00,   1.73721803e+00,
#               1.57106640e-01,  -1.23114242e-01,   1.07298418e+00,   7.51254725e-01] ) 

#
#p=asarray([ -6.90057247e+01,   3.00150652e+00,  -1.71869402e+02,   3.09456709e+00,
#             2.80153605e-01,   2.41600583e-02,   1.34205350e+00,  -9.07972171e+00,
#             3.06401744e+01,   9.15950512e+00,  -3.12792983e-01,   9.17167690e+00,
#             1.01324328e+00,   8.38454966e-01,   2.14176384e+00,   1.63739284e+00,
#             1.56816593e-01,  -1.21103695e-01,   4.45761096e+00,   2.41424542e+00])

p=asarray([ 6.04151604,  2.43018169, 14.79605951,  2.29450755,   
             0.17771438,   0.0241533 ,   3.8827558 , 0,   
             3.49625309,   9.16093547,  -0.29341299,   3.53438606,   
             1.00643321,   0.83846781,   1.91865226,    1.5718951,    
             1,  1,  1,   1])
             
#p=ones(shape(p))             
#[ -5.96721658e+00  -3.52604910e+00  -1.49843316e+01  -1.68222127e+00
#  -1.34736339e+01  -2.49967624e+00  -3.93781238e+00  -1.51898059e-06
#   2.68678174e+01   3.04649512e+01   4.61577583e+00  -3.71788456e-01
#   1.00000901e+00   1.00000380e+00   1.00000040e+00   9.99999798e-01
#   7.00056069e-01   1.00853445e+00   8.93148737e-01   9.04963671e-02] 
 
 
Data_Enthalpy_x_h2=asarray([0.5,0.5,0.5,0.5,0.5,
                            0.5,0.5,0.5,0.5,0.5,
                            0.5,0.5,0.5,0.5,0.5,
                            0.5,0.5,0.5,0.5,0.5,
                            0.5,0.5,0.5,0.5,0.5,
                            0.5,0.5,0.5,0.5,0.5,
                            0.5,0.5,0.5,0.5,0.5,
                            0.5,0.5,0.5,0.5,0.5,
                            0.5,0.5,0.5,0.5,0.5,
                            0.5,0.5,0.5,0.5,0.5,
                            0.5,0.5,0.5,0.5,0.5,
                            0.5,0.5,0.5,0.5,0.5,
                            0.5,0.69,0.588,0.481,0.371,
                            0.263,0.693,0.591,0.482,0.37,
                            0.257])
Data_Enthalpy_x_h2o=1.0-Data_Enthalpy_x_h2
Data_Enthalpy_T=asarray([448.2,473.2,498.2,523.2,548.2,
                         573.2,598.2,598.2,448.2,473.2,
                         498.2,523.2,548.2,573.2,598.2,
                         598.2,448.2,473.2,498.2,523.2,
                         548.2,573.2,598.2,598.2,448.2,
                         473.2,498.2,523.2,548.2,573.2,
                         598.2,598.2,448.2,473.2,498.2,
                         523.2,548.2,573.2,598.2,598.2,
                         648.2,648.2,648.2,648.2,648.2,
                         648.2,648.2,648.2,648.2,648.2,
                         648.2,698.2,698.2,698.2,698.2,
                         698.2,698.2,698.2,698.2,698.2,
                         698.2,598.2,598.2,598.2,598.2,
                         598.2,698.2,698.2,698.2,698.2,
                         698.2])
Data_Enthalpy_PMPa=asarray([0.55,0.68,0.48,1.06,1.47,
                            1.84,2.13,3.88,0.62,0.93,
                            0.7,2,2.93,3.42,5.38,6.88,
                            0.69,0.98,1.22,2.9,3.91,
                            4.84,6.93,8.2,0.76,1.28,
                            2.04,3.31,4.42,6.28,9.27,
                            9.32,0.93,1.34,2.14,3.55,
                            5.13,7.08,9.38,10.69,2.14,
                            4.01,5.88,7.96,9.75,9.89,
                            11.38,11.86,11.89,11.9,12.04,
                            2.09,2.93,4.05,4.96,5.85,
                            7.8,9.2,9.84,11.13,11.9,10.51,
                            10.51,10.51,10.51,10.51,11.13,
                            11.13,11.13,11.13,11.13])

Data_Enthalpy_H=asarray([163,160,75,182,225,
                         247,241,480,194,210,
                         111,386,489,488,710,
                         983,226,254,254,609,
                         706,758,1020,1269,247,
                         349,490,736,873,1103,
                         1555,1646,283,337,524,
                         750,1072,1356,1586,1971,
                         179,365,557,813,1039,
                         1077,1345,1456,1436,1438,
                         1472,149,173,292,321,
                         433,594,687,785,891,
                         992,1482,1783,1990,2025,
                         1789,743,852,941,878,
                         804])
                                                  
Data_Enthalpy_PkPa=Data_Enthalpy_PMPa*MPa_to_kPa
h2_M_amu=2.01594
h2o_M_amu=float(8.314472/0.46151805)
Data_pXT_T=asarray([374.3,374.5,374.5,374.5,374.5,374.5,374.5,375.3,376.5,379.0,381.3,374.3,374.5,374.5,374.5,374.5,374.5,374.5, 375.3, 376.5,379.0, 381.3,0,0,0])+273.15 #Celcius-->Kelvin
Data_pXT_P=asarray([229,237, 270,320,372,428, 690, 1010, 1410,2020,2520,229,237, 270,320,372,428, 690, 1010, 1410,2020,2520,0,0,0])*Bars_to_kPa # Bars->kPa
Data_pXT_x_h2=asarray([0.005,0.01,0.03,0.06,0.09,0.12,0.20,0.25,0.30,0.35,0.38,0.005,0.01,0.03,0.06,0.09,0.12,0.20,0.25,0.30,0.35,0.38,0,0,0])
Data_pXT_x_h2o=1.0-Data_pXT_x_h2 
Data_pXT_Mmix=((Data_pXT_x_h2o*h2o_M_amu)+(Data_pXT_x_h2*h2_M_amu))
Data_pXT_Vol=asarray([56.88,56.67,55.80,54.51,53.34,52.14,43.50,37.41,33.60,30.91,29.66,56.88,56.67,55.80,54.51,53.34,52.14,43.50,37.41,33.60,30.91,29.66,0,0,0])*(1.0/Data_pXT_Mmix)*(1.0/1.0e6)*1000.0 # cm**3/mol ->m**3/kg
                      #(cm**3/mol)(1 m**3/1e6 cm**3)(1 mol/Mmix(g))(1000 g/ 1 kg)
Data5050_Enthalpy_T=Data_Enthalpy_T[0:61]
Data5050_Enthalpy_P=Data_Enthalpy_PkPa[0:61]/Bars_to_kPa
Data5050_Data_Enthalpy_H=Data_Enthalpy_H[0:61]

sorted_index=Data5050_Enthalpy_T.argsort(axis=0)
Enthalpy_T_sorted=Data5050_Enthalpy_T[sorted_index]

Enthalpy_P_sorted=Data5050_Enthalpy_P[sorted_index]
Enthalpy_H_sorted=Data5050_Data_Enthalpy_H[sorted_index]

h2_M_amu=2.01594
h2o_M_amu=float(8.314472/0.46151805)
Data_pXT_T=asarray([374.3,374.5,374.5,374.5,374.5,374.5,374.5,375.3,376.5,379.0,381.3,374.3,374.5,374.5,374.5,374.5,374.5,374.5, 375.3, 376.5,379.0, 381.3,0,0,0])+273.15 #Celcius-->Kelvin
Data_pXT_P=asarray([229,237, 270,320,372,428, 690, 1010, 1410,2020,2520,229,237, 270,320,372,428, 690, 1010, 1410,2020,2520,0,0,0])*Bars_to_kPa # Bars->kPa
Data_pXT_x_h2=asarray([0.005,0.01,0.03,0.06,0.09,0.12,0.20,0.25,0.30,0.35,0.38,0.005,0.01,0.03,0.06,0.09,0.12,0.20,0.25,0.30,0.35,0.38,0,0,0])
Data_pXT_x_h2o=1.0-Data_pXT_x_h2 
Data_pXT_Mmix=((Data_pXT_x_h2o*h2o_M_amu)+(Data_pXT_x_h2*h2_M_amu))
Data_pXT_Vol=asarray([56.88,56.67,55.80,54.51,53.34,52.14,43.50,37.41,33.60,30.91,29.66,56.88,56.67,55.80,54.51,53.34,52.14,43.50,37.41,33.60,30.91,29.66,0,0,0])*(1.0/Data_pXT_Mmix)*(1.0/1.0e6)*1000.0 # cm**3/mol ->m**3/kg
                      #(cm**3/mol)(1 m**3/1e6 cm**3)(1 mol/Mmix(g))(1000 g/ 1 kg)
Data_pXT_density=1.0/Data_pXT_Vol #kg/m**3                      
Data_pXT_T[22]=376.204932    
Data_pXT_T[23]=377.726416    
Data_pXT_T[24]=446.697607      
     
Data_pXT_P[22]=19.826*Bars_to_kPa
Data_pXT_P[23]=75.059974199647826*Bars_to_kPa
Data_pXT_P[24]=87.860145913789353*Bars_to_kPa


Data_pXT_x_h2[22]=18.663980319129138/(18.663980319129138+1.126616830674388)#g/L equivalent to kg/m^3
Data_pXT_x_h2[23]=74.965164072252165/(74.965164072252165+1.126616830674388)#g/L equivalent to kg/m^3
Data_pXT_x_h2[24]=74.965164072252165/(74.965164072252165+1.126616830674388)#g/L equivalent to kg/m^3


Data_pXT_x_h2o[22]=1.0-Data_pXT_x_h2[22]
Data_pXT_x_h2o[23]=1.0-Data_pXT_x_h2[23]
Data_pXT_x_h2o[24]=1.0-Data_pXT_x_h2[24]


Data_pXT_density[22]=(18.663980319129138*h2_M_amu+1.126616830674388*h2o_M_amu)/32.328408091385590#g/L equivalent to kg/m^3
Data_pXT_density[23]=(74.965164072252165*h2_M_amu+1.126616830674388*h2o_M_amu)/32.337938041626330#g/L equivalent to kg/m^3
Data_pXT_density[24]=(74.965164072252165*h2_M_amu+1.126616830674388*h2o_M_amu)/32.413598694129391#g/L equivalent to kg/m^3                      

CC_per_mol_to_Liters_per_mol=0.001
Data_B_B_vals=CC_per_mol_to_Liters_per_mol*asarray([-72.3,-72.86,-65.5,-60.75,-58.8,-55.52,-38.62,-23.49,-6.26,5.9,8.16,-62.34,-57.42,-62.99,-56.35,-55.78,-53.14,-40.79,-19.24,-3.42,4.05,9.94,-62.57,-59.42,-59.38,-53.98,-52.24,-47.49,-35.37,-19.17,-4.38,5.35,8.6,-58.82,-56.22,-55.97,-50.65,-50.99,-42.79,-30.66,-18.16,-5.2,6.4,8.16])
#in C converted to K
Data_B_T_vals=273.15+asarray([380,380,380,380,380,380,380,380,380,380,380,400,400,400,400,400,400,400,400,400,400,400,420,420,420,420,420,420,420,420,420,420,420,440,440,440,440,440,440,440,440,440,440,440])
Data_B_x_h2_vals=asarray([0.01,0.03,0.06,0.09,0.12,0.2,0.25,0.4,0.6,0.8,0.9,0.01,0.03,0.06,0.09,0.12,0.2,0.25,0.4,0.6,0.8,0.9,0.01,0.03,0.06,0.09,0.12,0.2,0.25,0.4,0.6,0.8,0.9,0.01,0.03,0.06,0.09,0.12,0.2,0.25,0.4,0.6,0.8,0.9])
Data_B_x_h2o_vals=1.0-Data_B_x_h2_vals
Data_C_vals=CC_per_mol_to_Liters_per_mol*CC_per_mol_to_Liters_per_mol*asarray([1580.57,1578.49,1439.30,1395.68,1399.75,1511.79,1047.75,806.92,605.14,503.93,565.32,1281.16,1188.67,1420.97,1296.81,1345.21,1478.18,1166.11,710.47,528.30,566.33,496.10,1381.01,1299.78,1369.50,1297.59,1289.93,1309.09,1029.32,735.58,589.89,545.00,542.15,1315.16,1258.02,1320.81,1239.98,1345.80,1150.23,914.58,745.09,656.09,508.69,557.31])                       

#B12_vals=asarray([-13.69,-10.21,-7.38,-5.05,-3.10,-1.44,-0.02,1.22,2.29,3.24,4.08])
#B12_T=asarray([200.0,220.0,240.0,260.0,280.0,300.0,320.0,340.0,360.0,380.0,400.0])

#B_rabinovich=Calc_B_from_B12(B12_T,B12_vals,


[optimized_p,message]=leastsq(residual_enthalpy_and_PVT,p[:],args=(Data_Enthalpy_H,Data_Enthalpy_PkPa,Data_Enthalpy_T,Data_Enthalpy_x_h2,Data_Enthalpy_x_h2o,Data_pXT_P,Data_pXT_T,Data_pXT_density,Data_pXT_x_h2,Data_pXT_x_h2o,Data_B_B_vals,Data_B_T_vals,Data_B_x_h2_vals,Data_B_x_h2o_vals,Data_C_vals))
p=optimized_p
x_h2_598_2=asarray(arange(0.01,1.0,0.01))
x_h2_698_2=asarray(arange(0.01,1.0,0.01))


P448=asarray(arange(0.1,1.3,0.1)*MPa_to_kPa)
P473=asarray(arange(0.1,2,0.1)*MPa_to_kPa)
P498=asarray(arange(0.1,3,0.1)*MPa_to_kPa)
P523=asarray(arange(0.1,4,0.1)*MPa_to_kPa)
P548=asarray(arange(0.1,6,0.1)*MPa_to_kPa)
P573=asarray(arange(0.1,8.5,0.1)*MPa_to_kPa)
P598=asarray(arange(0.1,11,0.1)*MPa_to_kPa)
P648=asarray(arange(0.1,13,0.1)*MPa_to_kPa)
P698=asarray(arange(0.1,13,0.1)*MPa_to_kPa)

P_598_2=10.51*MPa_to_kPa*ones(shape(x_h2_598_2))
P_698_2=11.13*MPa_to_kPa*ones(shape(x_h2_698_2))

T448=448.2*ones(shape(P448))
T473=473.2*ones(shape(P473))
T498=498.2*ones(shape(P498))
T523=523.2*ones(shape(P523))
T548=548.2*ones(shape(P548))
T573=573.2*ones(shape(P573))
T598=598.2*ones(shape(P598))
T648=648.2*ones(shape(P648))
T698=698.2*ones(shape(P698))

T_598_2=598.2*ones(shape(P_598_2))
T_698_2=698.2*ones(shape(P_698_2))

x_h2_448=0.5*ones(shape(P448))
x_h2_473=0.5*ones(shape(P473))
x_h2_498=0.5*ones(shape(P498))
x_h2_523=0.5*ones(shape(P523))
x_h2_548=0.5*ones(shape(P548))
x_h2_573=0.5*ones(shape(P573))
x_h2_598=0.5*ones(shape(P598))
x_h2_648=0.5*ones(shape(P648))
x_h2_698=0.5*ones(shape(P698))

x_h2o_698_2=1.0-x_h2_698_2
x_h2o_598_2=1.0-x_h2_598_2

Curve_448=calc_enthalpy(p,P448,T448,x_h2_448,x_h2_448)
Curve_473=calc_enthalpy(p,P473,T473,x_h2_473,x_h2_473)
Curve_498=calc_enthalpy(p,P498,T498,x_h2_498,x_h2_498)
Curve_523=calc_enthalpy(p,P523,T523,x_h2_523,x_h2_523)
Curve_548=calc_enthalpy(p,P548,T548,x_h2_548,x_h2_548)
Curve_573=calc_enthalpy(p,P573,T573,x_h2_573,x_h2_573)
Curve_598=calc_enthalpy(p,P598,T598,x_h2_598,x_h2_598)
Curve_648=calc_enthalpy(p,P648,T648,x_h2_648,x_h2_648)
Curve_698=calc_enthalpy(p,P698,T698,x_h2_698,x_h2_698)

Curve_698_2=calc_enthalpy(p,P_698_2,T_698_2,x_h2_698_2,x_h2o_698_2)
Curve_598_2=calc_enthalpy(p,P_598_2,T_598_2,x_h2_598_2,x_h2o_598_2)

figure(1)
plot(P448/Bars_to_kPa,Curve_448,'r--',label='448 K')
plot(P473/Bars_to_kPa,Curve_473,'g--',label='473 K')
plot(P498/Bars_to_kPa,Curve_498,'b',label='498 K')
plot(P523/Bars_to_kPa,Curve_523,'k',label='523 K')

plot(P548/Bars_to_kPa,Curve_548,'c', label='548 K')
plot(P573/Bars_to_kPa,Curve_573,'m', label='573 K')
plot(P598/Bars_to_kPa,Curve_598,'y', label='598 K')
plot(P648/Bars_to_kPa,Curve_648,'r', label='648 K')
plot(P698/Bars_to_kPa,Curve_698,'g', label='698 K')

errorbar(Enthalpy_P_sorted[0:5],Enthalpy_H_sorted[0:5],yerr=0.02*Enthalpy_H_sorted[0:5],ecolor='r',fmt=None)
errorbar(Enthalpy_P_sorted[5:10],Enthalpy_H_sorted[5:10],yerr=0.02*Enthalpy_H_sorted[5:10],ecolor='g',fmt=None)
errorbar(Enthalpy_P_sorted[10:15],Enthalpy_H_sorted[10:15],yerr=0.02*Enthalpy_H_sorted[10:15],ecolor='b',fmt=None)
errorbar(Enthalpy_P_sorted[15:20],Enthalpy_H_sorted[15:20],yerr=0.02*Enthalpy_H_sorted[15:20],ecolor='k',fmt=None)

errorbar(Enthalpy_P_sorted[20:25],Enthalpy_H_sorted[20:25],yerr=0.02*Enthalpy_H_sorted[20:25],ecolor='c',fmt=None)
errorbar(Enthalpy_P_sorted[25:30],Enthalpy_H_sorted[25:30],yerr=0.02*Enthalpy_H_sorted[25:30],ecolor='m',fmt=None)
errorbar(Enthalpy_P_sorted[30:40],Enthalpy_H_sorted[30:40],yerr=0.02*Enthalpy_H_sorted[30:40],ecolor='y',fmt=None)
errorbar(Enthalpy_P_sorted[40:51],Enthalpy_H_sorted[40:51],yerr=0.02*Enthalpy_H_sorted[40:51],ecolor='r',fmt=None)
errorbar(Enthalpy_P_sorted[51:61],Enthalpy_H_sorted[51:61],yerr=0.02*Enthalpy_H_sorted[51:61],ecolor='g',fmt=None)
legend(loc='upper left')
xlabel('Pressure (bars)')
ylabel('Excess Enthalpy (J/Mol)')
title('Excess Enthalpy 50/50 H$_2$/H$_2$O')
savefig('Excess_Enthalpy_5050.pdf',papertype='letter')

figure(2)
plot(1-x_h2o_598_2,Curve_598_2,'b',label='598.2 K, 10.51 MPa')
plot(1-x_h2o_698_2,Curve_698_2,'r',label='698.2 K, 11.13 MPa' )
errorbar(Data_Enthalpy_x_h2[61:66],Data_Enthalpy_H[61:66],yerr=0.02*Data_Enthalpy_H[61:66],ecolor='b',fmt=None)
errorbar(Data_Enthalpy_x_h2[66:len(Data_Enthalpy_x_h2)],Data_Enthalpy_H[66:len(Data_Enthalpy_x_h2)],yerr=0.02*Data_Enthalpy_H[66:len(Data_Enthalpy_x_h2)],ecolor='r',fmt=None)
legend(loc='upper right')
xlabel('Mole Fraction X$_{H_2}$')
ylabel('Excess Enthalpy (J/mol)')
title('Excess Enthalpy for Increasing X$_{H_2}$')
savefig('Excess_Enthalpy_mole_frac.pdf',papertype='letter')
residuals=residual_enthalpy_and_PVT(p[:],Data_Enthalpy_H,Data_Enthalpy_PkPa,Data_Enthalpy_T,Data_Enthalpy_x_h2,Data_Enthalpy_x_h2o,Data_pXT_P,Data_pXT_T,Data_pXT_density,Data_pXT_x_h2,Data_pXT_x_h2o,Data_B_B_vals,Data_B_T_vals,Data_B_x_h2_vals,Data_B_x_h2o_vals,Data_C_vals)

figure(3)
plot(Data_pXT_T,residuals[0:len(Data_pXT_P)],'kx')
xlabel('Temperature ($^{\circ}$K)')
ylabel('Error in Pressure (%)')
title('Residual Pressure from EOS Fit')
savefig('Residual_Pressure.pdf',papertype='letter')

begin=len(Data_pXT_P)+len(Data_Enthalpy_H)
end=begin+len(Data_B_T_vals)

figure(4)
plot(Data_B_T_vals,residuals[begin:end],'kx')
xlabel('Temperature ($^{\circ}$K)')
ylabel('Second Virial Coefficient Error (%)')
title('Residual Second Virial Coefficient from EOS Fit')
savefig('Residual_B.pdf',papertype='letter')

begin=end
end=begin+len(Data_B_T_vals)
figure(5)
plot(Data_B_T_vals,residuals[begin:end],'kx')
xlabel('Temperature ($^{\circ}$K)')
ylabel('Third Virial Coefficient Error (%)')
title('Residual Third Virial Coefficient from EOS Fit')
savefig('Residual_C.pdf',papertype='letter')


