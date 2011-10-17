#from mpl_toolkits.mplot3d import Axes3D
from pylab import *
from numpy import asarray,meshgrid

err=asarray([21.77855892,25.52909439,22.24444298,21.56800364,24.4307399,34.42853439,17.44370318,11.51954211,-75.48968604,188.86954016,140.85385185,
             1.73252564,-2.4471527,12.21311504,8.03455549,13.19146565,24.90801632,13.98620768,-20.2927057,-261.98795042,245.35561691,137.11501189,
             -5.8003426,-7.06379125,-0.84228769,-4.11926021,-0.69856761,8.24700249,-8.6885433,-33.64901993,-216.06298894,222.65611303,147.23348945,
             -21.34718412,-22.09205503,-15.56816577,-20.02854285,-11.76274356,-10.8141245,-36.85748222,-55.36739644,-195.72171344,213.58528675,154.56543618])
Data_B_x_h2_vals=asarray([0.01,0.03,0.06,0.09,0.12,0.2,0.25,0.4,0.6,0.8,0.9,
                          0.01,0.03,0.06,0.09,0.12,0.2,0.25,0.4,0.6,0.8,0.9,
                          0.01,0.03,0.06,0.09,0.12,0.2,0.25,0.4,0.6,0.8,0.9,
                          0.01,0.03,0.06,0.09,0.12,0.2,0.25,0.4,0.6,0.8,0.9])
Data_B_x_h2o_vals=1.0-Data_B_x_h2_vals
Data_B_T_vals=273.15+asarray([380,380,380,380,380,380,380,380,380,380,380,
                              400,400,400,400,400,400,400,400,400,400,400,
                              420,420,420,420,420,420,420,420,420,420,420,
                              440,440,440,440,440,440,440,440,440,440,440])
plot(Data_B_x_h2_vals[0:11],err[0:11],'rx')
plot(Data_B_x_h2_vals[11:22],err[11:22],'bx')
plot(Data_B_x_h2_vals[22:33],err[22:33],'cx')
plot(Data_B_x_h2_vals[33:44],err[33:44],'mx')
show()

