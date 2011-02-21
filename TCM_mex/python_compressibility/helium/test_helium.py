#from alpha_normal_h2 import alpha_normal_h2
from pylab import shape
from pylab import zeros
from enthalpy_helium import enthalpy_helium

	#Temperature (C),Pressure (MPa),Density (kg/m3),Volume (m3/kg),Internal Energy (kJ/kg),Enthalpy (kJ/kg),Entropy (J/g*K),Cv (J/g*K),Cp (J/g*K),Sound Spd. (m/s)	

reference_data=[[-200,0.1,0.656852883,1.522410918,385.3145946,20.69679408,3.117211515,5.196111897,-1280.896984,-1280.899975],
[-150,0.1,0.390447217,2.561165652,645.0435402,23.40265537,3.116502934,5.193727527,-2493.110034,-2493.115033],
[-100,0.1,0.277797953,3.59973855,904.7145858,25.17234509,3.116245295,5.193219801,-3813.850822,-3813.85783],
[-50,0.1,0.215597012,4.638283198,1164.370957,26.48976923,3.116119723,5.193062396,-5210.649368,-5210.658391],
[0,0.1,0.176154738,5.67682716,1424.022429,27.53968004,3.116047811,5.193005697,-6666.12389,-6666.134933],
[50,0.1,0.148911978,6.715376494,1683.672108,28.41259903,3.116002208,5.192985044,-8169.396917,-8169.409982],
[100,0.1,0.128966823,7.753932166,1943.321171,29.15968143,3.115971182,5.19297899,-9713.007171,-9713.022258],
[150,0.1,0.113733376,8.792493792,2202.97011,29.812679,3.115948959,5.192979242,-11291.51439,-11291.5315],
[200,0.1,0.101718424,9.831060675,2462.61914,30.39266086,3.115932402,5.192982251,-12900.77441,-12900.79355],
[250,0.1,0.091999434,10.86963211,2722.268353,30.91432607,3.115919678,5.192986394,-14537.52454,-14537.5457],
[300,0.1,0.083975695,11.90820746,2981.917785,31.38833809,3.115909651,5.192990899,-16199.12894,-16199.15212],
[350,0.1,0.077239246,12.94678619,3241.567443,31.82267898,3.115901584,5.192995392,-17883.41358,-17883.43879],
[400,0.1,0.071503303,13.98536785,3501.217321,32.22347935,3.11589498,5.192999699,-19588.55459,-19588.58182],
[450,0.1,0.066560383,15.02395206,3760.867408,32.59555004,3.115889493,5.193003741,-21312.99981,-21313.02906],
[500,0.1,0.06225666,16.06253851,4020.51769,32.94273575,3.115884875,5.193007495,-23055.4123,-23055.44358],
[-200,1,6.457280697,0.154863951,387.3921081,15.90631526,3.129214702,5.222930511,-931.0188042,-931.0240646],
[-150,1,3.862970962,0.258868112,647.7940785,18.61964275,3.122376607,5.199690343,-1904.083038,-1904.089492],
[-100,1,2.757164953,0.362691394,907.6260328,20.39047233,3.11985039,5.19462848,-2985.675646,-2985.682657],
[-50,1,2.143690747,0.466485197,1167.310908,21.70804868,3.118611405,5.19304006,-4143.325351,-4143.333587],
[0,1,1.753531913,0.570277617,1426.946141,22.75789558,3.117899497,5.192462021,-5359.650652,-5359.660498],
[50,1,1.483514352,0.674075043,1686.562991,23.63070465,3.117447113,5.192248474,-6623.774258,-6623.785916],
[100,1,1.28554763,0.777878607,1946.173402,24.37767593,3.117138901,5.192183542,-7928.23498,-7928.248561],
[150,1,1.134188036,0.881688017,2205.782399,25.03057305,3.116917909,5.192183262,-9267.592605,-9267.60817],
[200,1,1.014710646,0.98550262,2465.392194,25.61046725,3.116753141,5.192211528,-10637.703,-10637.72059],
[250,1,0.91800243,1.089321736,2725.003747,26.13205677,3.116626443,5.192251741,-12035.30349,-12035.32311],
[300,1,0.838121276,1.19314475,2984.617433,26.6060035,3.116526556,5.192295958,-13457.75822,-13457.7799],
[350,1,0.771027184,1.296971133,3244.233344,27.04028792,3.116446163,5.19234032,-14902.8932,-14902.91694],
[400,1,0.71387756,1.400800439,3503.851436,27.4410392,3.116380327,5.192382985,-16368.88454,-16368.91033],
[450,1,0.664614208,1.504632294,3763.4716,27.81306699,3.116325611,5.192423132,-17854.18009,-17854.20795],
[500,1,0.621710227,1.608466382,4023.093702,28.16021501,3.116279553,5.19246047,-19357.44292,-19357.47284],
[-200,10,55.19701306,0.018116922,409.1036874,11.06162602,3.225341125,5.411735808,-581.2234797,-581.1481544],
[-150,10,34.94666506,0.028615034,674.9342659,13.83446066,3.174084797,5.256142982,-1314.929902,-1315.063951],
[-100,10,25.65387249,0.03898047,936.4229805,15.61696025,3.152878523,5.211438829,-2157.458384,-2157.507485],
[-50,10,20.2789361,0.049312252,1196.530449,16.93675776,3.14194186,5.195351306,-3076.029562,-3076.008784],
[0,10,16.76797674,0.059637487,1456.112158,17.98640825,3.135474821,5.188852172,-4053.250129,-4053.186064],
[50,10,14.29298509,0.069964391,1715.477083,18.85837621,3.131288521,5.186141077,-5078.251101,-5078.161851],
[100,10,12.45401215,0.080295409,1974.753323,19.60438768,3.128399398,5.185089766,-6143.578032,-6143.474865],
[150,10,11.03374753,0.090631039,2233.998722,20.25637074,3.126308203,5.184812732,-7243.794944,-7243.684842],
[200,10,9.903823872,0.100971101,2493.240687,20.83544324,3.124737789,5.184907523,-8374.76029,-8374.64763],
[250,10,8.983499148,0.1113152,2752.492376,21.35630955,3.123523392,5.185179429,-9533.212961,-9533.100521],
[300,10,8.219432878,0.121662895,3011.759918,21.82962415,3.122561614,5.18553013,-10716.51812,-10716.40768],
[350,10,7.574967373,0.13201377,3271.045841,22.26335638,3.121784654,5.185908705,-11922.50238,-11922.39509],
[400,10,7.024077523,0.142367449,3530.350804,22.66362414,3.121146413,5.186288589,-13149.34228,-13149.23885],
[450,10,6.54777629,0.152723605,3789.674491,23.03522695,3.120614576,5.186656153,-14395.48593,-14395.38683],
[500,10,6.13188631,0.163081954,4049.016101,23.3819998,3.120165891,5.18700479,-15659.59658,-15659.50209],
[-200,100,233.6687271,0.004279563,666.819624,6.261411852,3.582597454,5.210786236,-219.1589362,-231.2722441],
[-150,100,188.1101751,0.005316034,932.3419323,9.024444808,3.428125288,5.353483468,-710.6217997,-726.0384096],
[-100,100,156.9525687,0.006371352,1199.610895,10.84616765,3.340884222,5.328090065,-1315.538193,-1329.332313],
[-50,100,134.7770079,0.007419663,1464.85012,12.19216734,3.28813505,5.281681364,-1997.798338,-2008.683981],
[0,100,118.229936,0.008458095,1727.906011,13.25597836,3.253453872,5.242263381,-2738.773952,-2746.711629],
[50,100,105.3909673,0.009488479,1989.247225,14.13465307,3.229170247,5.212974197,-3527.213832,-3532.537785],
[100,100,95.11921538,0.010513123,2249.343975,14.88305954,3.211350441,5.192112508,-4355.581981,-4358.701168],
[150,100,86.70077752,0.011533922,2508.562684,15.53499433,3.197796051,5.177519015,-5218.46236,-5219.761515],
[200,100,79.66680028,0.01255228,2767.170799,16.11266156,3.18719061,5.16743675,-6111.763037,-6111.574671],
[250,100,73.69632083,0.013569198,3025.3598,16.6313988,3.178701095,5.160573235,-7032.276327,-7030.877931],
[300,100,68.56181543,0.014585378,3283.266246,17.10223202,3.171776369,5.156005994,-7977.415846,-7975.035455],
[350,100,64.09719794,0.015601306,3540.987669,17.53334908,3.166038194,5.153081193,-8945.049411,-8941.873234],
[400,100,60.17818684,0.016617317,3798.593927,17.93099574,3.161218934,5.151334483,-9933.387544,-9929.567368],
[450,100,56.7097958,0.017633638,4056.135175,18.30004491,3.157124252,5.150434266,-10940.90609,-10936.56572],
[500,100,53.61809533,0.01865042,4313.647453,18.64437195,3.153609871,5.150142163,-11966.29072,-11961.53135]]
[m,n]=shape(reference_data)
absolute=zeros(m)
ideal=zeros(m)
#normalized=absolute
#ideal=absolute
#residual=absolute
h=zeros(m)
Z=zeros(m)
s=zeros(m)
cp=zeros(m)
cv=zeros(m)
speed_o_sound=zeros(m)
P_out=zeros(m)
residual=zeros(m)
M_He=4.00260
R=8.314310#8.314472#8.314472#
R_He=R/M_He
 
for i in range(0,m):
        [Z[i],h[i],s[i],absolute[i],cv[i],cp[i],speed_o_sound[i],ideal[i],residual[i]]=enthalpy_helium(float(reference_data[i][0])+273.15,float(reference_data[i][2]))
        print 1e3*reference_data[i][1]/(reference_data[i][2]*R_He*(reference_data[i][0]+273.15))-Z[i],'\t',reference_data[i][4]-h[i],'\t',s[i]-reference_data[i][5],'\t', cv[i]-reference_data[i][6],'\t',cp[i]-reference_data[i][7],'\t',absolute[i]-reference_data[i][8] 
#print '\n'        
#print 'Z\t\t\t','Enthalpy\t\t','Entropy\t\t','Cv\t\t\t','Cp\t\t\t','Helmholtz\n'         
#for i in range(0,m):
#        print Z[i],'\t',h[i],'\t',s[i],'\t', cv[i],'\t',cp[i],'\t',absolute[i],speed_o_sound[i],ideal[i],residual[i]