import TCM
dz =1
XHe_i =0.1351
XH2S_i =6.6254e-05
XNH3_i =3.9150e-04
XH2O_i =0.0022
XCH4_i =0.0018
XPH3_i =5.1627e-07
XCO =0
P_temp =1000
T_temp =1.5838e+03
g0_i =2330
R0_i =7.1492e+09
P0_i =1
T_targ_i =166
P_targ_i =1
P_term_i =0.141
use_lindal_in=1
n_lindal_pts_i =61
SuperSatSelf1_i=0.
SuperSatSelf2_i=0.
SuperSatSelf3_i=0.
SuperSatSelf4_i=0.
supersatNH3_i=0.
supersatH2S_i=0.
AutoStep_constant=8
fp =666
dP_init =10
dP_fine =0.2500
P_fine_start=13
P_fine_stop=1
use_dz =0
frain =3
select_ackerman=0

a=TCM.intoTheVoid( dz,
		  XHe_i,
		  XH2S_i,
                  XNH3_i,
                  XH2O_i,
                  XCH4_i,
                  XPH3_i,
                  XCO,
                  P_temp,
                  T_temp,
                  g0_i,
                  R0_i,
                  P0_i,
                  T_targ_i,
                  P_targ_i,
                  P_term_i,
                  use_lindal_in,
                  n_lindal_pts_i,
                  SuperSatSelf1_i,                                             
                  SuperSatSelf2_i,
                  SuperSatSelf3_i,
                  SuperSatSelf4_i,
                  supersatNH3_i,
                  supersatH2S_i,
                  AutoStep_constant,
                  fp,
                  dP_init,
                  dP_fine,
                  P_fine_start,
                  P_fine_stop,
                  use_dz,
                  frain,
                  select_ackerman)

print a
#print TCM.getFloatValues(0,0)
#TCM.getCloudFlags(a)

b=[]
c=[]
#print a
#for i in range(0, a):
#	print i
#	b.append(TCM.getCldVal(i))
for i in range(0, a):
	cc=[]
	for j in range(0,20):
		print TCM.getFloatValues(j,i)
		cc.append(TCM.getFloatValues(j,i))
	c.append(cc)
import numpy
print c[0][0]
a=numpy.asarray(c)
print a[:,0],a[0:,2]


