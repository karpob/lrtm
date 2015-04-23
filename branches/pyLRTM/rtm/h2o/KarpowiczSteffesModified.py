def KarpowiczSteffes(f_in,density_h2o,density_h2,density_he,T):
	"""
	f_in in GHz
	density_X mass density of h2o/h2/he
	T in Kelvin
	"""

# This function also requires the vvwlineshape function originally written by Jim Hoffman (with removal of df factor).
# Added shift terms from Rosenkranz.
        from numpy import asarray,ones,power,exp,shape
        from  vvwlineshape_modified import vvwlineshape_modified
        mbars_to_bars=0.001
        inv_km_to_dB=4.342945
        convert_to_km=1e-4
        To=300.0
        Theta=300.0/T
        n=len(f_in)  #returns the number of columns in f
        nones=ones([1,n])#Vector to get 2d arrays for multiplication
        
        NA=6.0221415e23#molecules/mol
        M_amu=float(8.314472/0.46151805)#g/mol
        isotope_partition=0.997317
        M_amu_he=4.0026
        M_amu_h2=2.01594
        f=f_in*nones
        P_h2o=(1.0/M_amu)*density_h2o*8.314472e-5*T
        P_He=(1.0/M_amu_he)*density_he*8.314310e-5*T
        P_H2=(1.0/M_amu_h2)*density_h2*8.314472e-5*T
        
        density_h2o=isotope_partition*(density_h2o/M_amu)*NA*(1.0/1e6) #need in molecules/cc
        
        
        # Center Frequencies in GHZ
        f_o=asarray([[22.2351, 183.3101, 321.2256, 325.1529, 380.1974, 439.1508, 443.0183, 448.0011, 470.8890, 474.6891, 488.4911, 556.9360, 620.7008, 752.0332, 916.1712]])

        # Line intensities
        I_o=asarray([[.1314E-13, .2279E-11, .8058E-13, .2701E-11, .2444E-10,.2185E-11, .4637E-12, .2568E-10, .8392E-12, .3272E-11, .6676E-12, .1535E-08, .1711E-10, .1014E-08, .4238E-10]])
        # Temperature coefficients
        E_o=asarray([[2.144, .668, 6.179, 1.541, 1.048, 3.595, 5.048, 1.405,3.597, 2.379, 2.852, .159, 2.391, .396, 1.441]])


        #self broadening parameters converted to bars
        w_s=1.0/mbars_to_bars*asarray([[.01349, .01466, .01057,.01381,.01454,  .009715,.00788, .01275, .00983, .01095,.01313, .01405, .011836,.01253, .01275]])
        x_s=asarray([[0.61,.85, .54, .74, .89,.62, .50,  .67, .65, .64, .72,1.0, .68, .84,.78]])
        
        #shift terms
        SR=asarray([[ 0, 0, 0, 0, 0, 0.0,0.0,0.0,0.0,0.0,0.0, 0, 0.0,0.0,0.0]])
        #Measured for 183, and 380 only
        w_H2=asarray([2.395, 2.4000, 2.395, 2.395, 2.390, 2.395, 2.395, 2.395, 2.395, 2.395,2.395,2.395,2.395,2.395,2.395])
        w_He=asarray([0.67, 0.71, 0.67,0.67, 0.63, 0.67, 0.67, 0.67, 0.67,0.67,0.67, 0.67, 0.67,0.67, 0.67])

        x_H2=asarray([0.900, 0.950, 0.900, 0.900, 0.850, 0.900, 0.900, 0.900, 0.900, 0.900,0.900,0.900,0.900,0.900,0.900])
        x_He=asarray([0.515, 0.490, 0.515, 0.490, 0.540, 0.515, 0.515,0.515, 0.515,0.515,0.515,0.515,0.515,0.515,0.515])
        
        expo=E_o*(1-Theta)
        S=I_o*power(Theta,2.5)*exp(expo)

        alpha_noshape_mat=S.T*nones

        #df aka gamma delta-nu change aka pressure broadening term
        df=w_s*P_h2o*power(Theta,x_s)+w_H2*P_H2*power(Theta,x_H2)+w_He*P_He*power(Theta,x_He)#w_f*P_f*power(Theta,x_f)#+SR*w_f*P_f*power(Theta,x_f)

        # calculate van-vleck sum it up...
        F = vvwlineshape_modified(f,f_o.T,df.T,SR.T)
        vvw_alpha_matrix=alpha_noshape_mat*F
        line_contribution=inv_km_to_dB*convert_to_km*density_h2o*vvw_alpha_matrix.sum(axis=0)

        #Continuum Terms Foreign, and self 
        Cf_He=((1.0/mbars_to_bars)**2)*1.03562010226e-10#  (dB/km)/((GHz x kPa)^2)->db/km((GHz bars)^2)  #eqn 10 Rosenkranz,1998
        Cf_H2=((1.0/mbars_to_bars)**2) *5.07722009423e-11  #
        Cs1=  3.1e-7*power(Theta,12)#   #equation 5 (correction applied from ^+6,^-6) Rosenkranz,1998
        Cs2=  2.10003048186e-26 *((1.0/mbars_to_bars)*P_h2o)**(6.76418487001)*power(Theta,0.0435525417274)  
        P_f=P_He+P_H2
        Foreign_Continuum_He=Cf_He*P_He*P_h2o*power(f,2)*power(Theta,3) #foreign continuum term from Eqn 6, Rosenkranz
        Foreign_Continuum_H2=Cf_H2*P_H2*P_h2o*power(f,2)*power(Theta,3)
        Foreign_Continuum=Foreign_Continuum_He+Foreign_Continuum_H2
        
        Self_Continuum=Cs1*(P_h2o*(1.0/mbars_to_bars))**(2.0)*power(f,2.0)#+Cs2*power(f,2.0)#*power(Theta,3) #self continuum term from Eqn 6, Rosenkranz
          
        alpha_h2o=line_contribution+inv_km_to_dB*Foreign_Continuum+inv_km_to_dB*Self_Continuum
        return alpha_h2o
