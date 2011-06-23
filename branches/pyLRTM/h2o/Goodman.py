def Goodman(f,P,T,XH2O,XH2,XHe):
	"""
	Goodman,1969 Model. for water vapor absorption
	Input:
	        ->f: frequency in GHz
	        ->P: Pressure in bars.
	        ->T: Temperature in Kelvin.
	        ->XH2O mole fraction water vapor
	        ->XH2 mole fraction hydrogen
	        ->XHe mole fraction helium
	Output:
	        <-alpha: absorption coefficient in dB/km        
	"""
        from numpy import power
        OpticaldepthstodB=434294.5
        nu=((f*1e9)/2.99792458e8)/100;
        PH2O=XH2O*P*750.06;
        P=P*750.06;

        #stuff from equation B13
        front_part=PH2O*(power((273.0/T),(13.0/3.0)))*power(nu,2);

        delta_nu_1=0.1*((P/760.0)*power((273.0/T),(2.0/3.0)))*(0.810*XH2+0.35*XHe);

        frac_minus_074=(delta_nu_1)/(power((nu-0.74),2) +power(delta_nu_1,2));

        frac_plus_074=(delta_nu_1)/(power((nu+0.74),2.0) +power(delta_nu_1,2));

        middle_part=(1.073e-8)*(frac_minus_074+frac_plus_074);

        end_part=(17.20e-8)*delta_nu_1;

        alpha=front_part*(middle_part+end_part)*OpticaldepthstodB
        return alpha

