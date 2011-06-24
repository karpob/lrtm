def DeBoerCorrected(f_in,T,P_H2,P_He,P_H2O):
        """
         function alpha_watervapor=fwatervaporalpha(P_H2O,P_H2,P_He,f,T)
         Input:
                ->f: Frequency (GHz)
                ->T: Temperature (K)
                ->P_H2: Partial Pressure Hydrogen (bars)
                ->P_He: Partial Pressure Helium (bars)
                ->P_H2O: Partial Pressure Water vapor (bars)
         Finds the attenuation due to water vapor
         Taken from DeBoer Thesis page 66 (and related)
         Modified version of Goodman (1969) with Ulaby (1981), Joiner Steffes (1991) etc
         Assuming units are bars and GHz, may need to change
         Check if Bars/GHz correct units
         Uses function watervapor_broadening
         
         Output:
                <-alpha_watervapor: absorption coefficient in dB/km
        """
        from numpy import asarray,power,exp,ones,shape,pi
        from broadening import broadening
        
        # For this calc, f should be a column vector [1
        #                                             2...]
        # since fo is a row vector[1,2,...]

        # Want to get an alpha for each 'f'
        # So alpha(f)
        #if (P_H2O<=0):
        #        alpha_watervapor=0;
        #        return alpha_watervapor

        
        #  print shape(P_H2),shape(P_He),shape(P_H2O),shape(T),shape(f)
        # Equation for alpha is sum of contributions to observation freq from each
        # line frequency, with contributions from each broadening agent
        # Broadening agents are mainly H2, He, and H2O

        To=300.0			# Reference Temperature 300 Kelvin
        Tdiv=To/T			# Ratio of Ref Temp to actual temp, used often


        fo=asarray([22.23515, 183.31012, 323, 325.1538, 380.1968, 390, 436, 438, 442, 448.0008])
        Ep=asarray([644, 196, 1850, 454, 306, 2199, 1507, 1070, 1507, 412])
        A=asarray([1.0, 41.9, 334.4, 115.7, 651.8, 127.0, 191.4, 697.6, 590.2, 973.1])
        g_H2=asarray([2.395, 2.4000, 2.395, 2.395, 2.390, 2.395, 2.395, 2.395, 2.395, 2.395])
        g_He=asarray([0.67, 0.71, 0.67, 0.67, 0.63, 0.67, 0.67, 0.67, 0.67, 0.67])
        g_H2O=asarray([10.67, 11.64, 9.59, 11.99, 12.42, 9.16, 6.32, 8.34, 6.52, 11.57])

        zeta_H2=asarray([0.900, 0.950, 0.900, 0.900, 0.850, 0.900, 0.900, 0.900, 0.900, 0.900])
        zeta_He=asarray([0.515, 0.490, 0.515, 0.490, 0.540, 0.515, 0.515, 0.515, 0.515, 0.515])
        zeta_H2O=asarray([0.626, 0.649, 0.420, 0.619, 0.630, 0.330, 0.290, 0.360, 0.332, 0.510])

        # Calls fwatervaporbroadening to get gammawater
        gammawater=broadening(P_H2O,P_H2,P_He,T)
        n=len(f_in)  #returns the number of columns in f
        f_ones=ones([n,1])#Vector to get 2d arrays for multiplication
        f=f_in
        			# Vector of ones Length of column vector 'f'
        # Generate fo matrix
        fo_m=f_ones*fo
        				#Matrix made up of repeated rows of fo- (size(f) columns)
        gammawater_m=f_ones*gammawater

        
        # Do same thing for 'f'
        fo_ones=ones([1,len(fo)])			# Vector of ones Length of row vector 'fo'
        
        # Generate fo matrix
        f_m=fo_ones*f			        # Matrix made up of repeated columns of f- (size(fo) rows)
        

        # alpha_H2O=5.34*10^5*PH2O*(300/T)^(3/2)*sum[i(1:10)*A(i)*exp(-Ep(i)/T)*(Shape)]+ffactor
        # A*exp part
        OpticaldepthstodB=434294.5 				#/* converts from cm^-1 to dB/km */
        front=(1134.5)*(P_H2O)*(power((300/T),(7.0/2.0)))
        
        Aexp_part=A*exp(-Ep/T);
        Aexp_part=f_ones*Aexp_part;
       
        shape_part=((4.0*power(f_m,2))*gammawater_m/pi)/(power((power(fo_m,2)-power(f_m,2.0)),2.0)+((4.0*power(f_m,2.0))*power(gammawater_m,2)));
       
        second_part=((3.39e-3)*P_H2O)*(power((300.0/T),(3.1)))*(0.81*P_H2+0.35*P_He)*power(f,2.0)

        prod_part=Aexp_part*shape_part
        sum_part=asarray([prod_part.sum(axis=1)]).T
        
        alpha_watervapor=(front*sum_part+second_part)#/OpticaldepthstodB
        
        return alpha_watervapor
