#Data from Seward and Franck pXT data


Data_pXT_T=asarray([374.3,   
           374.5,   
           374.5,   
           374.5,   
           374.5,   
           374.5,   
           374.5,   
           375.3,   
           376.5,   
           379.0,   
           381.3])+273.15 #Celcius-->Kelvin
Data_pXT_P=asarray([229,
            237, 
            270,
            320,
            372,
            428, 
            690, 
            1010, 
            1410,
            2020,
            2520])*Bars_to_kPa # Bars->kPa
Data_pXT_x_h2=asarray([0.005,
               0.01,
               0.03,
               0.06,
               0.09,
               0.12,
               0.20,
               0.25,
               0.30,
               0.35,
               0.38])
Data_pXT_x_h2o=1.0-Data_pXT_x_h2                

Data_pXT_Mmix=((Data_pXT_x_h2o*h2o_M_amu)+(Data_pXT_x_h2*h2_M_amu))
Data_pXT_Vol=asarray([56.88,
              56.67,
              55.80,
              54.51,
              53.34,
              52.14,
              43.50,
              37.41,
              33.60,
              30.91,
              29.66])*(1.0/Data_pXT_Mmix)*(1.0/1.0e6)*1000.0 # cm**3/mol ->m**3/kg
                      #(cm**3/mol)(1 m**3/1e6 cm**3)(1 mol/Mmix(g))(1000 g/ 1 kg)
      

Data_pXT_density=1.0/Data_pXT_Vol #kg/m**3
              
