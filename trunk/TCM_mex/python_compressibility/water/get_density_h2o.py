def get_density_h2o(Ph2o,T):
        from scipy.optimize import fsolve
        #Take a given measured water pressure, return density
        # in kg/m^3
        R=0.46151805;#Specific gas constant for water
        Ph2o=100*Ph2o;
        #solve for the root with initial guess being the ideal behavior       
	density=fsolve(res_press,(Ph2o/(R*T)),args=(Ph2o,T))

def res_press(x,Ph2o,T):
        #minimization function residual between given pressure and calculated pressure with value x
	from IAPWS95 import pressure
        return Ph2o-pressure(x,T)
