def brParameters(To,T,PH2,PHe,PH2S):
	"""
	To is reference Temperature-usually 300 k
	 T is actual temperature
	 Mol's is for extending functionality to use more than just H2,He,PH3
	 dP is the partial pressure of each molecule contribution

	 Calculates the linewidth for a gas mixture of mixtures (i) and produces
	 a delta_nu line broadening in GHz for each line (j)
	 delta_nu will be a vector of length (j)
	 delta_nu(j)=sum{i}(delta_nu(i,j)*dP(i)*(To/T)^-(epsilson)
	 where epsilon is (m+1)/2(m-1)   m is related to force law 
	 m=3 dipole, m=infi is hardball

	 CHANGE: ORDER MUST BE H2,He,H2S

	 LOAD in CORRECT PRESSURE BROADENING FILES	
	"""
	from scipy.io import loadmat
	# Spilker Ammonia
	pow1=0.7
	Tdiv=To/T


	LineParameters=loadmat('h2s/h2sLineParameters/h2slin.mat')				# H2S line catalog
	gH2=1.96
	gHe=1.20
	gH2S=LineParameters['gH2So']

	gammaH2S=((gH2*PH2+gHe*PHe+gH2S*PH2S)*(Tdiv)**pow1)
	psiH2S=gammaH2S

	return gammaH2S,psiH2S
