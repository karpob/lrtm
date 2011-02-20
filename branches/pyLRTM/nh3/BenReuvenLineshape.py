def BenReuvenLineshape(f,fo,df,ce,pst):
	 """
	 ****Not used directly, could be used to  make HanleySteffes more legible.****
	 This function determines the Ben-Reuven (BR) lineshape over
	 a range of freqencies which are user defined as is the frequency step
	 the resulting vector is passed to this subroutine as f.
	 fo is a vector containing the center frequencies (GHz) of the spectral 
	 lines of the given molecule (ie PH3) as supplied by the "Submillimeter,
	 Millimeter, and Microwave Spectral Catalog" aka Poynter-Pickett Line Catalog
	 but in GHz, not in MHz as provided by the catalog.
	 Actually converted to inverse cm
	 df is elsewhere called d_nu in the main program and is the line half-width at
	 half-maximum or simply linewidth
	 ce is the coupling element  which is calc by function
	 pst is the pressure shift term which is calc by function
	 (BR) Becomes (VVW) if ce=pst=0
	 (BR) Becomes Gross or kinetic if df=ce and pst=0


	 NOTE: All frequencies are or have been converted to inverse cm!!!

	 fbr(f,fo,df,ce,pst)=(2/pi)*(f/fo)^2 * ( (df-ce)f^2 + (df+ce)*[(fo+pst)^2 + df^2 - ce^2])
	                                        -------------------------------------------------
	                                        [f^2 - (fo+pst)^2 - df^2 + ce^2]^2 + 4f^2df^2

	 F =                          A      * (    B        +    C    *[           D          ])
	                                        ------------------------------------------------
	                                        [E  -     J       - G    +  H  ]^2 + I
	"""
	n=f.shape[0]
	m=fo.shape[1]
	# f(1) f(1) f(1) f(1) ....f(n) n times where n is the number of frequency steps
	# f(2) f(2) f(2) f(2)				in the observation range                            
	# ...
	# f(n) f(n) f(n) f(n)   n times where n is the number of frequency steps
	# m times where m is the number of spectral lines

	nones=numpy.ones(n,1)
	mones=numpy.ones(1,m)
	f_matrix=f*mones
	fo_matrix=nones*fo
	df_matrix=nones*df
	ce_matrix=nones*ce
	pst_matrix=nones*pst


	A=(2./numpy.pi)*(numpy.power((f_matrix/fo_matrix),2))			# f(1)*fo(1)  f(2)*fo(1) f(3)*fo(1)
                                          # f(1)*fo(2)  f(2)*fo(2) f(3)*fo(2)
	B=(df_matrix-ce_matrix)*(numpy.power(f_matrix,2))
	C=df_matrix+ce_matrix
	D=(numpy.power((fo_matrix+pst_matrix),2)) + (numpy.power(df_matrix,2))-(numpy.power(ce_matrix,2))
	E=numpy.power(f_matrix,2)
	J=numpy.power((fo_matrix+pst_matrix),2)
	G=numpy.power(df_matrix,2)
	H=numpy.power(ce_matrix,2)
	I=4.*(numpy.power(f_matrix,2))*(numpy.power(df_matrix,2))
	F=A*(B+C*D)/((numpy.power((E-J-G+H),2))+I)						# m x n matrix

	# Maybe not=this under here -been screwing 2/3/01
	# Note: the extra df factor comes in because of the expression for alpha has
	# an additional df factor which cancels with the alpha_max, but is added here
	# for computational reasons.  result is actually df*F
	return F
