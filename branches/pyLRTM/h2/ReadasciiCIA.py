def ReadasciiCIA(name):
	"""
	Reads Glenn Orton's CIA files. Add pickle dump
	to make this faster.
	"""
	import re
	import numpy
	
	file_handle=open(name,'r')
	line1=file_handle.readline()
	#use regular expressions to tame this file
	nTemps,maxTemp,minTemp=re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", line1)
	nTemps=int(nTemps)
	maxTemp=float(maxTemp)
	minTemp=float(minTemp)	
	line2=file_handle.readline()
	nFreq=re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", line2)
	nFreq=int(nFreq[0])
	Values=file_handle.read()
	vals=re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", Values)

	Wavenumbers=[]
	for i in range(0,nFreq):
		Wavenumbers.append(vals.pop(0))
	Wavenumbers=numpy.asarray(Wavenumbers,dtype=float)

	#equilibrium H2-H2 collisions
	alpha_eH2=[]	
	for i in range(0,nFreq*nTemps):
		alpha_eH2.append(vals.pop(0))
	alpha_eH2=numpy.asarray(alpha_eH2,dtype=float)
	alpha_eH2=alpha_eH2.reshape(nTemps,nFreq)

	#normal H2-H2 collisions

	alpha_nH2=[]
        for i in range(0,nFreq*nTemps):
                alpha_nH2.append(vals.pop(0))
        alpha_nH2=numpy.asarray(alpha_nH2,dtype=float)
        alpha_nH2=alpha_nH2.reshape(nTemps,nFreq)	
	
	#equilibrium H2-He collisions
	alpha_eH2_He=[]
        for i in range(0,nFreq*nTemps):
                alpha_eH2_He.append(vals.pop(0))
        alpha_eH2_He=numpy.asarray(alpha_eH2_He,dtype=float)
        alpha_eH2_He=alpha_eH2_He.reshape(nTemps,nFreq)
	
	#equilibrium nH2-He collisions
	alpha_nH2_He=[]
        for i in range(0,nFreq*nTemps):
                alpha_nH2_He.append(vals.pop(0))
        alpha_nH2_He=numpy.asarray(alpha_nH2_He,dtype=float)
        alpha_nH2_He=alpha_nH2_He.reshape(nTemps,nFreq)

	#equilibrium H2-CH4 collisions
        alpha_eH2_CH4=[]
        for i in range(0,nFreq*nTemps):
                alpha_eH2_CH4.append(vals.pop(0))
        alpha_eH2_CH4=numpy.asarray(alpha_eH2_CH4,dtype=float)
        alpha_eH2_CH4=alpha_eH2_CH4.reshape(nTemps,nFreq)
	
	#save time...don't use, dont parse.
	#normal H2-CH4 collisions
        #alpha_nH2_CH4=[]
        #for i in range(0,nFreq*nTemps):
        #        alpha_nH2_CH4.append(vals.pop(0))
        #alpha_nH2_CH4=numpy.asarray(alpha_nH2_CH4,dtype=float)
        #alpha_nH2_CH4=alpha_nH2_CH4.reshape(nTemps,nFreq)	
	
	return 	nTemps,minTemp,maxTemp,Wavenumbers,alpha_eH2,alpha_eH2_He,alpha_eH2_CH4
