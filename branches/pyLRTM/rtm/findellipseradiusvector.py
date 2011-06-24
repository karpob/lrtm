def findellipseradiusvector(ao,bo,co,major,minor):

	# function ellipses
	# 
	#    This function generates the elliptical shells for each pressure level (dz). This is computed from DeBoer_TCM.m, taking the
	#    values for altitude values taken for the equatoral and polar locations. 

	#    For simplicity, and consistency Karpowicz simply
	#    multiplies the values of dz by the ratio of polar to equatorial radius (at P=1 bar). This is slightly different from 
	#    what Hoffman did, however sensitivity tests concluded that there was no difference when running the TCM twice once with an
	#    Equitorial value of Ro, and another with a polar value of Ro. Changes in gravity could play a role, but are likely to be 
	#    second or third order which don't affect calculation of dz for values on the order of ~1km. At finer altitude steps, this
	#    may show up, and running the TCM twice may be necessary.


	#JPH  function ellipses=findellipseradiusvector(ao,bo,co,major,minor)
	#JPH from master ellipse (largest extent of the planet)
	#JPH finds the ellipsoidal shells for the elliptical-shell model
	#JPH ao,bo,co are the largest extent of each axis (x,y,z)
	#JPH From the TCM.out (which gives dR for each P relative to Ro(P=0.5)
	#JPH After running model (TCM) twice with Ro=major axis, minor axis) gives
	#JPH Vector of ellipsoid vertices, from that figures out the data points in between


	a=ao+major				# semi-minor axis is x
	b=bo+major				# semi-minor axis is also y
	c=co+minor				# semi-major axis is z
	ellipses={}
	ellipses['a']=a
	ellipses['b']=b
	ellipses['c']=c
	return ellipses
