def rayellipseint(Rayorigin,Raydirection,ellipse):

	# function [intercept,internormal,d]=rayellipseint(Rayorigin,Raydirection,ellipses)
	# Haines Ray Intercept Test (with Ellipse modification)
	# Gets passed ray-info and ellipse info for THE ellipse currently being done
	# IMPOSE
	# ELLIPSE MODIFICATION: Sr=1, elliptical equation for spherical, origin=0,0,0
	#
	# With inputs of ray location, ray vector, center of sphere, and sphere radius
	# This finds the intercept location (in cartesian), and the local normal, as well
	# as the distance (d) to the sphere from the observer (ray origin)

	# Can be used to find the original distance to the planet from the craft
	# along with the intercept point and distance
	# Also can find the pathlength through the spherical shells model
	# In this case the ray origin is the last intercept point and the ray direction
	# must be found a priori from snells law for spherical shells and refractive index
	# 
	# Modification by Karpowicz. Include "negation" when ray origin is within the
	# ellipsoid to a tolerance of 1e-6. (if the equation for an ellipse is <1 to
	# a tolerance of 1e-6 or 0.9999....) 
	#fred=1

	# Put into notation of the Haines text (Ray tracing, Glassner ed)
	X0=Rayorigin[0]
	Y0=Rayorigin[1]
	Z0=Rayorigin[2]

	Xd_notunit=Raydirection[0]		# Not normalized, need to
	Yd_notunit=Raydirection[1]
	Zd_notunit=Raydirection[2]


	Xc=0.						# ASSUMING Centered
	Yc=0.
	Zc=0.

	Sr=1						# For ellipse always equal to unity

	a=1./(ellipse['a']**2)	# ellipse x^2/a^2 + y^2/b^2 etc=> making a==1/a^2 instead (faster)
	b=1./(ellipse['b']**2)
	c=1./(ellipse['c']**2)

	#added by BMK
	#
	# Check to see if ray originates from within ellipsoid
	# as per Mohammed, 2005

	origin_ellipse_val=((X0**2)*a)+((Y0**2)*b)+((Z0**2)*c)
	tolerance_roundoff=1e-6
	flip_rnormal=0.
	if(origin_ellipse_val<1.-tolerance_roundoff):
    		flip_rnormal=1.
   		print origin_ellipse_val

	####
	# End conversion of inputs
	####	
	
	# Normalize Ray direction vector
	Raydirection_length=1        

	Xd=Xd_notunit/Raydirection_length
	Yd=Yd_notunit/Raydirection_length
	Zd=Zd_notunit/Raydirection_length
	# Normalized

	# Using solution in Haines-see text pages 35-37 for details
	# Using same notation
	#A=Xd**2+Yd**2+Zd**2					# Should equal '1' already (good error check)
	#B=2*(Xd*(X0-Xc)+Yd*(Y0-Yc)+Zd*(Z0-Zc))
	#C=(X0-Xc)**2+(Y0-Yc)**2+(Z0-Zc)**2-Sr**2

	# Modification for Ellipsoid
	A=a*(Xd**2)+b*(Yd**2)+c*(Zd**2)					# Should equal '1' already (good error check)
	B=2*(a*Xd*X0+b*Yd*Y0+c*Zd*Z0)
	C=a*(X0**2)+b*(Y0**2)+c*(Z0**2)-Sr**2

	# Quadratic solutions
	# If discriminant is negative (so its root is imaginary) ray misses sphere
	# A Check and precalculation for later
	limbflag=0				# Fo checking ray goes outside of atmopshere
	squarerootterm=B**2-4.*A*C
	if squarerootterm < 0
   		squarerootpart=numpy.sqrt(squarerootterm)
		to=(-B-squarerootpart)/(2.*A)		# debug values(-B-sqrt(B^2-4*C))*(0.5)
   		t1=(-B+squarerootpart)/(2.*A)		# debug values(-B+sqrt(B^2-4*C))*(0.5)
   		limbflag=1								# This ray does not hit planet or skips at 90 deg
   		d=0
   		#disp('Ray never intersects sphere')
   		intercept=numpy.array([0., 0., 0.])
   		internormal=numpy.array([0., 0., 0.])
   		return intercept,internormal,d,limbflag
	end

	squarerootpart=numpy.sqrt(squarerootterm)
	to=(-B-squarerootpart)/(2.*A)		# (-B-sqrt(B^2-4*C))*(0.5)
	t1=(-B+squarerootpart)/(2.*A)		# (-B+sqrt(B^2-4*C))*(0.5)

	# Is there a solution????
	if to<0.:
   		if t1<0.:
      			limbflag=1
      			intercept=numpy.array([0., 0., 0.])
      			internormal=numpy.array([0., 0., 0])
      			d=0.
      			print 'Both to, t1 less than zero'
      			return intercept,internormal,d,limbflag

	# Seems to be a solution


	# Must find the smallest positive real root  (previous check should have found im roots)
	# Check for smaller-then check if thats positive-if not check if other is positive
	# If so, thats the answer, if not error that should have been caught already

	# Easier way to check, if B is greater than squareroot then to is negative, only
	# check t1
	if to<t1:				# Which is lesser?									
   		if to>0							# here to is lesser, but is it positive?
      			d=to							# here to is also positive, so its the solution
   		elif t1>0:						# here to is negative so is t1 positive?
      			d=t1							# here t1 is positive, so its the solution
   		else:
      			print 'No solution-Problem'
   			sys.exit()
	elif to>t1:						# t1 was lesser magnitude
   		if t1>0:							# is t1 positive?
      			d=t1							# here t1 is positive so its the solution
   		elif to>0:						# here t1 is negative so is to positive?
      			d=to							# here to is positive, so its the solution
   		else:
      			print 'No solution2-Problem'
			sys.exit()
	elif t1==to:
   		if t1>0:
      			d=t1
      			limbflag=1		# ray skips along tangent to surface, so dz is up not down


	# To keep with Haines notation t=d (t=either t1,to) solved above
	# d is same as t

	# Find the intersection point
	xi=X0+Xd*d
	yi=Y0+Yd*d
	zi=Z0+Zd*d

	# Find the normal to the intersection point
	normalize=numpy.sqrt(ellipse.a^2+ellipse.b^2+ellipse.c^2)
	unita=ellipse.a/normalize
	unitb=ellipse.b/normalize
	unitc=ellipse.c/normalize
	# recall a,b,c===1/(ellipse.a^2) etc
	delnormal_num=[(xi/(ellipse['a']**2)) (yi/(ellipse['b']**2)) (zi/(ellipse['c']**2))]
	delnormal_den=numpy.sqrt( ((xi**2)/(ellipse['a']**4))+((yi**2)/(ellipse['b']**4))+((zi**2)/(ellipse['c']**4)) )
	xn=delnormal_num[0]/delnormal_den
	yn=delnormal_num[1]/delnormal_den
	zn=delnormal_num[2]/delnormal_den


	intercept=numpy.array([xi, yi, zi])			# This is a vector location
	internormal=numpy.array([xn, yn, zn])			# This is a vector direction
	#added by BMK to ensure conditions according to Mohammed
	if(flip_rnormal==1):
    		internormal=-1.*internormal
	isone=(xi**2)*a+(zi**2)*c+(yi**2)*b
	print 'is one?',isone
	return intercept,internormal,d,limbflag
