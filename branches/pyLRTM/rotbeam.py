def rotbeam(Raydirection,b1):
        """
        Take sampled antenna pattern and align it with the unit vector givng the antenna pointing direction.
        
                Input:
                        -->Raydirection: Unit vector of the antennas pointing 
                        -->b1: collection of rays which sample antenna pattern
                Output:
                        <-Vr1: b1 values aligned with Ray direction
                        <-Zr1: positive z rotated to Raydirection (debugging).                
        """ 
	import numpy
	import scipy.linalg
	x=numpy.array([1., 0., 0.])
	z=numpy.array([0., 0., 1.])

	c=Raydirection

	rvect_e=numpy.cross(z,c)
	rvect_e=rvect_e/numpy.linalg.norm(rvect_e)
	rang_e=numpy.arccos(numpy.dot(z,c))
	kx=rvect_e[0]
	ky=rvect_e[1]
	kz=rvect_e[2]
	v=numpy.linalg.norm(numpy.array([ky, kz]))
	Rx=numpy.array([[1., 0., 0.],
                        [0., kz/v, ky/v],
                        [0., -1.0*ky/v, kz/v]])
	Ry=numpy.array([[v, 0., kx],
   	               [0., 1., 0.],
                       [-kx, 0., v]])
	Rz=numpy.array([[numpy.cos(rang_e), -1.0*numpy.sin(rang_e), 0.],
   	               [numpy.sin(rang_e), numpy.cos(rang_e), 0.],
  	               [0., 0., 1.]])

	# To transform any vector to align with the unit vector c
	# Multiply that vector with this matrix
	#
	#  new_vector=Rx*Ry*Rz*inv(Ry)*inv(Rx)*old_vector
	
	Zr=numpy.asmatrix(Rx)*numpy.asmatrix(Ry)*numpy.asmatrix(Rz)*scipy.linalg.inv(numpy.asmatrix(Ry))*scipy.linalg.inv(numpy.asmatrix(Rx))*(numpy.asmatrix(z).T)
	Vr1=numpy.asmatrix(Rx)*numpy.asmatrix(Ry)*numpy.asmatrix(Rz)*scipy.linalg.inv(numpy.asmatrix(Ry))*scipy.linalg.inv(numpy.asmatrix(Rx))*numpy.asmatrix(b1)
        Zr=numpy.asarray(Zr)
        Vr1=numpy.asarray(Vr1)
        	return [Vr1,Zr]
