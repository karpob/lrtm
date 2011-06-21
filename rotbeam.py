def rotbeam(Raydirection,b1): 
	import numpy
	import scipy.linalg
	x=numpy.array([1., 0., 0.])
	z=numpy.array([0., 0., 1.])

	c=Raydirection

	rvect_e=numpy.cross(z,c)
	rvect_e=rvect_e/numpy.norm(rvect_e)
	rang_e=arccos(numpy.dot(z,c))
	kx=rvect_e[0]
	ky=rvect_e[1]
	kz=rvect_e[2]
	v=numpy.norm(numpy.array([ky, kz]))
	Rx=numpy.array([[1., 0., 0.],
                        [0., kz/v, ky/v],
                        [0., -1.0*ky/v, kz/v]])
	Ry=numpy.array([[v, 0., kx],
   	               [0., 1., 0.],
                       [-kx 0. v]])
	Rz=numpy.array([[numpy.cos(rang_e), -1.0*numpy.sin(rang_e), 0.]
   	               [numpy.sin(rang_e), numpy.cos(rang_e), 0.],
  	               [0., 0., 1.]])

	# To transform any vector to align with the unit vector c
	# Multiply that vector with this matrix
	#
	#  new_vector=Rx*Ry*Rz*inv(Ry)*inv(Rx)*old_vector
	#
	#
	# Example of what I think you're trying to do


	Zr=(Rx*Ry*Rz*linalg.inv(Ry)*linalg.inv(Rx)*numpy.transpose(z))
	Vr1=(Rx*Ry*Rz*linalg.inv(Ry)*linalg.inv(Rx)*b1)

	return [Vr1,Zr]
