def fwght(layers,masterindex):

	# function fwght
	#
	# This function checks to see that double contributions are included in the weighting functions (intended for multiple rays).
	# Right now we only worry about boresight weighting function.
        #
        import numpy
        
        masterindex=numpy.asarray(masterindex).T
	sl=layers.shape[0]
	wght=numpy.zeros([masterindex.shape[0]])
	for n in range(0,1):
   		for k in range(0,masterindex.shape[0]):
      			ch=masterindex[n]==masterindex[k]	# are there double contributions?
      			#matchindextemp=masterindex(:).*ch  # Zero out other layer indices
      			matchlayers=ch*layers # matched layer values
      			wght[k]=matchlayers.sum(axis=0)						# add those layers together 

	return wght
