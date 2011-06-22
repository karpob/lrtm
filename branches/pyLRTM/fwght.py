def fwght(layers,masterindex):

	# function fwght
	#
	# This function checks to see that double contributions are included in the weighting functions
	#


	sl=layers.shape[1]
	for n in range(0,sl):
   		for k in range(0,len(masterindex[:,n])):
      			ch=masterindex[:,n]==masterindex[k,n]	# are there double contributions?
      			#matchindextemp=masterindex(:).*ch  # Zero out other layer indices
      			matchlayers=ch*layers[:,n] # matched layer values
      			wght[k,n]=matchlayers.sum(axis=0)						# add those layers together 

	return wght
