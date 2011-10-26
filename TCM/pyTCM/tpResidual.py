def tpResidual(tcme,tcme2,tcme3,tcme4,title_text1,title_text2,plotName):
	import pylab
	from numpy import shape,nonzero
	import numpy
	pylab.figure()
	ax=pylab.subplot(121)
	[m,n]=shape(tcme2)
	#idx=tcme2[:,22]==0.0
	#tcme2[idx,22]=tcme2[idx,0]
	idx1,=numpy.nonzero(tcme[:,0]>8.0)
	idx2,=numpy.nonzero(tcme2[:,0]>8.0)
	idx3,=numpy.nonzero(tcme3[:,0]>8.0)
	idx4,=numpy.nonzero(tcme4[:,0]>8.0)
	pylab.plot(tcme2[idx2,1]-tcme[idx1,1],tcme[idx1,0],'k:')
	pylab.plot(tcme4[idx4,1]-tcme3[idx3,1],tcme3[idx3,0],'k')
	pylab.xlabel('$\Delta$T ($^{\circ}$K)',fontsize=18)	
      #change the yaxis to be log scale
	pylab.gca().set_yscale('log');
	pylab.ylim((8.,200))
        pylab.xlim((-30,5))
	pylab.yticks((10,100),('10','100'),fontsize=18)
	pylab.ylabel('Ideal Pressure (bars)',fontsize=18)
	pylab.xticks(fontsize=18)
	pylab.gca().set_ylim(pylab.gca().get_ylim()[::-1])
	pylab.title(title_text1)
	ax=pylab.subplot(122)
	[m,n]=shape(tcme2)
	#idx=tcme2[:,22]==0.0
	#tcme2[idx,22]=tcme2[idx,0]
	pylab.plot(tcme2[idx2,22]-tcme[idx1,0],tcme[idx1,0],'k:')
	pylab.plot(tcme4[idx4,22]-tcme3[idx3,0],tcme3[idx3,0],'k')
	pylab.xlabel('$\Delta$P (bars)',fontsize=18)	
      #change the yaxis to be log scale
	pylab.gca().set_yscale('log');
	pylab.ylim((8.,200))
        #pylab.xlim((0,800))
	pylab.yticks((),(),fontsize=18)
	#pylab.ylabel('Ideal Pressure (bars)',fontsize=18)
	pylab.xticks(fontsize=18)
	pylab.xlim((-6.,10.))
	pylab.gca().set_ylim(pylab.gca().get_ylim()[::-1])
	pylab.title(title_text2)
	pylab.legend(('Voyager','Galileo'),loc='upper right')
	#pylab.legend(('DeBoer-Steffes','This Work'),loc='upper right')
	
	pylab.savefig(plotName)


