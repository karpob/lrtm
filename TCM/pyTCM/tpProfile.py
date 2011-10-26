def tpProfile(tcme,tcme2,tcme3,tcme4,title_text,plotName):
	import pylab
	from numpy import shape
	pylab.figure()
	[m,n]=shape(tcme2)
	#idx=tcme2[:,22]==0.0
	#tcme2[idx,22]=tcme2[idx,0]
	pylab.plot(tcme[:,1],tcme[:,0],'r')
	pylab.plot(tcme2[:,1],tcme2[:,22],'k')
	pylab.plot(tcme3[:,1],tcme3[:,0],'r:')
	pylab.plot(tcme4[:,1],tcme4[:,22],'k:')
	pylab.xlabel('Temperature ($^{\circ}$K)',fontsize=18)	
      #change the yaxis to be log scale
	pylab.gca().set_yscale('log');
	pylab.ylim((0.1,200))
        pylab.xlim((0,800))
	pylab.yticks((0.1,1,10,100),('0.1','1','10','100'),fontsize=18)
	pylab.ylabel('Pressure (bars)',fontsize=18)
	pylab.xticks(fontsize=18)
	pylab.gca().set_ylim(pylab.gca().get_ylim()[::-1])
	pylab.legend(('DeBoer, 1995 w/ Galileo Profile','This Work w/ Galileo Profile','DeBoer, 1995 w/ Voyager Profile','This Work w/ Voyager Profile'),loc='upper right')
	pylab.title(title_text)
	pylab.savefig(plotName)
	


