def constituentPlot(tcme_1,tcme_2,title_text,plotName):
	import pylab
	from matplotlib.font_manager import fontManager,FontProperties
	# Just in case user decides to use LaTeX in title_text
	pylab.figure()	
 	pylab.rc('text',usetex=True)
	#set up the contour axes for easy adjustment
	pylab.plot(tcme_1[:,5]*1e6,tcme_1[:,0],'r',linewidth=1) #h2s
        pylab.plot(tcme_1[:,6]*1e6,tcme_1[:,0],'k',linewidth=1)#nh3
        pylab.plot(tcme_1[:,7]*1e6,tcme_1[:,0],'b',linewidth=1)	#h2o
	pylab.plot(tcme_1[:,8]*1e6,tcme_1[:,0],'g',linewidth=1)#ch4
	pylab.plot(tcme_2[:,5]*1e6,tcme_2[:,0],'r:',linewidth=1) #h2s
        pylab.plot(tcme_2[:,6]*1e6,tcme_2[:,0],'k:',linewidth=1)#nh3
        pylab.plot(tcme_2[:,7]*1e6,tcme_2[:,0],'b:',linewidth=1)	#h2o
	pylab.plot(tcme_2[:,8]*1e6,tcme_2[:,0],'g:',linewidth=1)#ch4
	#pylab.plot(tcme_1[::-1,9]*1e6,tcme_1[:,0],'m')
	

        #change the yaxis to be log scale
	pylab.gca().set_yscale('log');
	# change the xaxis to be log scale
	pylab.gca().set_xscale('log');

	#Set the tick values
	# 10^(-1) 10^0 10^1 10^2 10^3 == Ugly
	# instead make the axis read 0.1, 1,10,100,1000
	pylab.ylim((0.1,1000))
	pylab.xlim((1e-3,1e6))
	
	pylab.yticks((0.1,1,10,100,1000),('0.1','1','10','100','1000'),fontsize=18)

	#reverse the order of the y axis (change so pressure goes up as you go down in y)
	pylab.gca().set_ylim(pylab.gca().get_ylim()[::-1])

	#Label x and y axis
	pylab.xlabel('Mole Fraction (ppm)',fontsize=18)
	pylab.ylabel('Pressure (bars)',fontsize=18)
	#pylab.title(title_text,fontsize=18)
	#fontstuff=FontProperties(size='smaller')
	
	pylab.legend(('H$_2$S','NH$_3$','H$_2$O','CH$_4$'),loc='lower left')
	
	#pylab.xticks(fontsize=18)
        #ax2=pylab.twiny()
	#pylab.gca().set_xscale('linear')
	
	#pylab.plot(tcme_2[:,1],tcme_2[:,0],'r:')
	
	
	
	#pylab.xlabel('Temperature ($^{\circ}$K)',fontsize=18)
	
	
	
      #change the yaxis to be log scale
	#pylab.gca().set_yscale('log');
	#pylab.ylim((0.1,200))
        #pylab.xlim((0,800))
	#pylab.yticks((0.1,1,10,100),('0.1','1','10','100'),fontsize=18)
	#pylab.ylabel('Pressure (bars)',fontsize=18)
	#pylab.xticks(fontsize=18)
	#pylab.gca().set_ylim(pylab.gca().get_ylim()[::-1])
	#pylab.legend(('T'),loc='upper right')

	pylab.savefig(plotName)

