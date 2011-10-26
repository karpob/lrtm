def cloudPlot(tcme1,tcme2,title_text,plotName):
	import pylab
	# Just in case user decides to use LaTeX in title_text
	pylab.figure()	
 	pylab.rc('text',usetex=True)
 	ax=pylab.subplot(111)
      #set up the contour axes for easy adjustment
	pylab.plot(tcme1[:,11]*1e6,tcme1[:,0],'c',linewidth=1)
	pylab.plot(tcme1[:,13]*1e6,tcme1[:,0],'k',linewidth=1) 
	pylab.plot(tcme1[:,14]*1e6,tcme1[:,0],'g',linewidth=1)
	pylab.plot(tcme1[:,17]*1e6,tcme1[:,0],'m',linewidth=1)
	pylab.plot(tcme2[:,11]*1e6,tcme2[:,0],'c:',linewidth=1)
	pylab.plot(tcme2[:,13]*1e6,tcme2[:,0],'k:',linewidth=1) 
	pylab.plot(tcme2[:,14]*1e6,tcme2[:,0],'g:',linewidth=1)
	pylab.plot(tcme2[:,17]*1e6,tcme2[:,0],'m:',linewidth=1)

		
      # change the xaxis to be log scale
	pylab.gca().set_xscale('log');
	
	pylab.xlim((1e-4,2e2))
	pylab.ylim((0.1,10))		
	
	#Set the tick values
	# 10^(-1) 10^0 10^1 10^2 10^3 == Ugly
	# instead make the axis read 0.1, 1,10,100,1000
	pylab.yticks((0.1,1,5,10),('0.1','1','5','10'),fontsize=18)
	

	#reverse the order of the y axis (change so pressure goes up as you go down in y)
	

	#Label x and y axis
	pylab.xlabel('Cloud Density (g/m$^3$)',fontsize=17)
	pylab.ylabel('Pressure (bars)', fontsize=17)
	
	pylab.legend(('NH$_4$SH','NH$_3$','H$_2$O','NH$_{3}$-H$_2$O'),loc='lower left')
	pylab.xticks(fontsize=18)
	
	pylab.gca().set_ylim(pylab.gca().get_ylim()[::-1])
	pylab.savefig(plotName)
	
