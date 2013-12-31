        rosenkranz_abs=rosenkranz(F_model,asarray([h2o_density[i,0]]),asarray([[0.0]]),Experiment[i]['T'][0])
        deboer_abs=deboer_corrected(asarray([F_model]).T,asarray([Experiment[i]['T'][0]]),asarray([0.0]),asarray([0.0]),asarray([Experiment[i]['Pave'][0]]))
        
        figure(i)
        errorbar(Experiment[i]['nh3PeakFa'][:,0]/1e9,Experiment[i]['Absorption'][:,0,0],yerr=Experiment[i]['Absorption_2Sigma'][:,0],ecolor='k',fmt=None)
        
        plot(F_model,rosenkranz_abs.T,'k',label='Rosenkranz')
        plot(F_model.T,deboer_abs,'b',label='Deboer')
        xlabel(r'Frequency (GHz)')
        ylabel(r'Absorption $\left(\frac{\mathrm{dB}}{\mathrm{km}}\right)$')
        semilogy()
        legend()
        title('T= %6.3f'%(Experiment[i]['T'][0])+r'$^{\circ}$K        $\rho=$ %6.3f'%(h2o_density[i,0])+r'$\frac{\mathrm{g}}{\mathrm{m}^3}$         '+r'P$_{H_2O}$= %6.3f'%(Experiment[i]['Pave'][0])+' bars')
        xlim([1,7])
        
        if(i!=18):
                ylim([1e-3,1e1])
        else:
                ylim([1e-6,1e1])        
        savefig('Experiment_'+str(i+1)+'_pure.pdf')
