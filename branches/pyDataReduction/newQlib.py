
def splineFun(freq,pp,peak):
        from scipy.interpolate import splev
        desired_point=peak
        
        calculated=splev(freq,pp)
        return calculated-desired_point
def processVSWR(S11,S21):
        import numpy
        S21=S21[S21==S21.max(axis=0)].mean()
        """
        idxs=[]
        idys=[]
        jj=0
        
        
        print maxV.shape
        for v in maxV:
                print jj
                idx=numpy.nonzero(v==S21[:,jj])
                idxs.append(idx)
                idys.append(jj)
                jj+=1
        """        
        o={} 
        
        RL=S11.min()
        VSWR=(1.+10**(RL/20))/(1-10**(RL/20))
        o['VSWR']=VSWR
        o['RL']=RL
        o['IL']=S21
        #o['VSWRe']=numpy.std(VSWR)
        return o
def lorentzFunction(p,f):
        
        
        gamma=p[0]
        fo=p[1]
        I=p[2]

        
        F=I*(gamma**2/((f-fo)**2+gamma**2))
        
        return F
def lorentzFunctionDf(p,f):
        gamma=p[0]
        fo=p[1]
        I=p[2]
        
        F=I*gamma**2*(-2*f+2*fo)/(gamma**2+(f-fo)**2)**2
        return F        
def lorentzFunctionPlusLinear(p,f):
        from scipy import polyval
        gamma=p[0]
        I=abs(p[1])
        fo=abs(p[0])
        
        
        
        #m=p[3]
        #b=p[4]
        
        F=I*(gamma**2/((f-fo)**2+gamma**2))+10**((polyval(p[1:len(p)],f)/10))

        return F



                     
def NAprocess(B,gamma,freq,Transmiss,r,experimentNumber,pressureLevel,oldB,oldF,oldQ,TypeProcess):
        from scipy.interpolate import splrep,splev,Rbf
        from numpy import mean,max,ones,asarray,zeros,linspace,nonzero,std,sqrt,interp,arange,log10,pi,diff,median,fliplr
        from scipy.optimize import fsolve,leastsq
        from scipy.interpolate import interp1d
        from scipy import polyfit,polyval
        import numpy as np
        import pylab
        
        """
        Algorithm for processing resonances taken on a network analyzer or other 
        device that does not have any sweep on scan nulls

        It compares the effect of the smoothing to the original data to make sure
        that the data isn't over smoothed, which would decrease the amplitude of
        the peak.  This algorithm tries highly smoothed values first and
        decreases by a factor of 10 each iteration until the absolute difference
        of the average of the smoothed versus unsmoothed peak values is less than
        0.01 dB.

        Inputs 
        B = Amplitude data in a matrix format with each sweep stored in its own column
        numswp = number of sweeps (typically 30)
        freq = 1601-point vector of frequency values corresponding to B
        Transmiss = 1601-point vector of transmissivity values used for
        deconvolution, if necessary

        Outputs
        BW3dBa BW3dBe = 3dB bandwidth mean and 2-sigma standard deviation, each
        with three values: the normal width, twice the higher half-width, and
        twice the lower half-width
        PeakFa PeakFe = peak frequency mean and 2-sigma standard deviation
        Qa Qe = Q mean and 2-sigma standard deviation, each with three values: 
        the normal Q, twice the higher half-width, and twice the lower half-width
        Levela Levele = amplitude mean and 2-sigma standard deviation
        Asymma Asymme = asymmetry mean and 2-sigma standard deviation
        """
        #print max(max(-20*log10(1+gamma.reshape([1601,1])))),min(min(-20*log10(1+gamma.reshape([1601,1]))))
        #tt=20.*log10(1.+10**(S11[0,:]/20.))
        print gamma.shape, Transmiss.shape
        alphaC1=fsolve(cableBalance,-1.*np.ones(Transmiss.shape[0]),args=(20.*np.log10(gamma.reshape([1601,1])[:,0]),Transmiss[:,0]))
        
        RL=20.*np.log10(gamma.reshape([1601,1])[:,0])-alphaC1
        
        gamma1=10.**(RL/20.)
        print 'gamma1',gamma1.shape,alphaC1.shape
        #print min(B),max(B)
        #B=10*log0(B)
        Boriginal=B  # maintains original data
        tt=10*log10(1.-gamma1*gamma1)
        #B=B-tt##+tt.reshape([1601,1])
        #for n in range(0,Boriginal.shape[1]):
        #        Boriginal[:,n]=Boriginal[:,n]        

        x=1.   #smoothing control parameter
        Pass=0
        mm,nn=B.shape
        peak=zeros(nn)
        threshold=zeros(nn)
        PeakF=zeros(nn)
        RawPeakF=zeros(nn)
        RawLevel=zeros(nn)
        fb=zeros(nn)
        fa=zeros(nn)
        ffa=zeros(nn)
        ffb=zeros(nn)
        PeakFF=zeros(nn)
        QQ=zeros(nn)
        QQh=zeros(nn)
        QQl=zeros(nn)
        
        bw3db=zeros(nn)
        bw3dbh=zeros(nn)
        bw3dbl=zeros(nn)
        Q=zeros(nn)
        Qh=zeros(nn)
        Ql=zeros(nn)
        asym=zeros(nn)
        locations=[]
        Ponumpoints=zeros([mm,nn])
        sweeps=nn
        N=160100
        #print freq.shape
        minFreq=min(freq[:,0])
        if(r==0):
                gamma0=12700.
        elif(r==1):
                gamma0=16300.
        elif(r==2):
                gamma0=39000.
        elif(r==3):
                gamma0=19000.
        elif(r==4):
                gamma0=24000.
        elif(r==5):
                gamma0=28800.
        elif(r==6):
                gamma0=22700.
        elif(r==7):
                gamma0=42000.
        elif(r==8):
                gamma0=33000.
        elif(r==9):
                gamma0=51000.
        elif(r==10):
                gamma0=44500.
        else:
                gamma0=181000.
        meanB=B.mean(axis=1)
        peakMeanB=max(meanB)
        
        idxP,=nonzero(meanB==peakMeanB)
        idxP=idxP[0]
        pF=freq[idxP,0]
        p=zeros([3])
        p[0]=gamma0#/(pi*max(B[:,s]))
        p[1]=pF
        p[2]=10.**(peakMeanB/10.)        
        #fit Center Lorentz distribution for resonance
        if (minFreq>4.3e9 and minFreq<5.7e9):idx=meanB>max(meanB)-5.
        else:idx=meanB>max(meanB)-3.0
        errFunctionLorentzOnly=lambda p,x,y:lorentzFunction(p,x)-y        
        pAll,stat=leastsq(errFunctionLorentzOnly,p,args=(freq[idx,0],10.**(meanB[idx]/10.)))
              
        #fit 2 nearby lorentz distributions holding center distribution constant       
        errFunctionSides=lambda p,x,y,pp:lorentzFunction(pp[0:3],x)+lorentzFunction(p[0:3],x)+lorentzFunction(p[3:6],x)-y
        actualFunctionSides=lambda p,x,pp:lorentzFunction(pp[0:3],x)+lorentzFunction(p[0:3],x)+lorentzFunction(p[3:6],x)
                
        p2=zeros([6])
        #set initial conditions to the same as center distribution        
        p2[0:3]=pAll
        p2[3:6]=pAll
                
        #set intial conditions to 2 times peak lower "resonance" and 0.5 peak for upper "resonance"        
        p2[0]=2*p2[0]
        p2[3]=0.5*p2[0]
        #put the center of the distributions -2 center frequency +2 center frequency
        p2[1]=p2[1]-2.*p2[0]
        p2[4]=p2[4]+2.*p2[3]
                
        # fit the error function        
        pFinal,stat=leastsq(errFunctionSides,p2,args=(freq[:,0],10.**(meanB/10.),pAll))
        
        #pFirstGuess=pAll
        #pSecondGuess=pFinal
        
        #fit side and center frequency distributions together using previous fitting as initial conditions
        ppAll=zeros([9])
        ppAll[0:3]=pAll
        ppAll[3:6]=pFinal[0:3]
        ppAll[6:9]=pFinal[3:6]
                
        errFunctionAll=lambda p,x,y:lorentzFunction(p[0:3],x)+lorentzFunction(p[3:6],x)+lorentzFunction(p[6:9],x)-y
        actualFunctionAll=lambda p,x:lorentzFunction(p[0:3],x)+lorentzFunction(p[3:6],x)+lorentzFunction(p[6:9],x)
                
        #fit using entire spectra        
        ppFinal,stat=leastsq(errFunctionAll,ppAll,args=(freq[:,0],10.**(meanB/10.)))
        
        idx=meanB>max(meanB)-5.0
        FreqFineGrid=linspace(min(freq[idx,0]),max(freq[idx,0]),N)
        BfineGrid=10.*log10(actualFunctionAll(ppFinal,FreqFineGrid))
        try:PeakFall=FreqFineGrid[max(BfineGrid)==BfineGrid]#o[0]#
        except:
                idx,=nonzero(max(BfineGrid)==BfineGrid)
                PeakFall=FreqFineGrid[median(idx)]
                
        solveForFreq=lambda f,p,y:10.*log10(actualFunctionAll(p,f))-y
        o=fsolve(solveForFreq,FreqFineGrid[N/4],args=(ppFinal,max(BfineGrid)-3.0))
        ffa=o[0]
        o=fsolve(solveForFreq,FreqFineGrid[N-N/4],args=(ppFinal,max(BfineGrid)-3.0))
        ffb=o[0]
        Qall=PeakFall/(ffb-ffa)
        
        BBfinalZ=np.zeros([1601,30])
        fiG=pylab.figure()                                                                                                
        for s in range(0,sweeps):
                
                
                B[:,s]=B[:,s]#-tt
                        
                peak[s]=max(B[:,s])
                RawLevela=max(B[:,s])
                RawPeakF[s]=freq[B[:,s]==max(B[:,s])]
                RawLevel[s]=max(B[:,s])
                location2=nonzero(B[:,s]==peak[s]) 
                location=location2
                locations.append(location)
                # pins the amplitude of the transmissivity
                # deconvolution to zero at the peak of the resonance       
                ZeroVal=Transmiss[location]*ones([len(freq),1])
                ZeroVal2=tt[location]*ones([len(freq),1]) 
                #ZeroValMismatch=10.*np.log10(1.-gamma.reshape([1601,1])[location]**2)
                #Mismatch=10.*np.log10(1.-gamma.reshape([1601,1])[location]**2)        
                # deconvolves the effect of transmissivity of the cables
                print tt.shape,ZeroVal2.shape
                B[:,s]=B[:,s]-Transmiss[:,0]-tt+ZeroVal2[:,0]+ZeroVal[:,0]
                
                #start with reasonable smoothing factor, decrease if results in poor fit to peak 
                     
                #p[3]=0.0
                       
                        
                        
                peakIdx=nonzero(B[:,s]==max(B[:,s]))
                LH=10.**(B[0:len(B[:,s])/2,s]/10)
                FrequencyLH=freq[0:len(B[:,s])/2,0]
                FrequencyRH=freq[len(B[:,s])/2:len(freq),0]
                RH=10.**(B[len(B[:,s])/2:len(freq),s]/10.)
                
                
                #fit using entire spectra        
                ppFinal,stat=leastsq(errFunctionAll,ppFinal,args=(freq[:,0],10.**(B[:,s]/10.)))
                
                #fit using -3dB down spectra for good measure
                #if(minFreq<6.0e9):
                #idx=idx=idx=B[:,s]>max(B[:,s])-3.1
                #ppFinal,stat=leastsq(errFunctionAll,ppFinal,args=(freq[idx,0],10.**(B[idx,s]/10.)))
                
                #fill a vector of the fitted spectra
                BBfinal=10.*log10(actualFunctionAll(ppFinal,freq[:,0]))                      
               
                #uncomment block below for plotting/debugging
                
                idxX=B[:,s]>max(B[:,s])-3.2
                pylab.subplot(311)
                pylab.plot(freq[:,0]/1e9,B[:,s]-BBfinal,'kx')
                pylab.plot(freq[:,0]/1e9,B[:,s]-oldB[:,s],'bx')
                pylab.subplot(312)
                pylab.plot(freq[idxX,0]/1e9,B[idxX,s]-BBfinal[idxX],'kx')
                pylab.plot(freq[idxX,0]/1e9,B[idxX,s]-oldB[idxX,s],'bx')
                
                
                """
                pylab.figure(1)
                
                #pylab.plot(freq[idx,0],B[idx,s],'b')
                #pylab.plot(freq[idx,0],Boriginal[idx,s],'r')
                #pylab.plot(freq[idx,0],allZ[idx],'r')
                #pylab.plot(freq[idx,0],Bfinal[idx],'g')
                #pylab.plot(freq[idx,0],BBfinal[idx],'k')
                pylab.plot(freq,tt)
                pylab.figure(2)
                
                #pylab.plot(freq[idx,0],B[idx,s]-allZ[idx],'r')
                #pylab.plot(freq[idx,0],B[idx,s]-Bfinal[idx],'b')
                pylab.plot(freq[idx,0],B[idx,s]-BBfinal[idx],'k')
                """
                
                
                
                
                          
                FreqFineGrid=linspace(min(freq[idx,0]),max(freq[idx,0]),N)
                
                
                BfineGrid=10.*log10(actualFunctionAll(ppFinal,FreqFineGrid))
                
                #errPeakFind=lambda f,p:lorentzFunctionDf(p[0:3],f)+lorentzFunctionDf(p[3:6],f)+lorentzFunctionDf(p[6:9],f)
                #o=fsolve(errPeakFind,FreqFineGrid[max(BfineGrid)==BfineGrid],ppFinal)
                
                try:PeakF[s]=FreqFineGrid[max(BfineGrid)==BfineGrid]#o[0]#
                except:
                        idx,=nonzero(max(BfineGrid)==BfineGrid)
                        PeakF[s]=FreqFineGrid[median(idx)]
                #solve for 3dB points
                #idxPeak,=nonzero(BfineGrid==max(BfineGrid))
                #LH=Rbf(BfineGrid[0:idxPeak+1],FreqFineGrid[0:idxPeak+1],)
                #RH=Rbf(fliplr(BfineGrid[idxPeak:len(BfineGrid)]),fliplr(FreqFineGrid[idxPeak:len(BfineGrid)]))
                #fa[s]=LH(BfineGrid[idxPeak]-3.0)
                #fb[s]=RH(BfineGrid[idxPeak]-3.0)
                
                solveForFreq=lambda f,p,y:10.*log10(actualFunctionAll(p,f))-y
                o=fsolve(solveForFreq,FreqFineGrid[N/4],args=(ppFinal,max(BfineGrid)-3.0))
                fa[s]=o[0]
                o=fsolve(solveForFreq,FreqFineGrid[N-N/4],args=(ppFinal,max(BfineGrid)-3.0))
                fb[s]=o[0]
                
                peak[s]=max(BfineGrid)
                        
                bw3db[s]=fb[s]-fa[s] # calculate 3dB bandwidth
                
                        # don't need to include effect of RBW here, since that effect is not
                        # present with network analyzer measurements

                        # used in calculating asymmetry uncertainty
                bw3dbh[s]=2.*(fb[s]-PeakF[s])
                
                        # bandwidth if doubling higher frequency half BW

                bw3dbl[s]=2.*(PeakF[s]-fa[s])
                
                        # bandwidth if doubling lower frequency half BW

                Q[s]=PeakF[s]/bw3db[s]
                Qh[s]=PeakF[s]/bw3dbh[s]
                Ql[s]=PeakF[s]/bw3dbl[s]
                        
                asym[s]=100.*((fb[s]-PeakF[s])-(PeakF[s]-fa[s]))/(fb[s]-fa[s])
                        # from DeBoer and Steffes, 1996
                BBfinalZ[:,s]=BBfinal

                

        pylab.subplot(313)
        pylab.plot(freq[idxX,0]/1e9,B[idxX,:],'rx')
        pylab.plot(freq[idxX,0]/1e9,BBfinalZ[idxX,:].mean(axis=1),'k')
        pylab.plot(freq[idxX,0]/1e9,oldB[idxX,:].mean(axis=1),'b')         
        #pylab.show()
        # calculates the mean and sample standard deviation of the values to return
        o={} #create output dict
        o['Qall']=Qall
        o['PeakFall']=PeakFall
        o['BW3dBa']=asarray([mean(bw3db), mean(bw3dbh), mean(bw3dbl)])
        o['BW3dBe']=asarray([std(bw3db), std(bw3dbh), std(bw3dbl)])
        o['PeakFa']=mean(PeakF)
        o['PeakFe']=std(PeakF)
        o['RawPeakFa']=mean(RawPeakF)
        o['RawPeakFe']=std(RawPeakF)
        o['RawLevela']=mean(RawLevel)
        o['RawLevele']=std(RawLevel)
        o['Qa']=asarray([mean(Q), mean(Qh), mean(Ql)])
        o['Qe']=asarray([std(Q), std(Qh), std(Ql)])
        
        o['Levela']=mean(peak)
        o['Levele']=std(peak)
        o['Asymma']=mean(asym)
        o['Asymme']=std(asym)
        print gamma.reshape([1601,1])[idxP],idxP
        o['cableTrans']=Transmiss[idxP]
        o['reflectionGamma']=gamma.reshape([1601,1])[idxP]
        #print o['PeakFa']
        #print o['PeakFe']
        #print o['Qa']
        #print o['Qe']
        #print o['PeakFe']
        pylab.xlabel('Frequency (GHz)')
        pylab.ylabel('Residual (dB)')
        #print oldQ,o['Qa'][0],oldF,o['PeakFa']
        fiG.suptitle('Old-New deltaQ,deltaF %f,%f'%(float(oldQ)-float(o['Qa'][0]),float(oldF)-float(o['PeakFa'])))
        pylab.savefig(TypeProcess+'Residuals_experiment_%02d_pressure_%02d_r_%02d.png'%(experimentNumber,pressureLevel,r))
        pylab.close(fiG)
             
        return o

def unloadedQ(insertionLossFo,loadedQ,alpha,returnLossFo):
        """
        Get the unloaded Q  from Cavity resonance.(circuits not dogmatic PAL notation)
        insertionLossFo--> measured S21 at resonance. (dB)
        loadedQ --> measured fo/BW 
        """          
                      
        insertionLossFo=insertionLossFo-alpha
        import numpy as np
        logS21return=20.*np.log10(1.-10**(returnLossFo/20.))
        
        linS21=np.power(10.,(insertionLossFo+logS21return)/20.) #note that 20 is correct need voltage for Q and coupling stuff
        
        unloadedQ=loadedQ/(1.-linS21)
        return unloadedQ
def calculateAbsorption(Freq,insertionLossFo,loadedQgas,Qmatch):
        """
        get Absorption coefficient from unloaded measured Q and dielectrically matched Q.
        """
        import numpy as np     
        dBconvert=2.*10.*np.log10(np.exp(1.0))
        totalAbs=0.0
        lamb=299792.458/Freq
        #km/s
        ThresholdFlag=True
        
        while(ThresholdFlag):
                totalAbsPrev=totalAbs
                Qgas=unloadedQ(insertionLossFo,loadedQgas,alpha=-1.0*totalAbsPrev)  
                alpha=dBconvert*(np.pi/lamb)*(1./Qgas-1./Qmatch)
                EPL=(lamb*Qgas)/(2.*np.pi)
                
                if(alpha>0.0):
                        totalAbs=EPL*alpha
                else:
                        totalAbs=0.        
                if(abs((totalAbs-totalAbsPrev)/totalAbs)*100<0.1):
                        ThresholdFlag=False
        return alpha
def cableBalance(alpha,RL,Tran):
        import numpy as np
        return Tran-alpha-10.*np.log10(1.-np.power((np.power(10,((RL-alpha)/20.))),2))        
def calculateAbsorptionAndQgas(Freq,loadedQgas,loadedQmatch,insertionLossFoGas,insertionLossFoDiel,cableLossFoGas,cableLossFoDiel,gammaGas,gammaDiel):
        #(Freq,insertionLossFo,loadedQgas,Qmatch,returnLossFo):
        from scipy.optimize import fsolve
        """
        get Absorption coefficient from unloaded measured Q and dielectrically matched Q.
        """
        import numpy as np
        extra=0.5    
        dBconvert=2.*10.*np.log10(np.exp(1.0))
        totalAbs=1e-10
        lamb=299792.458/Freq
        alphaC1=fsolve(cableBalance,-1.,args=(20.*np.log10(gammaGas),cableLossFoGas))
        
        alphaC2=fsolve(cableBalance,-1.,args=(20.*np.log10(gammaDiel),cableLossFoDiel)) 
        RLgas=20.*np.log10(gammaGas)-alphaC1
        RLdiel=20.*np.log10(gammaDiel)-alphaC2
        gammaGas=10.**(RLgas/20.)
        gammaDiel=10.**(RLdiel/20.)
        
        mismatchLossGas=10.*np.log10(1.-gammaGas**2)
        mismatchLossDiel=10.*np.log10(1.-gammaDiel**2)
        print 'ac1',alphaC1,cableLossFoGas,gammaGas,gammaDiel,mismatchLossGas,mismatchLossDiel
        if (mismatchLossGas<-1.0): print mismatchLossGas,mismatchLossDiel,gammaGas,gammaDiel
        #linS21Gas=np.sqrt(np.power(10.,(insertionLossFoGas-alphaC1-mismatchLossGas)/10.)) #note that 20 is correct need voltage for Q and coupling stuff
        linS21Gas=np.sqrt(np.power(10.,(insertionLossFoGas-cableLossFoGas)/10.)) #note that 20 is correct need voltage for Q and coupling stuff
        linS21Diel=np.sqrt(np.power(10.,(insertionLossFoDiel-cableLossFoDiel)/10.)) #note that 20 is correct need voltage for Q and coupling stuff
        #linS21Diel=np.sqrt(np.power(10.,(insertionLossFoDiel-alphaC2-mismatchLossDiel)/10.)) #note that 20 is correct need voltage for Q and coupling stuff
        print 'mm',linS21Gas,linS21Diel,gammaDiel,gammaGas,mismatchLossGas,mismatchLossDiel
        alpha=dBconvert*(np.pi/lamb)*((1.-linS21Gas)/loadedQgas-(1.-linS21Diel)/loadedQmatch)
        Qgas=0.0#1./(1./loadedQgas-1./loadedQmatch)*(1.-linS21Gas)
        """
        ThresholdFlag=True
        Qgas=unloadedQ(insertionLossFo,loadedQgas,0.0,returnLossFo)  
        alpha=dBconvert*(np.pi/lamb)*(1./Qgas-1./Qmatch)
        
        while(ThresholdFlag):
                totalAbsPrev=totalAbs
                
                Qgas=unloadedQ(insertionLossFo,loadedQgas,-1.0*totalAbsPrev)  
                alpha=dBconvert*(np.pi/lamb)*(1./Qgas-1./Qmatch)
                EPL=(lamb*Qgas)/(2.*np.pi)
                
                if(alpha>0.0):
                        totalAbs=EPL*alpha
               
                if(abs((totalAbs-totalAbsPrev)/totalAbs)*100<0.1):
                        ThresholdFlag=False
        """                
        return alpha,Qgas
        
def caculateErrInst(Ql,Qu,sl,su,fol,fou,BW,numswp,timez):
        """
        function sigma_n=insterror(Ql,Qu,sl,su,Span,RBW,fol,fou,t_loaded,t_unloaded,BW,numswp,time,device)
        Error due to electrical noise and instrumentation uncertainity
        HP 8564E calibrated 12/5/05

        Courtesy of Priscilla Mohammed, modified for network analyzer by Tom Hanley

        Ql = Q loaded
        Qu = Q matched
        sl = Sample standard deviation of BW loaded readings in Hz
        su = Sample standard deviation of BW matched readings in Hz
        fol = Freq of loaded resonance Hz
        fou = Freq of matched resonance Hz
        
        BW = BW of loaded resonance Hz
        numswp = number of sweeps (for finding confidence coefficient)
        time = seconds from epoch of measurement
       
        critical values from t-table for 95% confidence
        """ 
        from calendar import timegm
        from numpy import sqrt,log10,exp,pi
        if numswp==5:confc=2.776
        elif numswp==3:confc=4.303
        elif numswp==7:confc=2.447    
        elif numswp==10:confc=2.262
        elif numswp==30:confc=2.045 

    
        sigma_nl_sqrd = (confc*(sl))**2/numswp #Sigma loaded BW weighted by 95% confidence coefficient
        sigma_nu_sqrd = (confc*(su))**2/numswp #Sigma matched BW weighted by 95% confidence coefficient

        lamb = 299792.458/fol;   #Wavelength of loaded resonance in km


        # Use 95% confidence 

        
        yrs = (timez-timegm((2007, 1, 7, 0, 0, 0)))/(60*60*24*365.24) # years since E5071C network analyzer calibrated (1/7/07)
        sigma_delta = (BW*sqrt(2)*(5e-8+5e-7*yrs))*2/3 #Hz
        sigma_o = fol*(5e-8+5e-7*yrs)*2/3


        gamma_u = 1.# - sqrt(t_unloaded)
        gamma_l = 1.# - sqrt(t_loaded)

        A = (gamma_u**2)/(fou**2)
        B = (sigma_o**2)/(Qu**2)
        C = (2*sigma_o*sigma_delta)/Qu
        Gamma_u_sqrd = A*(B + sigma_delta**2 + sigma_nu_sqrd + C)

        a = (gamma_l**2)/(fol**2)
        b = (sigma_o**2)/(Ql**2)
        c = (2*sigma_o*sigma_delta)/Ql
        Gamma_l_sqrd = a*(b + sigma_delta**2 + sigma_nl_sqrd + c)

        D = -(gamma_u*gamma_l)/(fol*fou);
        E = (sigma_o**2)/(Ql*Qu);
        F = (sigma_o*sigma_delta)/Ql;
        G = (sigma_o*sigma_delta)/Qu;
        Gamma_u_Gamma_l = D*(E + sigma_delta**2 + F + G)

        sigma_psi_sqrd = Gamma_u_sqrd + Gamma_l_sqrd - 2*(Gamma_u_Gamma_l);
        sigma_n = (20*log10(exp(1))*pi/lamb*(sqrt(sigma_psi_sqrd))); #dB/km for 2-sigma uncertainty.
        return sigma_n
def calculateErrDiel(Qvac,Qdiel,peakFvac,peakFdiel,peakFmeas,insertionLossFo,loadedQgas):
        import numpy
        dQdf=numpy.abs((Qvac-Qdiel)/(peakFvac-peakFdiel))
        dQ=numpy.max(dQdf)*(peakFmeas-peakFdiel)
        Err=calculateAbsorption(peakFmeas,insertionLossFo,loadedQgas,Qdiel+dQ)-calculateAbsorption(peakFmeas,insertionLossFo,loadedQgas,Qdiel-dQ)
        return Err
def calculateErrTrans(loadedQdiel,loadedQgas,insertionLossFoGas,insertionLossFoDiel,sigmaGas,sigmaDiel,peakFmeas):
                                
        QdielPlus=unloadedQ(insertionLossFoDiel+sigmaDiel,loadedQdiel)
        QdielMinus=unloadedQ(insertionLossFoDiel-sigmaDiel,loadedQdiel)
        AlphaPlus=calculateAbsorption(peakFmeas,insertionLossFoGas+sigmaGas,loadedQgas,QdielPlus)
        AlphaMinus=calculateAbsorption(peakFmeas,insertionLossFoGas-sigmaGas,loadedQgas,QdielMinus)
        Err=abs(AlphaPlus-AlphaMinus)
        return Err
def calculateErrAsym(QgasH,QgasL,QdielL,QdielH,insertionLossFoDiel,insertionLossFoGas,peakFmeas):
        QdielPlus=unloadedQ(insertionLossFoDiel,QdielH)
        QdielMinus=unloadedQ(insertionLossFoDiel,QdielL)
        
        AlphaPlus=calculateAbsorption(peakFmeas,insertionLossFoGas,QgasH,QdielPlus)
        AlphaMinus=calculateAbsorption(peakFmeas,insertionLossFoGas,QgasL,QdielMinus)
        Err=abs(AlphaPlus-AlphaMinus)
        return Err        
def calculateAbsorptionAndErrs(oGas,oDiel,oVac,timez,vswrGas,vswrDiel):
        from calendar import timegm
        import numpy
        Qvacs=[]
        Fvacs=[]
        Lvacs=[]
        #for item in oVac:
        #        Qvacs.append(item['QQa'][0])
        #        Fvacs.append(item['RawPeakFa'])
        #        Lvacs.append(item['RawLevela'])
        #print 'err?',Fvacs[0],oDiel['RawPeakFa']        
        #Qmatch=unloadedQ(oDiel['Levela'],oDiel['Qa'][0],0.0,vswr['RL'])
        #Qvac=unloadedQ(numpy.asarray(Lvacs),numpy.asarray(Qvacs),0.0)
        #print oGas['Levela']
        #alpha,Qgas=calculateAbsorptionAndQgas(oGas['PeakFa'],oGas['Levela'],oGas['Qa'][0],Qmatch,oDiel['Levela'])
        #alpha,Qgas=calculateAbsorptionAndQgas(oGas['PeakFa'],oGas['Qa'][0],oDiel['Qa'][0],oDiel['Levela'],vswr['RL'])
        
        alpha,Qgas=calculateAbsorptionAndQgas(oGas['PeakFa'],oGas['Qa'][0],oDiel['Qa'][0],oGas['Levela'],oDiel['Levela'],oGas['cableTrans'],oDiel['cableTrans'],oGas['reflectionGamma'],oDiel['reflectionGamma'])
        beta1=Qgas/oGas['Qa'][0]
        beta2=beta1
        #beta2=Qmatch/oDiel['Qa'][0]
        #ErrInst=caculateErrInst(Qgas,Qmatch,oGas['BW3dBe'][0],oDiel['BW3dBe'][0],oGas['PeakFa'],oDiel['PeakFa'],oGas['BW3dBa'][0],30.,timez)
        #ErrDiel=calculateErrDiel(Qvac,Qmatch,numpy.asarray(Fvacs),oDiel['PeakFa'],oGas['PeakFa'],oGas['Levela'],oGas['Qa'][0])
        #ErrTrans=calculateErrTrans(oDiel['Qa'][0],oGas['Qa'][0],oGas['Levela'],oDiel['Levela'],oGas['Levele'],oDiel['Levele'],oGas['PeakFa'])
        #ErrAsym=calculateErrAsym(oGas['Qa'][1],oGas['Qa'][2],oDiel['Qa'][1],oDiel['Qa'][2],oDiel['Levela'],oGas['Levela'],oGas['PeakFa'])
        
        o={}
        o['AbsorptionK']=alpha
        #o['Qgas']=Qgas
        #o['beta1']=beta1
        #o['beta2']=beta2
        #o['ErrInstK']=ErrInst
        #o['ErrDielK']=ErrDiel
        #o['ErrTransK']=ErrTrans
        #o['ErrAsymK']=ErrAsym
        #o['Absorption_2sigmaK']=numpy.sqrt(ErrInst**2+ErrDiel**2+ErrTrans**2+ErrAsym**2)
        return o                
