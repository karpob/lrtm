def main():
        
        
        dirList=['1_375_0.33_h2o_20bar_max',
                 '2_375_0.33_h2o_86bar_max',
                 '3_375_0.33_h2o_6bars',
                 '4_375_0.33_h2o_13.5bars',
                 '5_375_0.33_h2o_he_13.5bars_hygrometer',
                 '6_375_0.33_h2o_h2_only',
                 '7_450_1_2_bar_h2o_h2he_premix',
                 '8_450_1bar_h2o_he_13bars',
                 '9_525_3bar_h2o_h2he_premix',
                 '10_525_1_8_bar_h2o_he_13bars',
                 '11_500_2bar_h2o_pure_h2',
                 '12_500_1bar_13.5_he',
                 '13_500_2bar_6bar_he',
                 '14_450_1bar_h2o_h2_only',
                 '15_450_2bar_6bars_he',
                 '16_450_1bar_13bars_he',
                 '17_333_0.1bar_13bar_he',
                 '18_375_0.5bar_13bar_he',
                 '19_333_0.05bar_13bar_he',
                 '20_500K_2bH2O_14bHe',
		 '21_500K_1.5bH2O_6bHe'
                 ]

        GasFiles,GasDir,ArgonFiles,ArgonDir,TransFiles,TransDir,VacFiles,VacDir=getFileLists(dirList)
        
        VacDir,VacFiles=sortVacFiles(VacDir,VacFiles)
        makePickle(VacDir,VacFiles,'allVacSpectra.pkl','allVacAxis.pkl','allVacTimes.pkl','allVacS11.pkl')
        
                           
        ArgonDir,ArgonFiles,ArgonPressureCounts=sortArgonFiles(ArgonDir,ArgonFiles)                                      
        GasDir,GasFiles,PressureCounts=sortGasFiles(GasDir,GasFiles)
        TransDir,TransFiles,PressureCounts=sortTransFiles(TransDir,TransFiles)

        

        makePickle(GasDir,GasFiles,'allGasSpectra.pkl','allGasSpectraAxis.pkl','allGasTimes.pkl','allGasS11.pkl')
        makePickle(ArgonDir,ArgonFiles,'allDielSpectra.pkl','allDielSpectraAxis.pkl','allDielTimes.pkl','allDielS11.pkl')
        makePickleTrans(TransDir,TransFiles,'allTransSpectra.pkl')

def getFileLists(dirList):
        import os
        ArgonFiles=[]
        ArgonDir=[]
        GasFiles=[]
        GasDir=[]
        VacFiles=[]
        VacDir=[]
        TransFiles=[]
        TransDir=[]
        namez=[]
        dirz=[]
        for dirs in dirList:
                matlabFilez=os.listdir(dirs+'/'+'raw/')
                for filez in matlabFilez:
                        dirz.append(dirs)
                namez.extend(matlabFilez)
        

        ii=0
        for name in namez:
                a=name.split('_')
                dirVal=dirz[ii]
                #useful files start with resonance number
                try:
                        if(int(a[0])<13):
                                if('Ar' in a):
                                        ArgonFiles.append(name)
                                        ArgonDir.append(dirVal)
                                elif('trans' in name):
                                        TransFiles.append(name)
                                        TransDir.append(dirVal)
                                elif('vac' in name):
                                        VacFiles.append(name)
                                        VacDir.append(dirVal)
                                else:
                                        GasFiles.append(name)
                                        GasDir.append(dirVal)                                        
                except:pass
                ii+=1
        return GasFiles,GasDir,ArgonFiles,ArgonDir,TransFiles,TransDir,VacFiles,VacDir                
def sortArgonFiles(ArgonDir,ArgonFiles):
        import numpy
        ii=0    
        #open all measured Gas files, and get time stamp
        ResNumbers=[]

        crapolas=[]
        for filez in ArgonFiles:
                crapola=filez.split('_')
                ResNumbers.append(crapola[0])
                crapolas.append(crapola)
        usefulItems=[]        
        for crapola in crapolas:
                useful=[]
                for item in crapola[1:len(crapola)]:
			
                        if(item=='' or item=='Ar' or item=='match' or item=='h2o'):pass
                        else:
                                if(item[len(item)-4:len(item)]=='.mat'):
                                        item=item[0:len(item)-4]
                                if(item=='h2o' or item==''):pass        
                                else:
                                        if(item[len(item)-4:len(item)]=='bars'):
                                                item=item[0:len(item)-4]
                                        if(item=='pure'):
                                                item=1
                                        if(item=='backside' or item=='backwards'):
                                                item=-1

                                        if(item=='13bars2'):
                                                item=13
                                        if(item=='pure3'):
                                                item=1                
                                        if(item=='exactly'):
                                                item=5
                                        if(item==''):pass
                                        else:                                
                                                useful.append(item)
                usefulItems.append(useful)
        ii=0
        item=[]
        for filez in ArgonFiles:
                if(len(usefulItems[ii])>1):
                        item=usefulItems[ii]
                        if int(item[0])<0:
                                item=[100+100-int(item[1])]
                        elif int(item[1])<0:
                                item=[100+100-int(item[0])]
                        elif item[0]==5 or item[1]==5:
                                item=[int(item[0])+int(item[1])]
                                
                        usefulItems[ii]=item
                                        
        #print ArgonDir[ii],filez,usefulItems[ii]
                ii+=1
        ii=0
        pressureOrder=usefulItems
        keys=[]        
        for filez in ArgonFiles:
                ExperimentNumber=int(ArgonDir[ii].split('_')[0])
                #print ExperimentNumber,pressureOrder[ii],ResNumbers[ii]
                key=ExperimentNumber*1000000+1000*int(pressureOrder[ii][0])+int(ResNumbers[ii])
                keys.append(key)
                ii+=1
        idx=numpy.argsort(numpy.asarray(keys))
        ArgonFiles=numpy.asarray(ArgonFiles)[idx]
        ArgonDir=numpy.asarray(ArgonDir)[idx]
        keys=numpy.asarray(keys)[idx]
        ResNumbers=numpy.asarray(ResNumbers)[idx]
        PressureCounts=[]
        ii=0
        ll=[]
        for ii in range(0,len(ArgonDir)-1):
                if(ArgonDir[ii]==ArgonDir[ii+1]):
                        ll.append(pressureOrder[ii][0])
                                
                else:
                        PressureCounts.append((len(ll)+1)/12)
                        ll=[]
                ii+=1
        return ArgonDir,ArgonFiles,PressureCounts

        
def sortGasFiles(GasDir,GasFiles):
        import numpy
        ii=0    
        #open all measured Gas files, and get time stamp
        ResNumbers=[]

        crapolas=[]
        for filez in GasFiles:
                G=GasDir[ii]
                crapola=filez.split('_')
                
                #ResNumbers.append(crapola[0])
                
                crapolas.append(crapola)
        usefulItems=[]        
        for crapola in crapolas:
                useful=[]
                
                for item in crapola[1:len(crapola)]:
			
                        if(item[len(item)-4:len(item)]=='.mat'):
                                item=item[0:len(item)-4]
                        
                        if(item[len(item)-4:len(item)]=='bars'):
                                item=item[0:len(item)-4]
                        if(item=='pure'):
                                item=1
                        if(item=='backside' or item=='backwards'):
                                item=-1
				
                        if(item=='13bars2'):
                                item=13
                        if(item=='pure3'):
                                item=1                
                        if(item=='exactly'):
                                item=5
                        if(not item=='h2o'):
                                if(not item==''):
                                        useful.append(item)
                                        
                usefulItems.append(useful)
                ResNumbers.append(crapola[0])
                
        for iii in usefulItems:
                
                if(len(usefulItems[ii])>1):
                        item=usefulItems[ii]
                        #print item
                        if int(item[0])<0:
                                item=[100+100-int(item[1])]
                        elif int(item[1])<0:
                                item=[100+100-int(item[0])]
                        elif item[0]==5 or item[1]==5:
                                item=[int(item[0])+int(item[1])]
                                
                        usefulItems[ii]=item
                                        
                ii+=1
        ii=0
        pressureOrder=usefulItems
        keys=[]        
        for filez in GasFiles:
                ExperimentNumber=int(GasDir[ii].split('_')[0])
                key=ExperimentNumber*1000000+1000*int(pressureOrder[ii][0])+int(ResNumbers[ii])
                #print filez,key,pressureOrder[ii][0],pressureOrder[ii]
		print filez, ExperimentNumber, pressureOrder[ii][0],ResNumbers[ii]
                keys.append(key)
                ii+=1
               
        idx=numpy.argsort(numpy.asarray(keys))
        GasFiles=numpy.asarray(GasFiles)[idx]
        GasDir=numpy.asarray(GasDir)[idx]
        keys=numpy.asarray(keys)[idx]
        ResNumbers=numpy.asarray(ResNumbers)[idx]
        PressureCounts=[]
        ii=0
        ll=[]
        for ii in range(0,len(GasDir)-1):
                if(GasDir[ii]==GasDir[ii+1]):
                        ll.append(pressureOrder[ii][0])
                                
                else:
                        PressureCounts.append((len(ll)+1)/12)
                        ll=[]
                ii+=1    
        return GasDir,GasFiles,PressureCounts
        
def sortTransFiles(TransDir,TransFiles):
        import numpy
        ii=0    
        #open all measured Trans files, and get time stamp
        ResNumbers=[]

        crapolas=[]
        for filez in TransFiles:
                G=TransDir[ii]
                crapola=filez.split('_')
                
                ResNumbers.append(crapola[0])
                
                crapolas.append(crapola)
        usefulItems=[]        
        for crapola in crapolas:
                useful=[]
                for item in crapola[1:len(crapola)]:
                        if(item[len(item)-4:len(item)]=='.mat'):
                                item=item[0:len(item)-4]
                        
                        if(item[len(item)-4:len(item)]=='bars'):
                                item=item[0:len(item)-4]
                                
                        if(item=='pure'):
                                item=1
                        if(item=='backside' or item=='backwards'):
                                item=-1
                        if(item=='13bars2'):
                                item=13
                        if(item=='pure3'):
                                item=1                
                        if(item=='exactly'):
                                item=5
                        if(not item=='h2o'):
                                if(not item==''):
                                        useful.append(item)
                                        
                usefulItems.append(useful)
                
                
        for filez in TransFiles:
                
                if(len(usefulItems[ii])>1):
                        item=usefulItems[ii]
                        #print filez,item
                        if int(item[1])<0:
                                itemy=100+100-int(item[2])
                                
                        elif len(item)>2:
                                
                                if(int(item[2])<0):
                                        itemy=100+100-int(item[1])
                                if item[1]==5 or item[2]==5:
                                        itemy=int(item[1])+int(item[2])
                        
                        else:
                             itemy=int(item[1])
                        itemx=int(item[0][5])                
                        usefulItems[ii]=[itemx,itemy]
                                        
                ii+=1
        ii=0
        pressureOrder=usefulItems
        keys=[]        
        for filez in TransFiles:
                ExperimentNumber=int(TransDir[ii].split('_')[0])
                
                key=ExperimentNumber*1000000+1000*int(pressureOrder[ii][1])+int(ResNumbers[ii])*10+int(pressureOrder[ii][0])
                
                keys.append(key)
                ii+=1
               
        idx=numpy.argsort(numpy.asarray(keys))
        TransFiles=numpy.asarray(TransFiles)[idx]
        TransDir=numpy.asarray(TransDir)[idx]
        keys=numpy.asarray(keys)[idx]
        ResNumbers=numpy.asarray(ResNumbers)[idx]
        ii=0
        ll=[]
        TransCounts=[]
        for ii in range(0,len(TransDir)-1):
                if(TransDir[ii]==TransDir[ii+1]):
                        ll.append(pressureOrder[ii][0])
                                
                else:
                        TransCounts.append(max(ll))
                        ll=[]
                ii+=1
                                       
        return TransDir,TransFiles,TransCounts
def sortVacFiles(VacDir,VacFiles):
        import numpy
        ii=0
        keys=[]
        for vac in VacFiles:
                
                crapola=vac.split('_')
                crapola2=VacDir[ii].split('_')
                Experiment=crapola2[0]
                Res=crapola[0]
               
                vacNum=crapola[1][3]
                if(vacNum=='l'):vacNum=1
                keys.append(int(Experiment)*10000+int(vacNum)*100+int(Res))
                ii+=1
        idx=numpy.argsort(numpy.asarray(keys))
        VacFiles=numpy.asarray(VacFiles)[idx]
        VacDir=numpy.asarray(VacDir)[idx]
        return VacDir,VacFiles                
def makePickle(GasDir,GasFiles,fname,fname2,fname3,fname4):
        import cPickle
        from scipy.io import loadmat
        from calendar import timegm
        GasDict={}
        Experiment={}
        ll=[]
        expIdx=0
        pIdx=0
        ii=0 
        ExperimentAxis={}
        GasAxDict={}
        TimeDict={}
        S11Dict={}
        S12Dict={}
        S22Dict={}
        ExperimentTime={}
        ExperimentS11={}
        ExperimentS12={}
        ExperimentS22={}
        Ax=[]
        timez=[]
        mm=[]
        nn=[]
        oo=[]       
        for filez in GasFiles:
                print GasDir[ii]+'/'+'raw/'+filez
                Gas=loadmat(GasDir[ii]+'/'+'raw/'+filez)
                yr=int(Gas['timestamp'][0][0:4])
                mon=int(Gas['timestamp'][0][4:6])
                day=int(Gas['timestamp'][0][6:8])
                hour=int(Gas['timestamp'][0][8:10])
                minute=int(Gas['timestamp'][0][10:12])
                second=int(Gas['timestamp'][0][12:14])
                timez.append(timegm((yr,mon,day,hour,minute,second)))
                ll.append(Gas['S21'])
                Ax.append(Gas['faxis'])
                mm.append(Gas['S11'])
                nn.append(Gas['S12'])
                oo.append(Gas['S22'])
                if(len(ll)==12):
                        pIdx+=1
                        ExperimentAxis[pIdx]=Ax[:]
                        Experiment[pIdx]=ll[:]
                        ExperimentTime[pIdx]=timez[:]
                        ExperimentS11[pIdx]=mm[:]
                        ExperimentS12[pIdx]=nn[:]
                        ExperimentS22[pIdx]=oo[:]
                        ll=[]
                        Ax=[]
                        timez=[]
                        mm=[]
                        nn=[]
                        oo=[]        
                try:
                        if(GasDir[ii]!=GasDir[ii+1]):
                                expIdx+=1
                                GasDict[expIdx]=Experiment
                                GasAxDict[expIdx]=ExperimentAxis
                                TimeDict[expIdx]=ExperimentTime
                                S11Dict[expIdx]=ExperimentS11
                                S12Dict[expIdx]=ExperimentS12
                                S22Dict[expIdx]=ExperimentS22
                                Experiment={}
                                ExperimentAxis={}
                                ExperimentTime={}
                                ExperimentS11={}
                                ExperimentS12={}
                                ExperimentS22={}
                                pIdx=0
                except:
                        expIdx+=1
                        GasDict[expIdx]=Experiment
                        GasAxDict[expIdx]=ExperimentAxis
                        TimeDict[expIdx]=ExperimentTime
                        S11Dict[expIdx]=ExperimentS11
                        S12Dict[expIdx]=ExperimentS12
                        S22Dict[expIdx]=ExperimentS22
                        ExperimentTime={}
                        Experiment={}
                        ExperimentAxis={}
                        ExperimentS11={}
                        ExperimentS12={}
                        ExperimentS22={}                
                ii+=1
                
        file_handle=open(fname,'w')
        cPickle.dump(GasDict,file_handle)
        file_handle.close()
        file_handle=open(fname2,'w')
        cPickle.dump(GasAxDict,file_handle)
        file_handle.close()
        file_handle=open(fname3,'w')
        cPickle.dump(TimeDict,file_handle)
        file_handle.close()
        file_handle=open(fname4,'w')
        cPickle.dump(S11Dict,file_handle)
        file_handle.close()
        file_handle=open(fname+'S12.pkl','w')
        cPickle.dump(S12Dict,file_handle)
        file_handle.close()
        file_handle=open(fname+'S22.pkl','w')
        cPickle.dump(S22Dict,file_handle)
        file_handle.close()

                              
def makePickleTrans(TransDir,TransFiles,fname):
        import cPickle
        from scipy.io import loadmat
        TransDict={}
        TransDictS12={}
        Experiment={}
        ExperimentS12={}
        ExperimentS11={}
        ll=[]
        rr=[]
        ss=[]
        expIdx=0
        pIdx=0
        ii=0        
        for filez in TransFiles:
        
                
                Trans=loadmat(TransDir[ii]+'/'+'raw/'+filez)
                ffname=filez.split('_')
                tNumber=int(ffname[1][5])
                if(tNumber==2):
                        print TransDir[ii]+'/'+'raw/'+filez,ii,expIdx,pIdx
                        ll.append(Trans['S21'])
                        ss.append(Trans['S12'])
                #rr.append(Trans['S11'])
                if(len(ll)==12):
                        pIdx+=1
                        Experiment[pIdx]=ll[:]
                        ExperimentS12[pIdx]=ss[:]
                        ll=[]
                        ss=[]
                                
                try:
                        if(TransDir[ii]!=TransDir[ii+1]):
                                expIdx+=1
                                TransDict[expIdx]=Experiment
                                TransDictS12[expIdx]=ExperimentS12
                                Experiment={}
                                ExperimentS12={}
                                pIdx=0
                except:
                        expIdx+=1
                        TransDict[expIdx]=Experiment
                        TransDictS12[expIdx]=ExperimentS12
                        Experiment={}
                        ExperimentS12={}                
                ii+=1
        print fname        
        file_handle=open(fname,'w')
        cPickle.dump(TransDict,file_handle)
        file_handle.close()
        file_handle=open(fname+'S12.pkl','w')
        cPickle.dump(TransDictS12,file_handle)
        file_handle.close()
if __name__ == "__main__":
    main()





     
        
                                 
