def main():
        
        
        dirList=[
                 #'20_500K_ 2bH2O_14bHe',
                 '21_500K_1.5bH2O_6bHe'
                 ]
        for dirz in dirList:
                changeFileNames(dirz)
        
       

def changeFileNames(dirName):
        import os
        import shutil
        pwd=os.getcwd()
        p=os.path.join(pwd,dirName,'raw')        
        matlabFilez=os.listdir(p)
        
        
        
        
        for n in matlabFilez:
            a=n.split('b')
            prefix=a[0]
            b=n.split('_')
            pressureLevel=b[len(b)-2]
            lp=len(pressureLevel)
            if('H2O' in n):
                    nameVal=a[0][0:len(a[0])-lp]+'_'+'pure.mat'
                    
            elif('vac' in n):
                    nameVal=n
            else:        
                    nameVal=prefix+'bars.mat'
            path1=os.path.join(pwd,dirName,'raw',n)
            path2=os.path.join(pwd,dirName,'raw',nameVal)
            
            
            shutil.move(path1,path2)        
                   
        
if __name__ == "__main__":
    main()
