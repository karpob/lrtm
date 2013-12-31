def hist_new_model(p,density_h2o,density_h2,density_he,Data,flags):
        from python_mods.new_h2o_p import new_h2o_p
        from python_mods.goodman
        from numpy import zeros,shape,ma
        from scipy.stats import histogram
        abs_calc=ma.zeros([12,19,11])
        abs_meas=ma.zeros([12,19,11])
        err=ma.zeros([12,19,11])
        pure_flags=zeros([12,19,11])
        #print shape(density_h2o)
        Tval_index=0
        p=abs(p)
        flag_repeat=0.0
        experiments=[0,1,2,3,4,9,11,12,14,15]
        max_pressure_in_experiment=[3,7,6,7,6,6,6,7,6,6]
        #max_pressure_in_experiment=[3,7,6,7,6,6,5,5,6,6,6,6,7,6,6,6,4,6,5]
        ii=0
        for i in experiments:
                
                for j in range(1,max_pressure_in_experiment[ii]):#len(Data[i]['Pave'][Tval_index,:])):
                        
                        abs_calc[0:12,i,j]=new_h2o_p(p,Data[i]['nh3PeakFa'][:,j]/1e9,density_h2o[i,j],density_h2[i,j],density_he[i,j],max(Data[i]['T'][:,j]))
                        abs_meas[0:12,i,j]=Data[i]['Absorption'][0:12,j]
                        
                        #        else:
                        #                err[:,i,j]=ma.masked
                                     
                ii=ii+1                         
        
        err[err==0.0]=ma.masked
        err.shape=12*19*11
        abs_calc.shape=12*19*11
        #print err
        
        return err
        
                        #(i>17)|(i==8)|(i==0)
                #if((i!=3)|(i!=4)|(i!=9)|(i!=11)|(i!=12)|(i!=14)|(i!=15)|(i!=16)|(i!=17)):
                #        err[0:12,i]=ma.masked#ma.zeros(shape(abs_calc[:,i])).masked#(abs_calc[:,i]-abs_meas[:,i])
                #else:
                
