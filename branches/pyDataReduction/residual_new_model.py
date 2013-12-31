
def residual_new_model(p,density_h2o,density_h2,density_he,Data,flags):
        from python_mods.new_h2o_p import new_h2o_p
        from numpy import zeros,shape,ma,power,max,min,abs
        abs_calc=ma.zeros([12,19,11])
        abs_meas=ma.zeros([12,19,11])
        err=ma.zeros([12,19,11])
        pure_flags=zeros([12,19,11])
        #print shape(density_h2o)
        Tval_index=0
        
        flag_repeat=0.0
        #################
        #He optimization
        ####################################
        #experiments=[4,9,11,12,14,15]
        #score=[10,1.0,1.0,1.0,1.5,1.0]
        ##################################
        ###################
        #H2 optimization
        #################
        #experiments=[5,10,13]
        #max_pressure_in_experiment=[6,6,6]
        #score=[10,1.0,1.0]
        #Final optimization 1.03562010226e-10 5.07722009423e-11 1.0 1.0


        experiments=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]    
        max_pressure_in_experiment=[1, 1,     1,  1, 1,  1,    1,     1,  1,  1, 1,  1,   1,  1, 1, 1,   1 ]#[3,7,6,7,6,6,5,5,6,6,6,6,7,6,6,6]
        score=[0.0, #1
               0.0, #2
               0.0, #3
               0.0, #4
               0.0, #5
               0.0, #6
               0.0, #7
               0.0, #8
               1.0, #9
               0.0, #10
               0.0, #11
               0.0, #12
               1.0, #13
               0.0, #14
               0.0, #15
               0.0, #16
               1.0, #17
               1.0, #18
               1.0,] #19
        #[0.0, 0.0,     0.0,  0.0,  0.0,    0.0,     0.0,  0.0,  1.8, 2.0,  1.8,   1.0,  1.0,  1.0,  1.0,   1.0 ]
        #score=       [10000,  0,     0,  0,  10000,    0,     0,  10000,  15000.0, 10000,  10000,   10000,   0,  0,   50,   0 ]
             #[       0,     1,     2,  3,     4,     5,    6,   7,   8  ,  9,      10,     11,    12,   13,  14,   15]
        #score=      [2.0, 21110.0, 5.0,21120.0,1800.0,40.0,10.0,1.0,200000.0,1.0,200000.0,30.0,30.0,1.0,420.0,1.0]
        #score=      [2.0, 110.0, 5.0,120.0,800.0,40.0,10.0,1.0,1.0,1.0,2500.0,30.0,30.0,1.0,420.0,1.0]
        #                           [0,1,2,3,4,5,6,7,8,9,0,1,2,3,4,5,6,7,8]
        #max_pressure_in_experiment=[3,7,6,7,6,6,5,5,6,6,6,6,7,6,6,6,4,6,5]
        ii=0
        for i in experiments:
                #print 'experiment included',i,
                #print 'weight',score[i]
                for j in range(0,max_pressure_in_experiment[ii]):#len(Data[i]['Pave'][Tval_index,:])):
                        
                        abs_calc[0:12,i,j]=new_h2o_p(p,Data[i]['nh3PeakFa'][:,j]/1e9,density_h2o[i,j],density_h2[i,j],density_he[i,j],max(Data[i]['T'][:,j]))
                        abs_meas[0:12,i,j]=Data[i]['Absorption'][0:12,j]
                        #err[:,i,j]=score[ii]*power(density_h2o[i,j],1.0)*100.0*(abs_calc[:,i,j]-abs_meas[:,i,j])/abs_meas[:,i,j]
                        if (i==8 and j==0):
                                ff=1.0
                        else:
                                ff=1.0        
                        err[:,i,j]=ff*score[ii]*(abs_calc[:,i,j]-abs_meas[:,i,j])*(1.0/Data[i]['Absorption_2Sigma'][0:12,j])#*score[ii]
                              
                        #        else:
                        #                err[:,i,j]=ma.masked
                                     
                ii=ii+1                         
        
        err[err==0.0]=ma.masked
        err.shape=12*19*11
        print shape(err)
        print max(abs(err)),min(abs(err))
        return err
        
                        #(i>17)|(i==8)|(i==0)
                #if((i!=3)|(i!=4)|(i!=9)|(i!=11)|(i!=12)|(i!=14)|(i!=15)|(i!=16)|(i!=17)):
                #        err[0:12,i]=ma.masked#ma.zeros(shape(abs_calc[:,i])).masked#(abs_calc[:,i]-abs_meas[:,i])
                #else:
                
