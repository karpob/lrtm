def residual_enthalpy_and_PVT(p,Data_Enthalpy_H,Data_Enthalpy_PkPa,Data_Enthalpy_T,Data_Enthalpy_x_h2,Data_Enthalpy_x_h2o,Data_pXT_P,Data_pXT_T,Data_pXT_density,Data_pXT_x_h2,Data_pXT_x_h2o,Data_B_B_vals,Data_B_T_vals,Data_B_x_h2_vals,Data_B_x_h2o_vals,Data_C_vals):
        from residual_enthalpy_data import residual_enthalpy_data
        from residual_pXT_data import residual_pXT_data
        from residual_B_data import residual_B_data
        from residual_C_data import residual_C_data
        from numpy import zeros,shape
        err1=residual_pXT_data(p,Data_pXT_P,Data_pXT_T,Data_pXT_density,Data_pXT_x_h2,Data_pXT_x_h2o)
        err2=residual_enthalpy_data(p,Data_Enthalpy_H,Data_Enthalpy_PkPa,Data_Enthalpy_T,Data_Enthalpy_x_h2,Data_Enthalpy_x_h2o)
        err3=residual_B_data(p,Data_B_B_vals,Data_B_T_vals,Data_B_x_h2_vals,Data_B_x_h2o_vals)
        err4=residual_C_data(p,Data_C_vals,Data_B_T_vals,Data_B_x_h2_vals,Data_B_x_h2o_vals)
        [m1]=shape(err1)
        [m2]=shape(err2)
        [m3]=shape(err3)
        [m4]=shape(err4)
        err_total=zeros(int(m1+m2+m3+m4))
        #print m1,m2
        #print shape(err1)
        #print shape(Data_Enthalpy_H)
        #print shape(err_total[0:m1])
        err2[len(err2)-1]=err2[len(err2)-1]*10
        err2[len(err2)-6]=err2[len(err2)-6]*10
        err2[len(err2)-4]=err2[len(err2)-4]*10
        err2[len(err2)-9]=err2[len(err2)-9]*10
        err_total[0:m1]=100.0*(err1/Data_pXT_P)
        err_total[m1:m1+m2]=100.0*(err2/Data_Enthalpy_H)
        err_total[m1+m2:m1+m2+m3]=100.0*(err3/Data_B_B_vals)
        err_total[m1+m2+m3:m1+m2+m3+m4]=100.0*(err4/Data_C_vals)
        print 'error', err_total
        print 'params',p
        return err_total
