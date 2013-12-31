#standard libs
import matplotlib
matplotlib.rc('text', usetex = True)
matplotlib.use('Agg')
import cPickle
from pylab import *
import pylab
import numpy
from scipy.io import savemat
import calendar
from scipy import polyval,polyfit
from scipy.optimize import leastsq

#custom libs

from python_mods.gaugeError import gaugeError
from python_mods.calculateDensityValuesAndError import calculateDensityValuesAndError 
from python_mods.plot_experiment import plot_experiment
from python_mods.plot_refractivity import plot_refractivity
from python_mods.costFunctionRefractivity import costFunctionRefractivity

from python_mods.vacuumTemperatureCost import vacuumTemperatureCost
from python_mods.new_h2o_p import new_h2o_p,costH2O,jacobianTemperatureDepH2O
from python_mods.saveNetCDF import saveNetCDF
from python_mods.collectExperimentData import collectExperimentData
from python_mods.saveExcelFile import saveExcelFile
from python_mods.calculateArgonDensity import calculateArgonDensity
from python_mods.plot_refractivity_argon import plot_refractivity_argon
from python_mods.fitWaterSelfContinuum import fitWaterSelfContinuum
from python_mods.fitWaterHeContinuum import fitWaterHeContinuum
from python_mods.fitWaterH2Continuum import fitWaterH2Continuum
olderr = numpy.seterr(all='ignore') #hide div by zeros from compressibility code.


[experiments,Experiment,
h2_density,he_density,h2o_density,
max_length_of_T,h2_density_high_P,he_density_high_P,h2o_density_high_P,
h2_density_low_P,he_density_low_P,h2o_density_low_P,
h2_density_low_T,he_density_low_T,h2o_density_low_T,
h2_density_high_T,he_density_high_T,h2o_density_high_T,
Worksheet_names,All_Garbage_Flags,Masked_Data]           =          collectExperimentData()


#For each experiment calculate the uncertainties associated with the instrumentation

Sigma_T,Sigma_P,Sigma_e_h2o=gaugeError(experiments,Experiment,max_length_of_T)

[Experiment,h2_density,he_density,h2o_density,
h2_density_high_P,he_density_high_P,h2o_density_high_P,
h2_density_low_P,he_density_low_P,h2o_density_low_P,
h2_density_low_T,he_density_low_T,h2o_density_low_T,
h2_density_high_T,he_density_high_T,h2o_density_high_T]= calculateDensityValuesAndError(experiments,Experiment,h2_density,he_density,h2o_density,
                                                                                        h2_density_high_P,he_density_high_P,h2o_density_high_P,
                                                                                        h2_density_low_P,he_density_low_P,h2o_density_low_P,
                                                                                        h2_density_low_T,he_density_low_T,h2o_density_low_T,
                                                                                        h2_density_high_T,he_density_high_T,h2o_density_high_T,
                                                                                        Sigma_T,Sigma_P)

ArDensity=numpy.zeros([h2_density.shape[0],h2_density.shape[1]])
for i in experiments:
        ArDensity[i,0:len(Experiment[i]['Pave'])]=calculateArgonDensity(Experiment[i])
                                                                                    

#save excel
Masked_Data=saveExcelFile('h2o_data.xls',experiments,Experiment,Worksheet_names,
                          Sigma_T,Sigma_P,Sigma_e_h2o,h2_density,he_density,h2o_density,
                          All_Garbage_Flags,Masked_Data)
#create a netCDF file
saveNetCDF(experiments,Masked_Data,'h2o_data.nc')

powerWater,coefWater=fitWaterSelfContinuum(Masked_Data,h2o_density,h2_density,he_density,experiments)

coefHe=fitWaterHeContinuum(Masked_Data,h2o_density,h2_density,he_density,experiments,coefWater,powerWater)
coefH2=fitWaterH2Continuum(Masked_Data,h2o_density,h2_density,he_density,experiments,coefWater,powerWater)

print 'my fitted params', coefWater,powerWater,coefH2,coefHe
constants=numpy.asarray([coefWater,powerWater,coefH2,coefHe])
fitType='none'


p=constants#optimized_p

Max_deboer=zeros([19])
Min_deboer=zeros([19])
Mean_deboer=zeros([19])

Max_goodman=zeros([19])
Min_goodman=zeros([19])
Mean_goodman=zeros([19])

Max_rosenkranz=zeros([19])
Min_rosenkranz=zeros([19])
Mean_rosenkranz=zeros([19])

Max_rosenkranz_old=zeros([19])
Min_rosenkranz_old=zeros([19])
Mean_rosenkranz_old=zeros([19])

pp=[]

		          
#print optimized_pp
optimized_pp=numpy.asarray([1.0,1.0])


#experiments=[8]
Plot_Number=0
T_val_index=0
#Fit H2O continuum with pure H2O experiments
pure_h2o_abs=zeros([12,19])
pure_h2o_freq=zeros([12,19])
pure_h2o_T=zeros([19])
pure_h2o_P=zeros([19])
N_deboer=zeros([19])
N_goodman=zeros([19])
N_rosenkranz=zeros([19])
N_rosenkranz_old=zeros([19])
N_all=zeros([19])

#experiments=[8,12]
for i in experiments:
        
        
        [Plot_Number,N_deboer[i],N_goodman[i],N_rosenkranz[i],N_rosenkranz_old[i],
        N_all[i],Max_deboer[i],Min_deboer[i],Mean_deboer[i],
        Max_goodman[i],Min_goodman[i],Mean_goodman[i],
        Max_rosenkranz[i],Min_rosenkranz[i],Mean_rosenkranz[i],
        Max_rosenkranz_old[i],Min_rosenkranz_old[i],Mean_rosenkranz_old[i]]=plot_experiment(i,Plot_Number,Masked_Data[i],
                                                                                h2o_density[i,:],h2_density[i,:],he_density[i,:],
                                                                                h2o_density_high_P[i,:],h2_density_high_P[i,:],he_density_high_P[i,:],
										h2o_density_low_P[i,:],h2_density_low_P[i,:],he_density_low_P[i,:],
 										h2o_density_high_T[i,:],h2_density_high_T[i,:],he_density_high_T[i,:],
										h2o_density_low_T[i,:],h2_density_low_T[i,:],he_density_low_T[i,:],
										Sigma_T,Sigma_P,p,constants,fitType,'pdf')

                 
        plot_refractivity(i,Plot_Number,Masked_Data[i],h2o_density[i,:],h2_density[i,:],he_density[i,:],
                          h2o_density_high_P[i,:],h2_density_high_P[i,:],he_density_high_P[i,:],
       			  h2o_density_low_P[i,:],h2_density_low_P[i,:],he_density_low_P[i,:],
        		  h2o_density_high_T[i,:],h2_density_high_T[i,:],he_density_high_T[i,:],
       		          h2o_density_low_T[i,:],h2_density_low_T[i,:],he_density_low_T[i,:],
       		          Sigma_T,Sigma_P,optimized_pp,'pdf')
       	#print ArDensity[i,:]	          
       	#plot_refractivity_argon(i,Masked_Data[i],ArDensity[i,:],'pdf')	          
        
		          																										         
   
Total_deboer=N_deboer.sum()
Total_goodman=N_goodman.sum()
Total_rosenkranz=N_rosenkranz.sum()
Total_rosenkranz_old=N_rosenkranz_old.sum()

Total_all=N_all.sum()
Max_DeBoer=Max_deboer.max()
Min_DeBoer=Min_deboer.min()
Mean_DeBoer=Mean_deboer.mean()

Max_Goodman=Max_goodman.max()
Min_Goodman=Min_goodman.min()
Mean_Goodman=Mean_goodman.mean()

Max_Rosenkranz=Max_rosenkranz.max()
Min_Rosenkranz=Min_rosenkranz.min()
Mean_Rosenkranz=Mean_rosenkranz.mean()

Max_Rosenkranz_old=Max_rosenkranz_old.max()
Min_Rosenkranz_old=Min_rosenkranz_old.min()
Mean_Rosenkranz_old=Mean_rosenkranz_old.mean()
print 'Karpowicz Steffes, 2011',Total_rosenkranz_old,'&',Max_Rosenkranz_old,'&',Min_Rosenkranz_old,'&',Mean_Rosenkranz_old,"\\"
print 'This work &',Total_rosenkranz,'&',Max_Rosenkranz,'&',Min_Rosenkranz,'&',Mean_Rosenkranz,"\\"
print '\citet{DeBoer_thesis} &',Total_deboer,'&',Max_DeBoer,'&',Min_DeBoer,'&',Mean_DeBoer,"\\"
print '\citet{Goodman_a} &',Total_goodman,'&',Max_Goodman,'&',Min_Goodman,'&',Mean_Goodman,"\\"
print 'Total &',Total_all,"&  &  &  \\"

#prestine!
count=[         True, #1
                True, #2
                True, #3
                True, #4
                True, #5
               True, #6
               False,#7
               False, #8
               False, #9
               False, #10
               False, #11
               False, #12
               True, #13
               True, #14
               True, #15
               True, #16
               False, #17
               True, #18
               False] #19
prestine_index=[0,
                 1,
                 2,
                 3,
                 4,
                 5,
                 12,
                 13,
                 14,
                 15,
                 17]
Total_deboer=N_deboer[prestine_index].sum()
Total_goodman=N_goodman[prestine_index].sum()
Total_rosenkranz=N_rosenkranz[prestine_index].sum()
Total_rosenkranz_old=N_rosenkranz_old[prestine_index].sum()

Total_all=N_all[prestine_index].sum()
Max_DeBoer=Max_deboer[prestine_index].max()
Min_DeBoer=Min_deboer[prestine_index].min()
Mean_DeBoer=Mean_deboer[prestine_index].mean()

Max_Goodman=Max_goodman[prestine_index].max()
Min_Goodman=Min_goodman[prestine_index].min()
Mean_Goodman=Mean_goodman[prestine_index].mean()

Max_Rosenkranz=Max_rosenkranz[prestine_index].max()
Min_Rosenkranz=Min_rosenkranz[prestine_index].min()
Mean_Rosenkranz=Mean_rosenkranz[prestine_index].mean()
print 'Prestine Data set'
print 'Karpowicz Steffes, 2011',Total_rosenkranz_old,'&',Max_Rosenkranz_old,'&',Min_Rosenkranz_old,'&',Mean_Rosenkranz_old,"\\"
print 'This work &',Total_rosenkranz,'&',Max_Rosenkranz,'&',Min_Rosenkranz,'&',Mean_Rosenkranz,"\\"
print '\citet{DeBoer_thesis} &',Total_deboer,'&',Max_DeBoer,'&',Min_DeBoer,'&',Mean_DeBoer,"\\"
print '\citet{Goodman_a} &',Total_goodman,'&',Max_Goodman,'&',Min_Goodman,'&',Mean_Goodman,"\\"
print 'Total &',Total_all,"&  &  &  \\"



# Create a pickle with add data+uncertainties
output_file=open('all_processed_data_including_uncertainties.pkl','w')
cPickle.dump(Experiment,output_file)
cPickle.dump(Sigma_T,output_file)
cPickle.dump(Sigma_P,output_file)
cPickle.dump(Sigma_e_h2o,output_file)
output_file.close()
print constants
print optimized_pp

