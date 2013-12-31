from pylab import *
from scipy.io import loadmat
import cPickle
file_handle=open('all_data.pkl','r')
Experiment=cPickle.load(file_handle)
print Experiment
    #T_vector_in_Experiment=Experiment[i]['T']
    #max_length_of_T=max([len(T_vector_in_Experiment),max_length_of_T])
file_handle.close()
