#include <Python.h>
/* get_P_from_python->is used to interface between python_compressibility/calc_Cp/Specific_Heat_CAPI.py and the TCM.
              Inputs:
                     -->T    :temperature in deg K.
                     -->PH2  : Ideal pressure of Hydrogen in bars (density proxy)
                     -->PHe  : Ideal pressure of Helium in bars (density proxy)
                     -->PCH4 : Ideal pressure of Methane in bars (density proxy)
                     -->PH2O : Ideal pressure of Water in bars (density proxy)
              Outputs:
                     <--vals[0]: Real Pressure in bars.
                     <--vals[1]: Specific heat in erg/K/mol
*/
/* get_P_from_python->is used to interface between python_compressibility/calc_Cp/Specific_Heat_CAPI.py and the TCM.
              Inputs:
                     -->T    :temperature in deg K.
                     -->PH2  : Ideal pressure of Hydrogen in bars (density proxy)
                     -->PHe  : Ideal pressure of Helium in bars (density proxy)
                     -->PCH4 : Ideal pressure of Methane in bars (density proxy)
                     -->PH2O : Ideal pressure of Water in bars (density proxy)
              Outputs:
                     <--vals[0]: Real Pressure in bars.
                     <--vals[1]: Specific heat in erg/K/mol
*/
             
double* get_P_from_python(float T,float PH2,float PHe,float PCH4,float PH2O)
{   
    char *current_path;
    PyObject *mymod, *strfunc;
    PyObject *py_out;
    double P;
    double Cp;
    //FILE* output_T_P;
    double static vals[2];
    char *path,*newpath,*to_gaslib,*to_site_packages;

    to_gaslib="/TCM/python_compressibility/calc_Cp";
    to_site_packages="/usr/local/lib/python2.7/site-packages";// <-- make me an input arg

    current_path=getcwd(NULL,0);
    //output_T_P=fopen("python_compressibility/calc_Cp/output_T_P.txt","w");
    //fprintf(output_T_P,"%f\t%f\t%f\t%f\t%f",T,PH2,PHe,PCH4,PH2O);
    //fclose(output_T_P);
    /****************Begin Path Stuff************************/
    //equivalent to pwd
    Py_Initialize();

    //get paths that python interpreter sees 
    path=Py_GetPath();
    
    //make a blank character array the size of paths + ":"
    newpath= new char[strlen(path)+strlen(current_path)+4+strlen(to_gaslib)+4+strlen(to_site_packages)];
    
    //copy python path to newpath
    strcpy(newpath,path);
    
    //add a : separator
    strcat(newpath,":");
    
    //add current path
    strcat(newpath,current_path);
    strcat(newpath,to_gaslib);

    //add site packages
    strcat(newpath,":");
    strcat(newpath,to_site_packages);
    /*****************End Path stuff**************************/

    
    //tell python where your modules are
    //printf("path %s",newpath);

    PySys_SetPath(newpath);     
    
    //load your module (filename)
    mymod=PyImport_ImportModule("Pressure_CAPI");

    //load your function def
    strfunc=PyObject_GetAttrString(mymod,"PressureCAPI");

   //call your function with built python arg
    
    py_out=PyEval_CallFunction(strfunc,"fffff",T,PH2,PHe,PCH4,PH2O);
    
   //pull out the answer from the python function
    PyArg_Parse(py_out,"(dd)",&P,&Cp);
    Py_DECREF(py_out);
    Py_DECREF(strfunc);
    Py_DECREF(mymod);
    PySys_SetPath(current_path);
    
    
    /**************************Done Talking to Python*******************/
    vals[0]=P;
    vals[1]=Cp;
    
    return vals;    
}
