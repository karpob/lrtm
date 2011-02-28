//TCM.i
%module TCM
%{
#include <Python.h>
#include "model.h"

extern long getCloudFlags(int i);
extern double getFloatValues(int i,int j);
      
extern int intoTheVoid( float dz,
                  double XHe_i,
                  double XH2S_i,
                  double XNH3_i,
                  double XH2O_i,
                  double XCH4_i,
                  double XPH3_i,
                  double XCO,
                  float P_temp,
                  float T_temp,
                  float g0_i,
                  float R0_i,
                  float P0_i,
                  float T_targ_i,
                  float P_targ_i,
                  float P_term_i,
                  float use_lindal_in,
                  int n_lindal_pts_i,
                  float SuperSatSelf1_i,
                  float SuperSatSelf2_i,
                  float SuperSatSelf3_i,
                  float SuperSatSelf4_i,
                  float supersatNH3_i,
                  float supersatH2S_i,
                  int AutoStep_constant,
                  float Hydrogen_Curve_Fit_Select,
                  float dP_init,
                  float dP_fine,
                  float P_fine_start,
                  float P_fine_stop,
                  int use_dz,
                  float frain,
                  float select_ackerman);

%}
extern int intoTheVoid( float dz,
                  double XHe_i,
                  double XH2S_i,
                  double XNH3_i,
                  double XH2O_i,
                  double XCH4_i,
                  double XPH3_i,
                  double XCO,
                  float P_temp,
                  float T_temp,
                  float g0_i,
                  float R0_i,
                  float P0_i,
                  float T_targ_i,
                  float P_targ_i,
                  float P_term_i,
                  float use_lindal_in,
                  int n_lindal_pts_i,
                  float SuperSatSelf1_i,
                  float SuperSatSelf2_i,
                  float SuperSatSelf3_i,
                  float SuperSatSelf4_i,
                  float supersatNH3_i,
                  float supersatH2S_i,
                  int AutoStep_constant,
                  float Hydrogen_Curve_Fit_Select,
                  float dP_init,
                  float dP_fine,
                  float P_fine_start,
                  float P_fine_stop,
                  int use_dz,
                  float frain,
                  float select_ackerman);
extern long getCloudFlags(int i);
extern double getFloatValues(int i,int j);
