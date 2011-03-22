#include "model.h"
extern struct ATM_LAYER *layer;
extern float P_targ, T_targ, P_term, g0, R0, P0, zP0, supersatNH3, supersatH2S;
extern char use_lindal, line[];
extern int sol_cloud, NH4SH_cloud_base, AutoStep, CrossP0,jj;
extern double XCO,AutoStep_constant;
extern int use_dz;
static int n_lindal_pts,heronout;
static float Hydrogen_Curve_Fit_Select;
