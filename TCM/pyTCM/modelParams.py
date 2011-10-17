
MAXLAYERS=5000
R=8.3143E7    #/* Universal gas constant [erg/K/mol] */
AMU=1.66056E-24 #/* Atomic Mass Unit in grams */
AMU_H2=2.01594  #   /*updated karpowicz*/
AMU_He=4.0026   #   /*updated bk*/  /* constituent masses [amu] */
AMU_H2S=34.076
AMU_NH3=17.030
AMU_H2O=18.015486068204702 #/*bk*/ 
AMU_CH4=16.0428  #/*bk*/
AMU_PH3=33.997
AMU_NH4SH=51.110

#/* triple points */
TRIPLEPT_H2S=187.61
TRIPLEPT_NH3=195.5       
TRIPLEPT_CH4=90.7
TRIPLEPT_PH3=0.0
TRIPLEPT_H2O=273.16
GONE=1e-30       #/* value when a constituent is considered gone */
ZERO=1e-30       #/* zero cloud density (for plotting) */
COUNT_CLOUD=1e-3  #     /* count cloud opacity if greater than this number */
MAXTRIES=100      #    /* Maximum number of iterations to fit deep temperature */
TLIMIT=0.05       #  /* allowed % error in T */

