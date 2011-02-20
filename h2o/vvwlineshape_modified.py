def vvwlineshape_modified(f,fo,df,shift):
        from numpy import ones,power,pi,abs
        # This function determines the Van,Velck,Wiesskopf (VVW) lineshape over
        # a range of freqencies which are user defined as is the frequency step
        # the resulting vector is passed to this subroutine as f.
        # fo is a vector containing the center frequencies (GHz) of the spectral 
        # lines of the given molecule (ie PH3) as supplied by the "Submillimeter,
        # Millimeter, and Microwave Spectral Catalog" aka Poynter-Pickett Line Catalog
        # but in GHz, not in MHz as provided by the catalog.
        # df is elsewhere called d_nu in the main program and is the line half-width at
        # half-maximum or simply linewidth
        # 
        # Bryan Karpowicz: Removed df term so fvvw is actually fvvw, and not
        # df*fvvw. Also added cutoff at 750 GHz for consistency with Rosenkranz.
        #
        #
        # NOTE: All frequencies are or have been converted to inverse cm!!!

        # fvvw(f,fo,df)=(1/pi)*(f/fo)^2*(    df          +         df	   )
        #                              --------------     --------------      -       base
        #                              (fo-f-shift)^2 +df^2     (fo+f+shift)^2 +df^2

        # F =                A      *(     B          +         C       )

        n=len(f)
        m=len(fo)
        # f(1) f(2) f(3) f(4) ....f(n) n times where n is the number of frequency steps
        # f(1) f(2) f(3) f(4)				in the observation range                            
        # ...
        # m times where m is the number of spectral lines

        nones=ones([1,n])
        mones=ones([m,1])

        f_matrix=mones*f
        fo_matrix=fo*nones
        df_matrix=df*nones      
        shift_mat=shift*nones
        df_minus=f_matrix-fo_matrix-shift_mat
        df_plus=f_matrix+fo_matrix+shift_mat
        
        cut_off=750

        base=df/(562500.0+power(df,2))*nones
        
        A=(1.0/pi)*power((f_matrix/fo_matrix),2);			# f(1)*fo(1)  f(2)*fo(1) f(3)*fo(1)
                                                        # f(1)*fo(2)  f(2)*fo(2) f(3)*fo(2)
        B=df_matrix/(power((df_minus),2)+power(df_matrix,2))
        C=df_matrix/(power((df_plus),2)+power(df_matrix,2))

        filter1=abs(f_matrix-fo_matrix-shift_mat)<=cut_off;
        filter2=abs(f_matrix+fo_matrix+shift_mat)<=cut_off;
        
        


        F=A*(filter1*(B-base)+filter2*(C-base))#.*df_matrix;						% m x n matrix
        
        return F
        # Note: the extra df factor comes in because of the expression for alpha has
        # an additional df factor which cancels with the alpha_max, but is added here
        # for computational reasons.  result is actually df*Fvvw

