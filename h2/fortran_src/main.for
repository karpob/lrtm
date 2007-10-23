      SUBROUTINE mexFunction(nlhs, plhs, nrhs, prhs)
      integer plhs(*),prhs(*)
      integer nlhs, nrhs
      integer MXGETPR,MXCREATEDOUBLEMATRIXN
      integer temp,freq,alfatot,dnumin,dnumax,dnu
      integer M,N,NF,NNF
      REAL*8 Rtemp,Rdnumin,Rdnumax,Rdnu,Rfreq(601),Ralfatot(601)
      
c      IF (NRHS.NE.4) THEN
c            CALL MEXERRMSGTXT('I need 4 Arguements, Jerk.')
c      ELSEIF(NLHS.GT.1) THEN
c            CALL MEXERRMSGTXT('I output 1 value, Jerk.')
     
      temp=MXGETPR(PRHS(1))
      dnumin=MXGETPR(PRHS(2))
      dnumax=MXGETPR(PRHS(3))
      dnu=MXGETPR(PRHS(4))
      
      PLHS(1) = MXCREATEDOUBLEMATRIX(601,1,0)
      PLHS(2) = MXCREATEDOUBLEMATRIX(601,1,0)
           
     

      alfatot=MXGETPR(PLHS(1))
      freq=MXGETPR(PLHS(2))
     
      CALL MXCOPYPTRTOREAL8(temp, Rtemp, 1)
      CALL MXCOPYPTRTOREAL8(dnumin,Rdnumin,1)
      CALL MXCOPYPTRTOREAL8(dnumax,Rdnumax,1)
      CALL MXCOPYPTRTOREAL8(dnu,Rdnu,1)
     
      CALL ADDEM(Rtemp,Rdnumin,Rdnumax,Rdnu,Rfreq,Ralfatot)
      CALL MXCOPYREAL8TOPTR(Ralfatot,alfatot,601)
      CALL MXCOPYREAL8TOPTR(Rfreq,freq,601)
      RETURN
      END
