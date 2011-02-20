def ABWR (b,t,rho,rho_c):
        from pylab import zeros
        from math import sqrt,exp
        an=zeros(16)
#
#  compute residual Helmholtz free energy for MBWR equation of state
#
#c  based on Younglove & McLinden (1994), JPCRD 23:731-779
#c  equations A4 and B6
#c
#c  inputs:
#c    icomp--pointer specifying component (1..nc)
#c        t--temperature (K)
#c      rho--molar density (mol/L)
#c  output (as function value):
#c       Ar--(A-A0) (J/mol)
#c
#c  written by M. McLinden, NIST Thermophysics Division, Boulder, Colorado
#c  10-06-94  MM, original version
#c  07-20-95  MM, separate polynomial Cp0 routines from core_BWR
#c                and restructure coefficient arrays
#c  11-29-95  MM, variable lower limit on coefficient/constant arrays
#c                to accommodate ECS reference fluid
#c  02-27-96  MM, parameter n0=-ncmax to accommodate ECS-thermo model
#c
#      implicit double precision (a-h,o-z)
#      implicit integer (i-k,m,n)
#      implicit logical (l)
#c
#      parameter (ncmax=20)        !max number of components in mixture
#      parameter (nrefmx=10)       !max number of fluids for transport ECS
#      parameter (n0=-ncmax-nrefmx,nx=ncmax)
#      common /WCFBWR/ b(n0:nx,32),
#     &                pc(n0:nx),rhoc(n0:nx),tc(n0:nx),Rbwr(n0:nx),
#     &                ptr(n0:nx),rhotr(n0:nx),ttr(n0:nx),gamma(n0:nx),
#     &                tmax(n0:nx),pmax(n0:nx)
#      common /Gcnst/ R,tz(n0:nx),rhoz(n0:nx)
#      dimension an(15)
#c
#      j=icomp              !pointer to appropriate fluid
#      ABWR=0.d0
#      if (t.le.0.d0) RETURN
#c
#c  calculate a(n) terms in MBWR (Eq B2)
#c
        tinv=1.0/t
        tinv2=tinv*tinv
        tinv3=tinv2*tinv
        tinv4=tinv3*tinv
        #an[1]=Rbwr*t !not used!
        an[2]= b[1-1]*t+b[2-1]*sqrt(t)+b[3-1]+b[4-1]*tinv+b[5-1]*tinv2
        an[3]= b[6-1]*t+b[7-1]+b[8-1]*tinv+b[9-1]*tinv2
        an[4]=b[10-1]*t+b[11-1]+b[12-1]*tinv
        an[5]=b[13-1]
        an[6]=b[14-1]*tinv+b[15-1]*tinv2
        an[7]=b[16-1]*tinv
        an[8]=b[17-1]*tinv+b[18-1]*tinv2
        an[9]=b[19-1]*tinv2
        an[10]=b[20-1]*tinv2+b[21-1]*tinv3
        an[11]=b[22-1]*tinv2+b[23-1]*tinv4
        an[12]=b[24-1]*tinv2+b[25-1]*tinv3
        an[13]=b[26-1]*tinv2+b[27-1]*tinv4
        an[14]=b[28-1]*tinv2+b[29-1]*tinv3
        an[15]=b[30-1]*tinv2+b[31-1]*tinv3+b[32-1]*tinv4
#c  summation of terms 2-9 in Eq B6
        ar=0.0
        rhon=1.0
        for i in range(2,10):
                rhon=pow(rho,float(i)-1)
         # 	print i      
                ar=ar+(an[i]*rhon)/(float(i-1))    
          
#c  summation of term 10 [first exponential term)
        delsq=pow((rho/rho_c),2)
        rhoc2=pow(rho_c,2)
        rhocn=rhoc2
        expdel=exp(-delsq)
#	print 'what I want',-0.5*an[10]*rhocn*(expdel-1.0)
        ar=ar-0.5*an[10]*rhocn*(expdel-1.0)
	#print 10
#c  summation of terms 11-15 [remaining exponential terms) in Eq B7
        expmul=1.0
        expsub=1.0
        delnew=1.0
        for i in range(11,16):
	#	print i                
                rhocn=rhocn*rhoc2
                delnew=delnew*delsq
                expmul=delnew+expmul*float(i-10)
                expsub=expsub*float(i-10)
                ar=ar-0.5*an[i]*rhocn*(expdel*expmul-expsub)
        return ar*1000



#c
#c  convert from L-bar/mol to J/mol
#      ABWR=ar*[100.0d0)
#c     ABWR=ar*[R/Rbwr(j))
#c
#      RETURN
#      end                                                 !function ABWR
