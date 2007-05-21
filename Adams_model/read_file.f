      IMPLICIT NONE
      REAL*8 fun(15,486)
      INTEGER J,K
      CHARACTER*7 crap 
      open(unit=1,file='mixrat.jup_test',status='unknown')
      read(1,*)crap
      read(1,*)crap
      read(1,*)crap
      read(1,*)crap
      
      read(1,*)crap
      read(1,*)crap
c      read(1,*)crap
      
      DO J=1,486
            
            read(1,*)fun(1,J),crap,fun(3,J),fun(4,J),fun(5,J),fun(6,J),
     ;fun(7,J),fun(8,J),fun(9,J),fun(10,J),fun(11,J),fun(12,J),
     ;fun(13,J),fun(14,J),fun(15,J)
  
      ENDDO
      close(1)
      
      DO J=1,486
            DO K=1,1
                  print *, fun(K,J)
            enddo
      enddo
      stop
      end
      
                  
