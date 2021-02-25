REAL FUNCTION random(idum)

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! This is the random number generator taken from                  !
 ! "Numerical Recipes"                                             !
 ! generate uniform random number between [0,1)                    !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 ! Passed
 INTEGER idum
 !REAL*8 mbig, mseed, mz, fac

 ! Local
 !PARAMETER (mbig=1.0e12, mseed=1.61803398e11,mz=0.0d0,fac=1.0/mbig)
 REAL, PARAMETER :: mbig=4000000., mseed=1618033., mz=0., fac=1./mbig
 INTEGER i, iff, ii, inext, inextp, k
 REAL*8  mj, mk, ma(55)
 SAVE inext, inextp, ma, iff
 DATA iff /0/


  IF(idum.LT.0.OR.iff.EQ.0) THEN

    iff=1
    mj=ABS(mseed-ABS(idum))
    mj=MOD(mj,mbig)
    ma(55)=mj
    mk=1

    DO  i=1,54
      ii=MOD(21*i,55)
      ma(ii)=mk
      mk=mj-mk
      IF(mk.LT.mz)mk=mk+mbig
      mj=ma(ii)
    END DO

    DO  k=1,4
      DO  i=1,55
        ma(i)=ma(i)-ma(1+MOD(i+30,55))  
        IF(ma(i).LT.mz)ma(i)=ma(i)+mbig
      END DO
    END DO

    inext=0
    inextp=31
    idum=1

  ENDIF
      
  inext=inext+1
  IF(inext.EQ.56)inext=1
  inextp=inextp+1
  IF(inextp.EQ.56)inextp=1
  mj=ma(inext)-ma(inextp)
  IF(mj.LT.mz)mj=mj+mbig
  ma(inext)=mj
  random=mj*fac


 RETURN

 END FUNCTION random
