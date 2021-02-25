SUBROUTINE  seed_random(idum)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                                                                    !
  ! Initializes the seed of the random number generator                !
  !                                                                    !
  ! Updates global variables: idum                                     !
  !                                                                    !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IMPLICIT NONE

  ! Local
  INTEGER :: idum
  INTEGER, DIMENSION(8) :: val


  CALL DATE_AND_TIME(VALUES=val)
  idum=-(1000*val(7)+val(8))
!  idum = -200


  RETURN

END SUBROUTINE seed_random
