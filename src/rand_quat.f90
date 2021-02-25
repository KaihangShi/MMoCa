! ==========================================================
! This subroutine is used to generate a random orientation 
! for the molecule
! Created on 1/5/2017 by Kaihang Shi
!
! Reference:
! Vesely, Journal of Computational Physics 47.2 (1982): 291-296.
! ==========================================================

	  Subroutine rand_quat(q1i,q2i,q3i,q4i)

	  Use global

	  IMPLICIT NONE 

	  ! Passed 
	  Double Precision :: q1i, q2i, q3i, q4i

	  ! Local
	  Double Precision :: S1, S2

	  ! Initialize Local variables
	  S1 = 1.0d0
	  S2 = 1.0d0

	  ! Step 1
	  Do
	  	q1i = 2.0d0 * random(idum) - 1.0d0
	  	q2i = 2.0d0 * random(idum) - 1.0d0
	  	S1 = q1i**2 + q2i**2
	  	if(S1 .lt. 1.0d0) Exit
	  End do
	  
	  ! Step 2
	  Do
	  	q3i = 2.0d0 * random(idum) - 1.0d0
	  	q4i = 2.0d0 * random(idum) - 1.0d0
	  	S2 = q3i**2 + q4i**2
	  	if(S2 .lt. 1.0d0) Exit
	  End do

	  ! Set the random unit quaternion (Step 3)
	  q1i = q1i
	  q2i = q2i
	  q3i = q3i * dSQRT((1.0d0 - S1)/S2)
	  q4i = q4i * dSQRT((1.0d0 - S1)/S2)



	  	Return
	  
	  End Subroutine 