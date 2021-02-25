! ==========================================================
! This subroutine is used to update probabilities of translate/
! rotate/transfer/volume moves during the simulation 
! For some high density system (e.g., slit pore) this might be 
! necessary
! Created on 8-11-2017 by Kaihang Shi
! Last modified on 8-11-2017 by Kaihang Shi
! ==========================================================


	  Subroutine update_prob(iblock)

	  Use global

	  IMPLICIT NONE

	  ! Passed
	  Integer :: iblock

	  ! Local
	  Integer :: imove



	  If (DBLE(iblock) .le. DBLE(n_blocks_equil)/2.0d0) Then

	  	! Do nothing
	  	Continue

	  Else if ((DBLE(iblock) .gt. DBLE(n_blocks_equil)/2.0d0) .and. (iblock .le. n_blocks_equil)) Then

	  	! Update the probability of moves for later half part of equilibrium
	  	Do imove = 1,4

	  		move_type_prob(imove) = move_type_prob_update(1,imove)
	  		
	  	End Do

	  Else if (iblock .gt. n_blocks_equil) Then

	  	! Update the probability of moves for the production stage
	  	Do imove = 1,4

	  		move_type_prob(imove) = move_type_prob_update(2,imove)
	  		
	  	End Do
	  	
	  End If

	  
	  	
	  
	  	Return
	  
	  End Subroutine 
















