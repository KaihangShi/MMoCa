! ==========================================================
! This subroutine is used to perform a MC trial move
! Created on 12-29-2016 by Kaihang Shi
! ==========================================================

	  Subroutine trial_move(iblock, istep)

	  Use global

	  IMPLICIT NONE

	  ! Passed
	  Integer :: iblock, istep

	  ! Local
	  Double Precision :: ran

	  ! Generate a random number
	  ran = random(idum)

	  ! Choose ensemble
	  Select Case (ensmbl)
	  	
	  	! NVT ensemble
	  	Case(ENS_NVT)
	  		
	  		! Choose move types
	  		If (ran .lt. move_type_prob(TRANSLATION)) Then
	  			! Translational move
	  			Call translate

	  		Else if(ran .lt. move_type_prob(ROTATION)) Then
	  			! Rotational move
	  			Call rotate

	  		Else 
	  			Write(*,*) 'FATAL ERROR: RANDOM NUMBER = 1'
	  			STOP
	  		End If


	  	! NPT ensemble
	  	Case(ENS_NPT)

	  		! Choose move types
	  		If (ran .lt. move_type_prob(TRANSLATION)) Then
	  			! Translational move
	  			Call translate

	  		Else if(ran .lt. move_type_prob(ROTATION)) Then
	  			! Rotational move
	  			Call rotate

	  		Else if(ran .lt. move_type_prob(VOLCHANGE)) Then
	  			! Volume change, only for box 1
	  			Call volchange_npt(1)
	  		Else
	  			Write(*,*) 'FATAL ERROR: RANDOM NUMBER = 1'
	  			STOP
	  		End If

	  	! uVT ensemble (Grand Canonical)
	  	Case(ENS_uVT)

	  		! Choose move types
	  		If (ran .lt. move_type_prob(TRANSLATION)) Then
	  			! Translational move
	  			Call translate

	  		Else if(ran .lt. move_type_prob(ROTATION)) Then
	  			! Rotational move
	  			Call rotate

	  		Else if(ran .lt. move_type_prob(TRANSFER)) Then
	  			
	  			! Insertion/deletion move
	  			If (random(idum) .lt. 0.5d0 ) Then
	  				! Insert a molecule
	  				Call insert(1)
	  			Else 
	  				! remove a molecule
	  				Call remove(1)
	  			End If

	  		Else
	  			Write(*,*) 'FATAL ERROR: RANDOM NUMBER = 1'
	  			STOP
	  		End If

	  		
	  	Case default
	  		Write(*,*) 'FATAL ERROR: INVALID ENSEMBLE TYPE IN TRIAL_MOVE SUBROUTINE'
	  		STOP
	 
	  End Select




	  
	  	
	  
	  	Return
	  
	  End Subroutine 