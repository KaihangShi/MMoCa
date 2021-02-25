! ==========================================================
! This subroutine is used to gradually change chemical potential
! of uVT system for dense system. 
! This algorithm is recommended by several paper
!
! Created on 8-13-2017 by Kaihang Shi
! Ref.
! [1] Kowalczyk, Piotr, et al. "Thermodynamics of hydrogen adsorption in slit-like 
! carbon nanopores at 77 K. Classical versus path-integral Monte Carlo simulations." Langmuir 23.7 (2007): 3666-3672.
! 
! !!! Read before Proceed!!!
! Right now, the starting chempot is -2084.716617d0 (Argon@87.3K, 2e-5bar)
! Change it for other specific system
! ==========================================================



	  Subroutine update_chempot(iblock)
	  
	  Use global

	  IMPLICIT NONE

	  ! Passed
	  Integer :: iblock

	  ! Local
	  Integer :: nchange, itype


	  ! Return if not uVT system
	  If (ensmbl .ne. ENS_uVT) Return

	  ! Change chempot 10 times successively
	  ! Only change chempot during the first 1000 blocks
	  nchange = INT(DBLE(iblock)/100.0d0)

	  Do itype = 1, n_mol_types

	  	mu(itype) = -2084.716617d0 + DBLE(nchange)*dmu(itype)
	  	Write(*,'(A,I0,A)') 'Chemical potential for block ', iblock, ': '
	  	Write(*,'(A,I0,A,F10.4)') 'Molecule type ', itype, ': ', mu(itype)

		! Calculate activity constant * volume (GCMC)
		mol_act(itype) = vol_pore(1)*dEXP(mu(itype)/temp)/(lambda(itype)**3)
	  	
	  End Do


	  	
	  
	  	Return
	  
	  End Subroutine 













