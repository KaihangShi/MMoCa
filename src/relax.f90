! ==========================================================
! This subroutine is used to create a relaxed system
! Using biased algorithm to get to the quasi-equilibrium faster
! For some high density system (e.g., slit pore) this might be 
! necessary
! Created on 8-12-2017 by Kaihang Shi
! Adapted from Jeremy Palmer's code
! ==========================================================

      Subroutine relax

      Use global

      IMPLICIT NONE

      ! Local 
      Integer :: ibox, iblock, istep, nstep
      Double Precision, Dimension(n_box) :: engtst

      ! Return if no relaxation blocks
      If (n_blocks_relax .eq. 0) Return

      ! Initialize local variable
      nstep = 0

      ! Print banner and progress
	  WRITE(*,'(A,I0,A)') '== Begin relaxing ',n_box,' simulation boxes... =='

	  ! Print initial energies for each box
	  WRITE(*,'(A)')'Total energies before relaxation: '
	  DO ibox = 1,n_box
	    WRITE(*,'(A,I0,A,E20.12,A)') 'Box ',ibox,': ', R*energy(ibox)/DBLE(n_mol_tot(ibox)),' J/mol'
	  ENDDO

	  ! Loop over the blocks
	  DO iblock = 1,n_blocks_relax

	  	! Flush I/O unit (write information in buffer immediately to file)
	  	! unit=6 is screen
	  	Call FLUSH(6)
	 
	    ! Loop over the steps per block
	    DO istep = 1,block_size

	    	! Update nstep
	    	nstep = nstep + 1

	    	! Recalibrate total energy every so often to mitigate error
	  		If (MOD(nstep,1000) .eq. 0) Then
	  			Do ibox = 1, n_box
	  				Call eng_total(ibox,energy(ibox))
	  			End Do
	  		End if

	  		! Select ensemble
	  		Select Case (ensmbl)

	  			! Canonical
	  			Case(ENS_NVT)

	  				Call translate
		  			Call rotate
		  		! NPT
		  		Case(ENS_NPT)

		  			Call translate
		  			Call rotate
	  			
	  			! Grand canonical ensemble 
	  			Case(ENS_uVT)

		  			Call translate
		  			Call rotate
		  			! Insertion/deletion move
		  			If (random(idum) .lt. 0.5d0 ) Then
		  				! Insert a molecule
		  				Call insert(1)
		  			Else 
		  				! remove a molecule
		  				Call remove(1)
		  			End If


		  			
	  				
	  			Case default
	  				
	  		
	  		End Select

	  	! End loop over steps per block
	    ENDDO

	    ! Print temporary results after each block for reference
	    WRITE(*,'(A,I0)') 'Relax block ', iblock
	  	DO ibox = 1,n_box
	    	WRITE(*,'(A,I0,A,I0)') 'Total number of molecules in box ',ibox,': ', n_mol_tot(ibox)
	    	WRITE(*,'(A,I0,A,E20.12,A)') 'Total energies in box ',ibox,': ', R*energy(ibox)/DBLE(n_mol_tot(ibox)),' J/mol'
	  	ENDDO
	  
	  ! End loop over iblock
	  ENDDO

	  ! Calculate energy and print configuration after relax
	  Do ibox = 1, n_box
	  	Call eng_total(ibox,energy(ibox))
	  	Call print_xyz(ibox,nstep,'sites')

	  	! Check if overlap
	  	If (OVERLAP) Then
	  		Write(*,'(A,I0)') 'FATAL ERROR: INITIAL CONFIGURATION OVERLAPS FOR BOX: ', ibox
	  		STOP
	  	End If

	  End Do

	  ! Print final energies for each box
	  WRITE(*,*) '== End relaxation stage =='


            
      	Return
      
      End Subroutine relax












