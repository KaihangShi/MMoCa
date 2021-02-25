! ==========================================================
! This subroutine is used to reset/accumulate/calculate averages
! and standard deviations of thermo statistics
! Created on 12-27-2016 by Kaihang Shi
! 
! Note: Widom insertion method is only performed during the 
!		production run. So the average chemical potential is 
!		only computed in production stage.
! ==========================================================


	  Subroutine stats(selector, iblock, istep)

	  Use global

	  IMPLICIT NONE

	  ! Passed
	  Integer :: selector, iblock, istep

	  ! Local
	  Integer :: ibox, itype, istp, iblk
	  Double Precision, Dimension(n_mol_types) :: nmolavg, rhoavg, nmolstd, rhostd
	  Double Precision :: engavg, engstd, n_total
	  Double Precision :: engdispavg, engdispstd
	  Double Precision :: ewldself, ewldintra
	  ! Define variables for widom insertion
	  Double Precision, Dimension(n_mol_types) :: muavg, mustd


	  ! Select statistics operation
	  Select Case (selector)


	  	! Reset the statistics
	  	! Input: selector
	  	Case(1)

	  		! Reset the statistics arrays
	  		! Move counter must be reset every block
	  		vol_stat = 0
		    swap_ident_stat = 0
		    trans_stat = 0
		    transfer_stat = 0
		    rotat_stat = 0


	  		! Initialize block averages
	  		blk_nmol = 0
	  		blk_rho = 0.0d0
	  		blk_eng = 0.0d0
	  		blk_eng_disp = 0.0d0
	  		blk_vol = 0.0d0

	  		! Reset Widom Insertion Counter
	  		if(lwdm) widom_stat = 0


	  	! Accumulate current block statistics
	  	! Input: istep
	  	Case(2)

	  		! Loop over the boxes
	  		Do ibox = 1, n_box

	  			! Loop over molecule type
	  			Do itype = 1, n_mol_types

	  				! Update number molecule and number density
	  				blk_nmol(istep,itype,ibox) = n_mol(itype,ibox)
	  				blk_rho(istep,itype,ibox) = DBLE(n_mol(itype,ibox))/vol(ibox)

	  			End Do

	  			! Update volume (no need this time)
!	  			blk_vol(istep,ibox) = vol(ibox)

				! Calculate the total energy per molecule [K/molecule]
				If (n_mol_tot(ibox) .gt. 0) Then
					blk_eng(istep,ibox) = energy(ibox)/DBLE(n_mol_tot(ibox))
				Else
					blk_eng(istep,ibox) = 0.0d0
				End If

				! Calculate the dispersion energy per molecule [K/molecule]
				If (n_mol_tot(ibox) .gt. 0) Then
					blk_eng_disp(istep,ibox) = energy_disp(ibox)/DBLE(n_mol_tot(ibox))
				Else
					blk_eng_disp(istep,ibox) = 0.0d0
				End If



	  		End Do

	  	
	  	! Calculate block averages 
	  	! Input: iblock
	  	Case(3)

	  		! Loop over every step in iblock
	  		Do istp = 1, block_size

	  			! Loop over each box
	  			Do ibox = 1, n_box

	  				! Loop over each molecule type
	  				Do itype = 1, n_mol_types

	  					avg_nmol(iblock,itype,ibox) = avg_nmol(iblock,itype,ibox) + &
	  												& DBLE(blk_nmol(istp,itype,ibox))
	  					avg_rho(iblock,itype,ibox) = avg_rho(iblock,itype,ibox) + &
	  												& blk_rho(istp,itype,ibox)
	  					
	  				End Do

	  				avg_eng(iblock,ibox) = avg_eng(iblock,ibox) + &
	  									& blk_eng(istp,ibox)

	  				avg_eng_disp(iblock,ibox) = avg_eng_disp(iblock,ibox) + &
	  									    & blk_eng_disp(istp,ibox)

	  			End Do
	  			
	  		End Do

	  		! Compute the average
	  		! Loop over boxes
	  		Do ibox = 1, n_box

	  			! Loop over molecule types
	  			Do itype = 1, n_mol_types
	  				
	  				avg_nmol(iblock,itype,ibox) = avg_nmol(iblock,itype,ibox)/DBLE(block_size)
	  				avg_rho(iblock,itype,ibox) = avg_rho(iblock,itype,ibox)/DBLE(block_size)

	  			End Do

	  			avg_eng(iblock,ibox) = avg_eng(iblock,ibox)/DBLE(block_size)
	  			avg_eng_disp(iblock,ibox) = avg_eng_disp(iblock,ibox)/DBLE(block_size)

	  		End Do

	  	! Print block statistics 
	  	! Input: iblock
	  	Case(4)

	  		Write(*,'(A,I0,A)') '* Steps: ', iblock*block_size, ' * '
	  		Write(*,'(A,I0)') 'Block Average for Block ', iblock


	  		! Print thermo stats
	  		! Loop over boxes
	  		Do ibox = 1, n_box

	  			! Calculate average total number of molecules in ibox for iblock
	  			! Initialze variable
	  			n_total = 0.0d0
	  			Do itype = 1, n_mol_types
	  				! Update total number of molecules for ibox
	  				n_total = n_total + avg_nmol(iblock,itype,ibox)
	  			End Do

	  			Write(*,'(A,I0)') 'Box: ', ibox
	  			Write(*,'(A,F8.3,2X,F8.3,2X,F8.3)') 'Box Length: ', box(1,ibox), box(2,ibox), box(3,ibox)
	  			Write(*,'(A,F20.4)') 'Box Volume: ', vol(ibox)
	  			Write(*,'(A,F20.4)') 'Pore Volume: ', vol_pore(ibox)
	  			Write(*,'(A,E15.7,A,4X,E15.7,A)') 'Total Energy: ', R*avg_eng(iblock,ibox), ' [J/mol]', &
	  																& n_total*avg_eng(iblock,ibox),' [K]'
!	  			Write(*,'(A,E15.7,A,4X,E15.7,A)') 'Inter vdW   : ', R*avg_eng_disp(iblock,ibox), ' [J/mol]', &
!	  																& n_total*avg_eng_disp(iblock,ibox),' [K]'

				Write(*,'(A,F10.7)') 'rcelect: ', rcelect
				Write(*,'(A,F10.7)') 'alpha: ', alpha

	  			! Loop over molecule types
	  			Do itype = 1, n_mol_types

	  				Write(*,'(A,I0,2X,A)') 'Molecule type: ', itype, mol_type_name(itype)
	  				Write(*,'(A,F10.3)') 'Molecule Number: ', avg_nmol(iblock,itype,ibox)
	  				Write(*,'(A,E15.7,A,2X,E15.7,A)') 'Number Density: ', avg_rho(iblock,itype,ibox), '[A^-3]', &
	  					& (mol_mass(itype)/(Na*1.0d-24))*avg_rho(iblock,itype,ibox), '[g/ml]'
	  				Write(*,'(A,F15.7,A)') 'Mu:              ', widom_sample(iblock,itype,ibox), ' [K]'

	  				! Write stats of translational move
	  				If (trans_stat(1,itype,ibox) .gt. 0) Then
	  					Write(*,*) ' Translational Move '
	  					Write(*,'(A,I0,2x,A,I0,2x,A,F7.3,A)') 'Attempted: ', trans_stat(1,itype,ibox), &
	  					& 'Accepted: ', trans_stat(2,itype,ibox),'Accepted Ratio: ', &
	  					& 100.0d0*DBLE(trans_stat(2,itype,ibox))/DBLE(trans_stat(1,itype,ibox)), ' %'
	  					Write(*,'(A,F6.2)') 'Max Displacement: ', max_trans(itype,ibox)
	  				End If

	  				! Write stats of rotational move
	  				If(rotat_stat(1,itype,ibox) .gt. 0) Then
	  					Write(*,*) ' Rotational Move '
	  					Write(*,'(A,I0,2x,A,I0,2x,A,F7.3,A)') 'Attempted: ', rotat_stat(1,itype,ibox), &
	  					& 'Accepted: ', rotat_stat(2,itype,ibox),'Accepted Ratio: ', &
	  					& 100.0d0*DBLE(rotat_stat(2,itype,ibox))/DBLE(rotat_stat(1,itype,ibox)), ' %'
	  					Write(*,'(A,F6.2)') 'Max Rotation: ', max_rotat(itype,ibox)
	  				End If

	  				! Write stats of volume change move
	  				If (vol_stat(1,ibox) .gt. 0) Then
	  					Write(*,*) ' Volume Change Move '
	  					Write(*,'(A,I0,2x,A,I0,2x,A,F7.3,A)') 'Attempted: ', vol_stat(1,ibox), &
	  					& 'Accepted: ', vol_stat(2,ibox),'Accepted Ratio: ', &
	  					& 100.0d0*DBLE(vol_stat(2,ibox))/DBLE(vol_stat(1,ibox)), ' %'
	  					Write(*,'(A,F6.2)') 'Max Volume Change: ', max_vol(ibox)
	  				End If

	  				! Write stats of insert/remove move
	  				If ((in_stat(1,itype) .gt. 0) .or. (rem_stat(1,itype) .gt. 0)) Then
	  					Write(*,*) ' Insertion Move '
	  					Write(*,'(A,I0,2x,A,I0,2x,A,F7.3,A)') 'Attempted: ', in_stat(1,itype), &
	  					& 'Accepted: ', in_stat(2,itype),'Accepted Ratio: ', &
	  					& 100.0d0*DBLE(in_stat(2,itype))/DBLE(in_stat(1,itype)), ' %'
	  					Write(*,*) ' Deletion Move '
	  					Write(*,'(A,I0,2x,A,I0,2x,A,F7.3,A)') 'Attempted: ', rem_stat(1,itype), &
	  					& 'Accepted: ', rem_stat(2,itype),'Accepted Ratio: ', &
	  					& 100.0d0*DBLE(rem_stat(2,itype))/DBLE(rem_stat(1,itype)), ' %'
	  				End If

	  			! End loop over molecule types
	  			End Do	
	  		! End loop over boxes
	  		End Do


	  	! Print final results 
	  	! Input: selector
	  	Case(5)

	  		! Write MC move stats
	  		Write(*,*) 
	  		Write(*,*) '+++++ End of Markov Chain +++++'
	  		Write(*,*)


	  		! Print energies of final configuration to the screen
	  		Write(*,*)
	  		! Loop over boxes
	  		Do ibox = 1, n_box
	  			! Calculate total self correction and INTRA-molecular contributions
	  			ewldself = 0.0d0 
	  			ewldintra = 0.0d0

	  			DO itype = 1,n_mol_types

	  				! Check if external structure
	    			if(initstyle(itype,ibox) .eq. 'coords') Then

	    				! Check direct Coulombic interaction for ad-surf feature
	    				If (lcoulsc) CYCLE

	    				! Exclude self/intra term of external structure
	    				If (mix_rule .eq. EXPLICIT) CYCLE

	    			End if

	    			ewldself = ewldself + n_mol(itype,ibox)*ewld_self(itype)
	    			ewldintra = ewldintra + n_mol(itype,ibox)*ewld_intra(itype)
	  			ENDDO
	  			! Concert to [K]
	  			ewldself = ewldself * EETOK
	  			ewldintra = ewldintra * EETOK


	  			Write(*,'(A,I0)') 'Final Energy for Box: ', ibox
	  			Write(*,'(A,4X,E15.7,A)') 'Total Energy                   ', energy(ibox),' [K] '
	  			Write(*,'(A,4X,E15.7,A)') 'Inter vdW                      ', energy_disp(ibox),' [K] '
	  			Write(*,'(A,4X,E15.7,A)') 'vdW Tail Correction            ', energy_tail(ibox),' [K] '
	  			Write(*,'(A,4X,E15.7,A)') 'Total Coulomb                  ', ewld_tot(ibox),' [K] '
	  			Write(*,'(A,4X,E15.7,A)') '   Real Space (Inter-molecular)', ewld_real(ibox),' [K] '
	  			Write(*,'(A,4X,E15.7,A)') '   Fourier Space               ', ewld_fourier(ibox),' [K] '
	  			Write(*,'(A,4X,E15.7,A)') '   Self Correction             ', ewldself,' [K] '
	  			Write(*,'(A,4X,E15.7,A)') '   Intra-molecular             ', ewldintra,' [K] '
	  			Write(*,'(A,4X,E15.7,A)') '   Slab Correction             ', EETOK*ewld_slab(ibox),' [K] '
	  			Write(*,'(A,4X,E15.7,A)') 'External Field                 ', energy_field(ibox),' [K] '


	  		End Do


	  		

	  		! Block averages for equilibrium stage
	  		Write(*,*) 
	  		Write(*,'(A,I0,A)') 'Block Averages for Equilibrium Stage (', n_blocks_equil, ' blocks)'
	  		Write(*,*) '                   Units     Type    Box     Average           Standard Deviation'

	  		! Loop over the boxes
	  		Do ibox = 1, n_box

	  			! Initalize local variables
	  			nmolavg = 0.0d0
	  			nmolstd = 0.0d0
	  			rhoavg = 0.0d0
	  			rhostd = 0.0d0
	  			engavg = 0.0d0
	  			engstd = 0.0d0
	  			engdispstd = 0.0d0
	  			engdispavg = 0.0d0
	  			n_total = 0.0d0

	  			! Loop over blocks
	  			Do iblk = 1, n_blocks_equil

	  				! Update sum
	  				engavg = engavg + avg_eng(iblk,ibox)
	  				engdispavg = engdispavg + avg_eng_disp(iblk,ibox)

	  				! Loop over molecule types
	  				Do itype = 1, n_mol_types

	  					! Update sum
	  					nmolavg(itype) = nmolavg(itype) + avg_nmol(iblk,itype,ibox)
	  					rhoavg(itype) = rhoavg(itype) + avg_rho(iblk,itype,ibox)
	  					
	  				End Do
	  				
	  			End Do

	  			! Compute the averages
	  			engavg = engavg/DBLE(n_blocks_equil)
	  			engdispavg = engdispavg/DBLE(n_blocks_equil)

	  			Do itype = 1, n_mol_types
	  				nmolavg(itype) = nmolavg(itype)/DBLE(n_blocks_equil)
	  				rhoavg(itype) = rhoavg(itype)/DBLE(n_blocks_equil)
	  				n_total = n_total + nmolavg(itype)
	  			End Do

	  			! Start Computing standard deviation
	  			! Loop over blocks
	  			Do iblk = 1, n_blocks_equil
	  				
	  				! Update sum
	  				engstd = engstd + (avg_eng(iblk,ibox) - engavg)**2
	  				engdispstd = engdispstd + (avg_eng_disp(iblk,ibox) - engdispavg)**2

	  				! Loop over molecule types
	  				Do itype = 1, n_mol_types

	  					! Update sum
	  					nmolstd(itype) = nmolstd(itype) + (avg_nmol(iblk,itype,ibox) - nmolavg(itype))**2
	  					rhostd(itype) = rhostd(itype) + (avg_rho(iblk,itype,ibox) - rhoavg(itype))**2
	  				End Do
	  			End Do

	  			! Compute standard deviations
	  			engstd = dSQRT(engstd/DBLE(n_blocks_equil))
	  			engdispstd = dSQRT(engdispstd/DBLE(n_blocks_equil))

	  			Do itype = 1, n_mol_types
	  				nmolstd(itype) = dSQRT(nmolstd(itype)/DBLE(n_blocks_equil))
	  				rhostd(itype) = dSQRT(rhostd(itype)/DBLE(n_blocks_equil))
	  			End Do
	  			
	  			! Dump stats to the screen
	  			Write(*,'(A,8X,A,14X,I0,5X,E15.7,3X,E15.7)') 'Total energy','J/mol', ibox, R*engavg,R*engstd
	  			Write(*,'(A,21X,A,16X,I0,5X,E15.7,3X,E15.7)') ' ','K', ibox, n_total*engavg, n_total*engstd
!	  			Write(*,'(A,8X,A,14X,I0,5X,E15.7,3X,E15.7)') 'Inter vdW   ','J/mol', ibox, R*engdispavg,R*engdispstd
!	  			Write(*,'(A,21X,A,16X,I0,5X,E15.7,3X,E15.7)') ' ','K', ibox, n_total*engdispavg, n_total*engdispstd

	  			Do itype = 1, n_mol_types
	  				Write(*,'(A,16X,I2,6X,I0,5X,F10.3,8X,E15.7)') &
	  				&'Molecule Number',itype,ibox,nmolavg(itype),nmolstd(itype)
	  				Write(*,'(A,6X,A,7X,I2,6X,I0,5X,E15.7,3X,E15.7)') &
	  				&'Number Density', 'A^-3', itype,ibox,rhoavg(itype), rhostd(itype)
	  				Write(*,'(A,6X,A,7X,I2,6X,I0,5X,E15.7,3X,E15.7)') &
	  				&'Density       ', 'g/ml', itype,ibox,(mol_mass(itype)/(Na*1.0d-24))*rhoavg(itype),&
	  				& (mol_mass(itype)/(Na*1.0d-24))*rhostd(itype)
	  			End Do

	  			
	  		End Do



	  		! Block averages for production stage
	  		Write(*,*) 
	  		Write(*,'(A,I0,A)') 'Block Averages for Production Stage (', n_blocks_prod, ' blocks)'
	  		Write(*,*) '                   Units     Type    Box     Average           Standard Deviation'

	  		! Loop over the boxes
	  		Do ibox = 1, n_box

	  			! Initalize local variables
	  			nmolavg = 0.0d0
	  			nmolstd = 0.0d0
	  			rhoavg = 0.0d0
	  			rhostd = 0.0d0
	  			engavg = 0.0d0
	  			engstd = 0.0d0
	  			engdispavg = 0.0d0
	  			engdispstd = 0.0d0
	  			n_total = 0.0d0

	  			! Initialize widom insertion variables (chemical potential)
	  			muavg = 0.0d0
	  			mustd = 0.0d0


	  			! Loop over blocks
	  			Do iblk = n_blocks_equil+1, n_blocks_tot

	  				! Update sum
	  				engavg = engavg + avg_eng(iblk,ibox)
	  				engdispavg = engdispavg + avg_eng_disp(iblk,ibox)

	  				! Loop over molecule types
	  				Do itype = 1, n_mol_types

	  					! Update sum
	  					nmolavg(itype) = nmolavg(itype) + avg_nmol(iblk,itype,ibox)
	  					rhoavg(itype) = rhoavg(itype) + avg_rho(iblk,itype,ibox)
	  					muavg(itype) = muavg(itype) + widom_sample(iblk,itype,ibox)
	  					
	  				End Do
	  				
	  			End Do

	  			! Compute the averages
	  			engavg = engavg/DBLE(n_blocks_prod)
	  			engdispavg = engdispavg/DBLE(n_blocks_prod)

	  			Do itype = 1, n_mol_types
	  				nmolavg(itype) = nmolavg(itype)/DBLE(n_blocks_prod)
	  				rhoavg(itype) = rhoavg(itype)/DBLE(n_blocks_prod)
	  				n_total = n_total + nmolavg(itype)
	  				muavg(itype) = muavg(itype)/DBLE(n_blocks_prod)
	  			End Do

	  			! Start Computing standard deviation
	  			! Loop over blocks
	  			Do iblk = n_blocks_equil+1, n_blocks_tot
	  				
	  				! Update sum
	  				engstd = engstd + (avg_eng(iblk,ibox) - engavg)**2
	  				engdispstd = engdispstd + (avg_eng_disp(iblk,ibox) - engdispavg)**2

	  				! Loop over molecule types
	  				Do itype = 1, n_mol_types

	  					! Update sum
	  					nmolstd(itype) = nmolstd(itype) + (avg_nmol(iblk,itype,ibox) - nmolavg(itype))**2
	  					rhostd(itype) = rhostd(itype) + (avg_rho(iblk,itype,ibox) - rhoavg(itype))**2	  					
	  					mustd(itype) = mustd(itype) + (widom_sample(iblk,itype,ibox) - muavg(itype))**2
	  					
	  				End Do
	  			End Do

	  			! Compute standard deviations
	  			engstd = dSQRT(engstd/DBLE(n_blocks_prod))
	  			engdispstd = dSQRT(engdispstd/DBLE(n_blocks_prod))

	  			Do itype = 1, n_mol_types
	  				nmolstd(itype) = dSQRT(nmolstd(itype)/DBLE(n_blocks_prod))
	  				rhostd(itype) = dSQRT(rhostd(itype)/DBLE(n_blocks_prod))
	  				mustd(itype) = dSQRT(mustd(itype)/DBLE(n_blocks_prod))
	  				
	  			End Do
	  			
	  			! Dump stats to the screen
	  			Write(*,'(A,8X,A,14X,I0,5X,E15.7,3X,E15.7)') 'Total energy','J/mol', ibox, R*engavg,R*engstd
	  			Write(*,'(A,21X,A,16X,I0,5X,E15.7,3X,E15.7)') ' ','K', ibox, n_total*engavg, n_total*engstd
!	  			Write(*,'(A,8X,A,14X,I0,5X,E15.7,3X,E15.7)') 'Inter vdW   ','J/mol', ibox, R*engdispavg,R*engdispstd
!	  			Write(*,'(A,21X,A,16X,I0,5X,E15.7,3X,E15.7)') ' ','K', ibox, n_total*engdispavg, n_total*engdispstd

	  			Do itype = 1, n_mol_types
	  				Write(*,'(A,16X,I2,6X,I0,5X,F10.3,8X,E15.7)') &
	  				&'Molecule Number',itype,ibox,nmolavg(itype),nmolstd(itype)

	  				Write(*,'(A,6X,A,7X,I2,6X,I0,5X,E15.7,3X,E15.7)') &
	  				&'Number Density', 'A^-3', itype,ibox,rhoavg(itype), rhostd(itype)

	  				Write(*,'(A,6X,A,7X,I2,6X,I0,5X,E15.7,3X,E15.7)') &
	  				&'Density       ', 'g/ml', itype,ibox,(mol_mass(itype)/(Na*1.0d-24))*rhoavg(itype),&
	  				& (mol_mass(itype)/(Na*1.0d-24))*rhostd(itype)

	  				Write(*,'(A,6X,A,7X,I2,6X,I0,5X,F15.7,3X,E15.7)') &
	  				&'Mu            ', ' K  ', itype,ibox,muavg(itype),mustd(itype)

	  				
	  			End Do
	
	  		End Do



	  		! Dump all blocks stats to the screen
	  		Write(*,*)
	  		Write(*,*) '---------- Block Averages ----------'

	  		! Loop over the boxes
	  		Do ibox = 1, n_box

	  			Write(*,'(A,I0)') 'Box: ',ibox
	  			Write(*,*) 'Block    Energy[J/mol]           Mol Number'

	  			! Loop over all blocks
	  			Do iblk = 1, n_blocks_tot

	  				Write(*,'(I5,4X,E15.7,6X,4F10.3)') &
	  				& iblk, R*avg_eng(iblk,ibox), (avg_nmol(iblk,itype,ibox), itype = 1, n_mol_types)
	  			End Do

	  		End Do




	  	! Throw an error
	  	Case default

	  		Write(*,*) 'FATAL ERROR: INVALID CASE FLAG IN STATS.F90'
	  		STOP
	  
	  End Select
	  
	  	
	  
	  	Return
	  
	  End Subroutine 









