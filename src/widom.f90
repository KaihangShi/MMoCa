! ==========================================================
! This subroutine is used to perform Widom Insertion method 
! and compute Widom statistics for block average
! Created on 1-26-2017 by Kaihang
!
! Notes on 1-31-2017 by Kaihang:
! Modified the code, take all trial moves into account. We 
! must accept all the trial moves instead of reinsert it.
! ==========================================================


      Subroutine widom(selector,iblock,istep)

      Use global

      IMPLICIT NONE

      ! Passed
      Integer :: selector, iblock, istep

      ! Local
      Integer :: ibox, itype, imol, iitype, jjtype
      Double Precision, Dimension(n_sites_max) :: rx_s_old, ry_s_old, rz_s_old
      Double Precision :: eng_new, eng_change, ewldchange, engtailold, engtailnew



      Select Case (selector)
      	
      	! Perform one-time insertion
      	! Input: iblock, istep
      	Case(1)

      		! Loop over the boxes
	        DO ibox = 1,n_box 
	        	! Loop over the molecule types
	        	DO itype = 1,n_mol_types

	          		! Check if a widom trial should be attempted
	          		IF((widom_freq(itype,ibox) .EQ. 0) .OR. MOD(istep,widom_freq(itype,ibox)) .NE. 0) CYCLE

			        ! Set the molecule's array index 
			        imol = n_mol_tot(ibox) + 1

			        ! Set the molecule's type (used in ewld_change subroutine)
			        mol_type(imol,ibox) = itype
			 
			        ! Generate random position
			        rx(imol,ibox) = random(idum)*box(1,ibox)
			        ry(imol,ibox) = random(idum)*box(2,ibox)
			        rz(imol,ibox) = random(idum)*box(3,ibox)

			        ! Generate orientation
			        CALL rand_quat(q1(imol,ibox),q2(imol,ibox),q3(imol,ibox),q4(imol,ibox))
			  
			        ! Calculate site coordinates
			        CALL site_coords(imol,ibox) 

			        ! Calculate the energy of the new configuration
			        CALL eng_mol(ibox,imol,eng_new) 

			        ! Update the number of widom samples
			        ! Must include all trials
			        widom_stat(itype,ibox) = widom_stat(itype,ibox) + 1

			        ! Check for overlap 
			        If (OVERLAP) Then

			        	! if overlap, boltzmann factor is zero
			        	! Set the value to zero to avoid overflow
			        	widom_sample(iblock,itype,ibox) = widom_sample(iblock,itype,ibox) + 1.0d-300

			        	! Cycle to next trial
			        	CYCLE

			        End If

			        ! Calculate the change in energy and virial
			        eng_change = eng_new - 0.0d0


			        ! Check for Ewald summation
			        IF(lewld) THEN

			        	! Initialize dummy arrays (redundant but compiler may complain)
			            rx_s_old = 0.0d0
			            ry_s_old = 0.0d0
			            rz_s_old = 0.0d0

			            ! Calculate the Ewald contribution (including reciprocal, self and intra term change)
			            CALL ewld_change(2,ibox,imol,rx_s_old,ry_s_old,rz_s_old,ewldchange)

			            ! Update the energy change [K]
			            eng_change = eng_change + ewldchange

			        ENDIF

			        ! Check for vdW tail correction
			        ! Because n_mol changes, so energy_tail changes accordingly
			        If (ltailc) Then

			        	! Initialize new tail energy
			        	engtailnew = 0.0d0

			       		! Store old vdW tail energy
			       		engtailold = energy_tail(ibox)

			       		! Temporarily increase itype's mol number
			       		n_mol(itype,ibox) = n_mol(itype,ibox) + 1

				 	  	! Loop over all molecule types
				 	  	Do iitype = 1, n_mol_types
				 	  		Do jjtype = 1, n_mol_types

				 	  			! Apply prefactor N*rho and convert to unit of [K]
				 	  			engtailnew = engtailnew + &
				 	  				& n_mol(iitype,ibox)*n_mol(jjtype,ibox)/vol(ibox)*vdw_tail(iitype,jjtype)
				 	  		End Do
				 	  	End Do

				 	  	! Update the energy change [K]
				 	  	eng_change = eng_change + engtailnew - engtailold

				 	  	! Restore itype's mol number
				 	  	n_mol(itype,ibox) = n_mol(itype,ibox) - 1

				 	End If

			  
			        ! NPT
			        IF(ensmbl .EQ. ENS_NPT) widom_sample(iblock,itype,ibox) = &
			        						& widom_sample(iblock,itype,ibox) +  vol(ibox)*dEXP(-eng_change/temp)

			          
	        	! End loop over the molecule types
	        	ENDDO	   
	      	! End loop over the boxes
	      	ENDDO 

	    ! Compute the block average
	    ! Input: iblock
	    Case(2)

	    	! Loop over boxes
	    	Do ibox = 1, n_box
	    		! Loop over each molecule type
	    		Do itype = 1, n_mol_types

	    			! Compute block average
          			IF(widom_stat(itype,ibox) .GT. 0) THEN
            			widom_sample(iblock,itype,ibox) = widom_sample(iblock,itype,ibox)/DBLE(widom_stat(itype,ibox))
          			ELSE
            			widom_sample(iblock,itype,ibox) = 0.0d0
          			ENDIF       

            		! NPT [K]
            		IF(ensmbl .EQ. ENS_NPT) widom_sample(iblock,itype,ibox) = muid(itype) - &
            			& temp*dLOG((PVTOK*press*widom_sample(iblock,itype,ibox))/(temp*DBLE((n_mol_tot(ibox)+1))))

	    			
	    		End Do
	    		
	    	End Do
      		

      		
      	Case default
      		Write(*,*) 'FATAL ERROR: INVALID SELECTOR IN WIDOM SUBROUTINE'
      		STOP
      		
      
      End Select
      
      	
      
      	Return
      
      End Subroutine 