! ==========================================================
! This subroutine is used to calculate self correction part
! and intramolecular part (from k-space) in Ewald sum
! Both of the terms should be subtracted from the Ewald total 
! to give correct Inter-molecular configurational energy
!
! Created on 1-9-2017 by Kaihang Shi
! Last modified on 2-2-2017 by Kaihang Shi:
! 	Added minimum image convension when calculating intramolecule term.
! ==========================================================

	  Subroutine set_ewldself

	  Use global 

	  IMPLICIT NONE

	  ! Local
	  INTEGER :: itype, isite, isitetype, jsite, jsitetype
	  DOUBLE PRECISION :: rijs, rxijs, ryijs, rzijs


	  ! Initialize the energy array
	  ewld_self = 0.0d0
	  ewld_intra = 0.0d0

	  ! Calculate the self energy contribution 
	  ! Loop over the molecule types
	  DO itype = 1,n_mol_types

	    ! Loop over the charges on the molecule
	    DO isite = 1, n_sites(itype)
	    
	    	! Get isite type
	    	isitetype = site_type(isite,itype)

	    	! Calculate the self interaction
	    	if(qsq(isitetype,itype,isitetype,itype) .ne. 0.0d0) &
	    	& ewld_self(itype) = ewld_self(itype) + qsq(isitetype,itype,isitetype,itype)

	    ENDDO

	    ! Apply the prefactor
	    ewld_self(itype) = alpha*ewld_self(itype)/dSQRT(Pi)
	    
	    ! Calculate the intramolecular contribution 
	    ! Loop over all charge pairs
	    DO isite = 1, n_sites(itype)-1
	      isitetype = site_type(isite,itype)
	      DO jsite = isite+1, n_sites(itype)
	        jsitetype = site_type(jsite,itype)

	        If (qsq(isitetype,itype,jsitetype,itype) .ne. 0.0d0) Then
	        	! Calculate the separation vector between the two sites
		        rxijs = rx_i(jsite,itype) - rx_i(isite,itype)
		        ryijs = ry_i(jsite,itype) - ry_i(isite,itype)
		        rzijs = rz_i(jsite,itype) - rz_i(isite,itype)

		        ! Apply minimum image convension
		        ! It seems like Towhee also apply minimum image convension 
		        rxijs = rxijs - dNINT(rxijs/box(1,1))*box(1,1)
		        ryijs = ryijs - dNINT(ryijs/box(2,1))*box(2,1)
		        rzijs = rzijs - dNINT(rzijs/box(3,1))*box(3,1)

		        ! Calculate scalar separation distance
		        rijs = dSQRT(rxijs**2 + ryijs**2 + rzijs**2)
		      
		        ! Calculate the intramolecular energy of itype
		        ewld_intra(itype) = ewld_intra(itype) + qsq(isitetype,itype,jsitetype,itype)*dERF(alpha*rijs)/rijs

	        End If
	        
	      ENDDO
	    ENDDO

!	    Write(*,'(A,E15.7)') 'Self_term: ', ewld_self(itype)*EETOK
!	    Write(*,'(A,E15.7)') 'Intra_term: ', ewld_intra(itype)*EETOK


	  ! End loop over molecule types
	  ENDDO

	  
	  
	  
	  	
	  
	  	Return
	  
	  End Subroutine 