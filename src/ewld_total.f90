! ==========================================================
! This subroutine is used to calculate Fourier space &
! self/intra term (from Fourier space) & Surface term 
! in Ewald Sum
! Created on 1-10-2017 by Kaihang Shi
!
! Noted by Kaihang on 2-3-2017:
! Surface term is added as slab geometry correction
! ==========================================================

	  Subroutine ewld_total(ibox,engewld)

	  Use global

	  IMPLICIT NONE

	  ! Passed
	  Integer :: ibox
	  ! engewld = k-space - self_correction + surface term
	  Double Precision :: engewld

	  ! Local
	  INTEGER :: n_k,ksq,kx,ky,kz,k
	  INTEGER :: imol, jmol, itype, jtype, isite, isitetype, jsite, jsitetype
	  DOUBLE PRECISION,DIMENSION(3) :: two_Pibox
	  DOUBLE PRECISION :: ewldself, ewldintra

	  
	  ! Calculate the total reciprocal energy  
	  ! Initialize variables
	  engewld = 0.0d0
	  ewldself = 0.0d0
	  ewldintra = 0.0d0
	  two_Pibox(:) = two_Pi/box(:,ibox)
	  ewld_fourier(ibox) = 0.0d0
	  ewld_surfc(ibox) = 0.0d0
	  ewld_slab(ibox) = 0.0d0


	  ! Loop over all of the molecules in ibox
	  DO imol = 1,n_mol_tot(ibox)
	  
	    ! Get imol's type
	    itype = mol_type(imol,ibox)

	    ! Check if external structure
	    if(initstyle(itype,ibox) .eq. 'coords') Then

	    	! Check direct Coulombic interaction for ad-surf feature
	    	If (lcoulsc) Then
	    		! Loop over the sites on the molecule
	    		DO isite = 1,n_sites(itype)

	      			! Get the site type 
	      			isitetype = site_type(isite,itype)

	      			If (q(isitetype,itype,isitetype,itype) .ne. 0.0d0) Then
	      	
		    			! Calculate surface correction term (for slab geometry)
		    			If (lslabc) Then

		    				ewld_surfc(ibox) = ewld_surfc(ibox) + &
		    					& q(isitetype,itype,isitetype,itype)*rz_s(isite,imol,ibox)
		    
		    			End If
	      			End If	    	      
	    		ENDDO

	    		CYCLE

	    	End if

	    End if


	    ! Loop over the sites on the molecule
	    DO isite = 1,n_sites(itype)

	      ! Get the site type 
	      isitetype = site_type(isite,itype)

	      If (q(isitetype,itype,isitetype,itype) .ne. 0.0d0) Then
	      	! Calculate the reciprocal vector components for k = 0
		    eikx(ibox,imol,isite,0) = (1.0d0, 0.0d0)
		    eiky(ibox,imol,isite,0) = (1.0d0, 0.0d0) 
		    eikz(ibox,imol,isite,0) = (1.0d0, 0.0d0)

		    ! Calculate the reciprocal vector components for k = 1
		    eikx(ibox,imol,isite,1) = &
		    	& dCMPLX(dCOS(two_Pibox(1)*rx_s(isite,imol,ibox)),dSIN(two_Pibox(1)*rx_s(isite,imol,ibox)))
		    eiky(ibox,imol,isite,1) = &
		    	& dCMPLX(dCOS(two_Pibox(2)*ry_s(isite,imol,ibox)),dSIN(two_Pibox(2)*ry_s(isite,imol,ibox)))
		    eikz(ibox,imol,isite,1) = &
		    	& dCMPLX(dCOS(two_Pibox(3)*rz_s(isite,imol,ibox)),dSIN(two_Pibox(3)*rz_s(isite,imol,ibox)))

		    ! Calculate the reciprocal vector components for k = -1 (NOTE : kx > 0)
		    eiky(ibox,imol,isite,-1) = dCONJG(eiky(ibox,imol,isite,1))
		    eikz(ibox,imol,isite,-1) = dCONJG(eikz(ibox,imol,isite,1))

		    ! Calculate surface correction term (for slab geometry)
		    If (lslabc) Then

		    	ewld_surfc(ibox) = ewld_surfc(ibox) + q(isitetype,itype,isitetype,itype)*rz_s(isite,imol,ibox)
		    
		    End If

	      End If	    	      
	      
	    ENDDO
	  ENDDO

	  ! Calculate the remaining reciprocal vector components by recursion
      DO imol = 1,n_mol_tot(ibox)
      	itype = mol_type(imol,ibox)

      	! Check if external structure
	    if(initstyle(itype,ibox) .eq. 'coords') Then

	    	! Check direct Coulombic interaction for ad-surf feature
	    	If (lcoulsc) CYCLE

	    End if

      	DO isite = 1,n_sites(itype)
      		isitetype = site_type(isite,itype)

      		If (q(isitetype,itype,isitetype,itype) .ne. 0.0d0) Then

      			! x-component 
      			Do k = 2, k_max(1,ibox)
      				eikx(ibox,imol,isite,k) = eikx(ibox,imol,isite,k-1)*eikx(ibox,imol,isite,1) 
      			End Do
        		
         
       			! y-component
       			Do k = 2, k_max(2,ibox)
       				eiky(ibox,imol,isite,k) = eiky(ibox,imol,isite,k-1)*eiky(ibox,imol,isite,1)
        			eiky(ibox,imol,isite,-k) = dCONJG(eiky(ibox,imol,isite,k))
       			End Do
        		

        		! z-component
        		Do k = 2, k_max(3,ibox)
        			eikz(ibox,imol,isite,k) = eikz(ibox,imol,isite,k-1)*eikz(ibox,imol,isite,1)
        			eikz(ibox,imol,isite,-k) = dCONJG(eikz(ibox,imol,isite,k))
        		End Do
        		
      		End If      
      	ENDDO
      ENDDO
	  


	  ! Calculate the structure factor and reciprocal space energy
	  ! Initialze the vector
	  n_k = 0

	  DO kx = 0,k_max(1,ibox)
	    DO  ky = -k_max(2,ibox),k_max(2,ibox)
	      DO kz = -k_max(3,ibox),k_max(3,ibox)

	        ! Calculate the square of the vector
	        ksq = kx**2 + ky**2 + kz**2

	        ! Check the bounds of the vector
	        IF((ksq .LT. ksq_max(ibox)) .AND. (ksq .NE. 0)) THEN

	          ! Update the number of vectors
	          n_k = n_k + 1

	          ! Initialize the structure factor
	          skewld(ibox,n_k) = (0.0d0, 0.0d0)
	          skewld_ex(:,ibox,n_k) = (0.0d0, 0.0d0)

	          ! Loop over the molecules
	          DO imol = 1,n_mol_tot(ibox)

	            ! Get imol's type
	            itype = mol_type(imol,ibox)

	            ! Check if external structure
	    		if(initstyle(itype,ibox) .eq. 'coords') Then

	    			! Check direct Coulombic interaction for ad-surf feature
	    			If (lcoulsc) CYCLE

	    			! Check mixing rules
	    			If (mix_rule .eq. EXPLICIT) Then

	    				! Calculate structure factors between adsorbate sites and surface sites
	    				! Explicitly exclude the self/intra correction
	    				
	    				! Loop over imol's charges
	            		DO isite = 1,n_sites(itype)

	              			! Get the site type
	              			isitetype = site_type(isite,itype)

	              			! Calculate the structure factor
          					If (qex(AS,isitetype,itype,isitetype,itype) .ne. 0.0d0) Then


          						! Calculate structure factor for surface sites
          						skewld_ex(1,ibox,n_k) = skewld_ex(1,ibox,n_k) + qex(AS,isitetype,itype,isitetype,itype)*&
          						& eikx(ibox,imol,isite,kx)*eiky(ibox,imol,isite,ky)*eikz(ibox,imol,isite,kz)


          					End If

          				! End loop over surface sites	
	            		ENDDO	                   		

	            		! Cycle to next molecule
	            		CYCLE

	    			End If

	    		End if


	            ! Loop over adsorbate's charges
	            DO isite = 1,n_sites(itype)

	              ! Get the site type
	              isitetype = site_type(isite,itype)

	              ! Calculate the structure factor for ad-ad interactions (original charge)
	              If (q(isitetype,itype,isitetype,itype) .ne. 0.0d0) Then

	              	skewld(ibox,n_k) = skewld(ibox,n_k) + q(isitetype,itype,isitetype,itype)*&
	              		& eikx(ibox,imol,isite,kx)*eiky(ibox,imol,isite,ky)*eikz(ibox,imol,isite,kz)


	              End If

	              ! Calculate the separate structure factor for adsorbates (perturbed charge)
	              If ((qex(AS,isitetype,itype,isitetype,itype) .ne. 0.0d0) .and. (mix_rule .eq. EXPLICIT)) Then

	              	skewld_ex(2,ibox,n_k) = skewld_ex(2,ibox,n_k) + qex(AS,isitetype,itype,isitetype,itype)*&
	              		& eikx(ibox,imol,isite,kx)*eiky(ibox,imol,isite,ky)*eikz(ibox,imol,isite,kz)


	              End If
	              
	            ENDDO
	            
	          ! End loop over molecules
	          ENDDO

 
	          ! Calculate the reciprocal space energy (Explicit adsorbates-surface energy, perturbed)
              engewld = engewld + k_vec(ibox,n_k) * (skewld_ex(1,ibox,n_k)* dCONJG(skewld_ex(2,ibox,n_k)) + &
              & dCONJG(skewld_ex(1,ibox,n_k))* skewld_ex(2,ibox,n_k))
	          
	          ! Calculate the reciprocal space energy (add adsorbate-adsorbate part, unperturbed)
	          engewld = engewld + k_vec(ibox,n_k) * dCONJG(skewld(ibox,n_k))*skewld(ibox,n_k)

	        ENDIF

	      ENDDO
	    ENDDO
	  ENDDO


	  ! Update reciprocal space energy
	  ewld_fourier(ibox) = EETOK*engewld  


	  ! Calculate the spurious k-space contribution (self correction + intramolecular)
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

	  ! Calculate surface correction term (for slab geometry) [EE]
	  If (lslabc) Then

	  	ewld_slab(ibox) = two_Pi/vol(ibox)*ewld_surfc(ibox)**2
	  
	  End If

	  ! Calculate the final reciprocal space energy & self correction term 
	  ! and convert to units of [Kelvin]
	  engewld = EETOK*(engewld - ewldself - ewldintra + ewld_slab(ibox)) 
	  
	  
	  Return
	  
	  End Subroutine 






