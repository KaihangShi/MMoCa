! ==========================================================
! This subroutine is used to calculate the reciprocal energy
! change, self/intra contribution change and slab correction 
! energy change.
! 
! Final output energy in unit of [Kelvin]
!
! Created on 1-11-2017 by Kaihang Shi
! ==========================================================

	  Subroutine ewld_change(selector,ibox,imol,rx_s_old,ry_s_old,rz_s_old,ewldchange)

	  Use global

	  IMPLICIT NONE

	  ! Passed 
	  Integer :: selector, ibox, imol
	  Double Precision, Dimension(n_sites_max) :: rx_s_old, ry_s_old, rz_s_old
	  Double Precision :: ewldchange

	  ! Local
	  Integer :: itype, isite, isitetype, jtype, jsite, jsitetype
	  Integer :: k, n_k, kx, ky, kz, ksq
	  DOUBLE PRECISION,DIMENSION(3) :: two_Pibox
	  COMPLEX*16  :: eikr_new, eikr_old 
	  Double Precision :: ewldintrachange, ewldselfchange, ewldsurfchange



	  ! Initialize variables
	  ewldchange = 0.0d0
	  ewldintrachange = 0.0d0
	  ewldselfchange = 0.0d0
	  ewldsurfchange = 0.0d0
	  two_Pibox(:) = two_Pi/box(:,ibox)
	  eikx_mol = 0.0d0
	  eiky_mol = 0.0d0
	  eikz_mol = 0.0d0
	  ewld_surfc_new(ibox) = ewld_surfc(ibox)

	  ! Get imol's type
	  itype = mol_type(imol,ibox)

	  ! If we are not deleting imol, set new reciprocal vectors (using NEW coordinates)
	  If (selector .NE. 3) Then
	  
	    ! Loop over the sites on the molecule
	    DO isite = 1,n_sites(itype)

	      ! Get the site type 
	      isitetype = site_type(isite,itype)

	      If (q(isitetype,itype,isitetype,itype) .ne. 0.0d0) Then

	      	! Calculate the reciprocal vector components for k = 0
		    eikx_mol(isite,0) = (1.0d0, 0.0d0)
		    eiky_mol(isite,0) = (1.0d0, 0.0d0) 
		    eikz_mol(isite,0) = (1.0d0, 0.0d0)

		    ! Calculate the reciprocal vector components for k = 1
		    eikx_mol(isite,1) = &
		    	& dCMPLX(dCOS(two_Pibox(1)*rx_s(isite,imol,ibox)),dSIN(two_Pibox(1)*rx_s(isite,imol,ibox)))
		    eiky_mol(isite,1) = &
		    	& dCMPLX(dCOS(two_Pibox(2)*ry_s(isite,imol,ibox)),dSIN(two_Pibox(2)*ry_s(isite,imol,ibox)))
		    eikz_mol(isite,1) = &
		    	& dCMPLX(dCOS(two_Pibox(3)*rz_s(isite,imol,ibox)),dSIN(two_Pibox(3)*rz_s(isite,imol,ibox)))

		    ! Calculate the reciprocal vector components for k = -1 (NOTE : kx > 0)
		    eiky_mol(isite,-1) = dCONJG(eiky_mol(isite,1))
		    eikz_mol(isite,-1) = dCONJG(eikz_mol(isite,1))

	      End If	    	      
	      
	    ENDDO

	    ! Calculate the remaining reciprocal vector components by recursion
		
		DO isite = 1,n_sites(itype)
			isitetype = site_type(isite,itype)

			If (q(isitetype,itype,isitetype,itype) .ne. 0.0d0) Then

				! x-component 
				Do k = 2, k_max(1,ibox)
					eikx_mol(isite,k) = eikx_mol(isite,k-1)*eikx_mol(isite,1) 
				End Do
				
	         
				! y-component
				Do k = 2, k_max(2,ibox)
					eiky_mol(isite,k) = eiky_mol(isite,k-1)*eiky_mol(isite,1)
					eiky_mol(isite,-k) = dCONJG(eiky_mol(isite,k))
				End Do
				

	        	! z-component
	        	Do k = 2, k_max(3,ibox)
	        		eikz_mol(isite,k) = eikz_mol(isite,k-1)*eikz_mol(isite,1)
	        		eikz_mol(isite,-k) = dCONJG(eikz_mol(isite,k))
	        	End Do
	        	
	      	End If      

	  	ENDDO
		

	  End If

	  ! Translation or rotation
	  If (selector .eq. 1) Then
		
		! Calculate surface correction term (for slab geometry)
	  	If (lslabc) Then
	  		
		  	! Loop over the sites on the molecule
		    DO isite = 1,n_sites(itype)

		      ! Get the site type 
		      isitetype = site_type(isite,itype)

		      If (q(isitetype,itype,isitetype,itype) .ne. 0.0d0) Then

		      	ewld_surfc_new(ibox) = ewld_surfc_new(ibox) + &
		      		& q(isitetype,itype,isitetype,itype)*(rz_s(isite,imol,ibox) - rz_s_old(isite))
		      	
		      End If	    	

		    ENDDO

		    ! Calculate surface term change
		    ewldsurfchange = two_Pi/vol(ibox)*(ewld_surfc_new(ibox)**2 - ewld_surfc(ibox)**2)
	  	End If

	  	
	  ! Insertion
	  Else if (selector .eq. 2) Then

	  	! Calculate surface correction term for slab geometry
	  	If (lslabc) Then
	  		
		  	! Loop over the sites on the molecule
		    DO isite = 1,n_sites(itype)

		      ! Get the site type 
		      isitetype = site_type(isite,itype)

		      If (q(isitetype,itype,isitetype,itype) .ne. 0.0d0) Then

		      	ewld_surfc_new(ibox) = ewld_surfc_new(ibox) + &
		      		& q(isitetype,itype,isitetype,itype)*(rz_s(isite,imol,ibox) - 0.0d0)
		      	
		      End If	    	

		    ENDDO

		    ! Calculate surface term change
		    ewldsurfchange = two_Pi/vol(ibox)*(ewld_surfc_new(ibox)**2 - ewld_surfc(ibox)**2)

	  	End If
	  	
	  	! Calculate the change in self contribution
	  	ewldselfchange = ewld_self(itype)

	  	! Calculate the change in intra contribution
	  	ewldintrachange = ewld_intra(itype)

	  	! Make sure the old vectors are zero for newly inserted molecule
	  	eikx(ibox,imol,:,:) = 0.0d0
	  	eiky(ibox,imol,:,:) = 0.0d0
	  	eikz(ibox,imol,:,:) = 0.0d0

	  ! Deletion
	  Else if (selector .eq. 3) Then

	  	! Calculate surface correction term for slab geometry
	  	If (lslabc) Then
	  		
		  	! Loop over the sites on the molecule
		    DO isite = 1,n_sites(itype)

		      ! Get the site type 
		      isitetype = site_type(isite,itype)

		      If (q(isitetype,itype,isitetype,itype) .ne. 0.0d0) Then

		      	ewld_surfc_new(ibox) = ewld_surfc_new(ibox) + &
		      		& q(isitetype,itype,isitetype,itype)*(0.0 - rz_s_old(isite))
		      	
		      End If	    	

		    ENDDO

		    ! Calculate surface term change
		    ewldsurfchange = two_Pi/vol(ibox)*(ewld_surfc_new(ibox)**2 - ewld_surfc(ibox)**2)

	  	End If
	  	
	  	! Calculate the change in self contribution
	  	ewldselfchange = 0.0d0 - ewld_self(itype)

	  	! Calculate the change in intra contribution
	  	ewldintrachange = 0.0d0 - ewld_intra(itype)

	  	! 'New' eikr_new are already zero for the deleted molecules

	  End If


	  ! Calculate the new structure factor and the change in reciprocal space energy

	  ! Initialize the vector index and structure factor
	  n_k = 0
	  skewld_new = skewld
	  skewld_ex_new = skewld_ex

	  DO kx = 0,k_max(1,ibox)
	    DO  ky = -k_max(2,ibox),k_max(2,ibox)
	      DO kz = -k_max(3,ibox),k_max(3,ibox)

	        ! Calculate the square of the vector
	        ksq = kx**2 + ky**2 + kz**2

	        ! Check the bounds of the vector
	        IF((ksq .LT. ksq_max(ibox)) .AND. (ksq .NE. 0)) THEN

	        	! Update the number of vectors
	        	n_k = n_k + 1


    			! Check mixing rules
    			! Explicit implies there is a external structure
    			If (mix_rule .eq. EXPLICIT) Then

    				! Calculate structure factors change between adsorbate sites and surface sites
    				! Explicitly exclude the self/intra correction

          			! Loop over imol's sites
          			DO isite = 1,n_sites(itype)

        				! Get the site type
          				isitetype = site_type(isite,itype)

          				! Calculate the structure factor (perturbed for adsorbates)
      					If (qex(AS,isitetype,itype,isitetype,itype) .ne. 0.0d0) Then

      						skewld_ex_new(2,ibox,n_k) = skewld_ex_new(2,ibox,n_k) + &
      						& qex(AS,isitetype,itype,isitetype,itype)*&
      						& (eikx_mol(isite,kx)*eiky_mol(isite,ky)*eikz_mol(isite,kz) - &
      						& eikx(ibox,imol,isite,kx)*eiky(ibox,imol,isite,ky)*eikz(ibox,imol,isite,kz))

      					End If

      				! End loop over imol's sites	
      				End do
  

            		! Calculate the change in reciprocal space energy between adsorbates-surface
            		ewldchange = ewldchange + k_vec(ibox,n_k) * &
            			& ((skewld_ex_new(1,ibox,n_k)* dCONJG(skewld_ex_new(2,ibox,n_k)) + &
              			&  dCONJG(skewld_ex_new(1,ibox,n_k))* skewld_ex_new(2,ibox,n_k)) - &
              			&   (skewld_ex(1,ibox,n_k)* dCONJG(skewld_ex(2,ibox,n_k)) + &
              			&  dCONJG(skewld_ex(1,ibox,n_k))* skewld_ex(2,ibox,n_k)))

    			End If

            	! Loop over imol's charges
            	DO isite = 1,n_sites(itype)

            		! Get the site type
              		isitetype = site_type(isite,itype)

              		! Calculate the structure factor (unperturbed for adsorbates)
              		If (q(isitetype,itype,isitetype,itype) .ne. 0.0d0) Then

              			! Calculate the old exponential factor
              			eikr_old = eikx(ibox,imol,isite,kx)*eiky(ibox,imol,isite,ky)*eikz(ibox,imol,isite,kz)

              			! Calculate the new exponential factor
              			eikr_new = eikx_mol(isite,kx)*eiky_mol(isite,ky)*eikz_mol(isite,kz)

              			! Calculate the structure factor
              			skewld_new(ibox,n_k) = skewld_new(ibox,n_k) + &
              				& q(isitetype,itype,isitetype,itype)*(eikr_new - eikr_old)


              		End If
              
              	ENDDO
	          
	          	! Calculate the reciprocal space energy
          		ewldchange = ewldchange + k_vec(ibox,n_k) * &
          		& ((dCONJG(skewld_new(ibox,n_k))*skewld_new(ibox,n_k)) - (dCONJG(skewld(ibox,n_k))*skewld(ibox,n_k)))

	        ENDIF

	      ENDDO
	    ENDDO
	  ENDDO


	  ! Calculate the final change in Ewald part and convert units from [EE] to [K]
	  ewldchange = (ewldchange - ewldselfchange - ewldintrachange + ewldsurfchange) * EETOK

	  
	  	Return
	  
	  End Subroutine 








