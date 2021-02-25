! ==========================================================
! This subroutine is used to set up Ewald sum parameters
! Created on 1-9-2017 by Kaihang Shi
! Borrow from Jeremy Palmer's set_ewld subroutine
! 
! Note on 1-28-2017: We modified ewald scheme using same method 
! as that used in Towhee for 'ewald_fix_kmax'
! ==========================================================

	  Subroutine set_ewld(selector,ibox)

	  Use global

	  IMPLICIT NONE

	  ! Passed 
	  Integer :: selector, ibox

	  ! Local
	  INTEGER :: n_k,ksq,kx,ky,kz,k
	  INTEGER :: itype,ichrg,isite,ibx
	  DOUBLE PRECISION :: rksq,rkx,rky,rkz
	  DOUBLE PRECISION :: b
	  DOUBLE PRECISION,DIMENSION(3) :: two_Pibox
	  INTEGER :: ierr

	  DOUBLE PRECISION :: testeng


	  ! Select a case
	  SELECT CASE (selector)

	    ! Perform one time allocation tasks
	    CASE(0)


	      If (ewld_style .eq. ewald_fix_kmax) Then
	      	! Loop over the boxes
		      	DO ibx = 1,n_box

		        	! Calculate maxk
		        	maxk(ibx) = 3*MAXVAL(k_max(:,ibx))**3

		        	! Calculate ksq_max
		        	ksq_max(ibx)= MAXVAL(k_max(:,ibx))**2 + 2

		      	ENDDO
	      End If

	      ! Initialize the error hanlder
	      ierr = 0

	      ! Allocate the arrays 
	      ALLOCATE(k_vec(n_box,MAXVAL(maxk)), STAT=ierr)
	      ALLOCATE(eikx(n_box,n_mol_max,n_sites_max,0:MAXVAL(k_max)), STAT=ierr)
	      ALLOCATE(eiky(n_box,n_mol_max,n_sites_max,-MAXVAL(k_max):MAXVAL(k_max)), STAT=ierr)
	      ALLOCATE(eikz(n_box,n_mol_max,n_sites_max,-MAXVAL(k_max):MAXVAL(k_max)), STAT=ierr)
	      ALLOCATE(eikx_mol(n_sites_max,0:MAXVAL(k_max)), STAT=ierr)
	      ALLOCATE(eiky_mol(n_sites_max,-MAXVAL(k_max):MAXVAL(k_max)), STAT=ierr)
	      ALLOCATE(eikz_mol(n_sites_max,-MAXVAL(k_max):MAXVAL(k_max)), STAT=ierr)
	      ALLOCATE(skewld(n_box,MAXVAL(maxk)), STAT=ierr)
	      ALLOCATE(skewld_new(n_box,MAXVAL(maxk)), STAT=ierr)
	      ALLOCATE(skewld_ex(2,n_box,MAXVAL(maxk)), STAT=ierr)
	      ALLOCATE(skewld_ex_new(2,n_box,MAXVAL(maxk)), STAT=ierr)
	      Allocate(ewld_surfc(n_box), STAT=ierr)
	      Allocate(ewld_surfc_new(n_box), STAT=ierr)

	      

	      ! Check for allocation error
	      IF(ierr .NE. 0) THEN
	        WRITE(*,*) 'FATAL ERROR: ALLOCATION OF ARRAYS FAILED IN SUBROUTINE SET_EWLD '
	        STOP
	      ENDIF


	    ! Set prefactor in the sum of k vectors in reciprocal space
	    CASE(1)

	      
	      ! Setup calculate reciprocal vectors
	 
	      ! Calculate the denominator for the k-space exponential
	      b = (1.0d0/(4.0d0*alpha**2))
	 
	      ! Initialize the number of vectors
	      n_k = 0
	     
	      ! Loop over all the vectors
	      ! vectors follow Stanford Ewald notes 
	      DO kx = 0,k_max(1,ibox)
	        rkx = two_Pi*DBLE(kx)/box(1,ibox)
	        DO ky = -k_max(2,ibox),k_max(2,ibox)
	          rky = two_Pi*DBLE(ky)/box(2,ibox)
	          DO kz = -k_max(3,ibox),k_max(3,ibox)
	            rkz = two_Pi*DBLE(kz)/box(3,ibox)

	            ! Calculate the square of the vector
	            ksq = kx**2 + ky**2 + kz**2

	            ! Check the bounds of the vector
	            IF((ksq .LT. ksq_max(ibox)) .AND. (ksq .NE. 0)) THEN

	              ! Update the number of vectors
	              n_k = n_k + 1

	              If (n_k .gt. MAXVAL(maxk)) Then
	              	Write(*,*) 'FATAL ERROR: TOO SMALL OF maxk'
	              	STOP
	              End If

	              ! Calculate rksq
	              rksq = rkx**2 + rky**2 + rkz**2

	              ! Calculate the the k-vector 
	              IF (kx .GT. 0) THEN
	                k_vec(ibox,n_k) = 2.0d0*two_Pi*(1.0d0/vol(ibox))*dEXP(-b*rksq)/rksq
	              ELSE
	                k_vec(ibox,n_k) = two_Pi*(1.0d0/vol(ibox))*dEXP(-b*rksq)/rksq
	              ENDIF      

	            ENDIF
	          ENDDO
	        ENDDO
	      ENDDO

	      ! Print number of vectors to screen
!	      WRITE(*,'(A,I0,A,I0)') 'Box ',ibox, ' k-vectors : ', n_k


		! Set/reset rcelect, alpha, self/intra values
		Case(2)

			If (ewld_style .eq. ewald_fix_kmax) Then
				! ewald_fix_kmax scheme follows the same one in Towhee
				! Calculate corresponding convergence parameter 
		      	alpha = kalp/MINVAL(box(:,ibox))

		   		! Set cut_off raidus in real space to half of the box length
		   		rcelect = MINVAL(box(:,ibox))/2.0d0
		   		rcelectsq = rcelect**2

			  	! Compute the self correction & intramolecular contribution in k-space for each molecule type
		      	CALL set_ewldself

		      	
		    Else 
		    	Write(*,*) 'SET_EWLD: INVALID EWLD_STYLE'
		    	STOP
			End If

			
	     
	      
	    CASE DEFAULT

	      WRITE(*,*) 'FATAL ERROR: INVALID FLAG TYPE IN SET_EWALD SUBROUTINE'
	      STOP


	  END SELECT

	    
	  	Return
	  
	  End Subroutine 