! ==========================================================
! This subroutine is used to calculate the energy of imol with
! external field 
! Created on 2-2-2017 by Kaihang Shi
! 
! Note: Right now we are assuming the center position of hard wall
!	is at z = 0.0d0. PLEASE CHANGE IT ACCORDING TO SIMULATION PURPOSE
! ==========================================================


	  Subroutine eng_field(ibox,imol,engfield)

	  Use global

	  IMPLICIT NONE

	  ! Passed 
	  Integer :: ibox, imol
	  Double Precision :: engfield

	  ! Local
	  Integer :: itype, isite, isitetype, ibin, xbin, cbin, zbin
	  Double Precision :: upbound, lobound, zz
	  Double Precision :: rxcos, rxx, rzz, zlimit, rsq, dacc

	  ! Initialize
	  engfield = 0.0d0
	  ! Assume no OVERLAP initially
	  OVERLAP = .FALSE.
	  
	  ! Get imol's type
	  itype = mol_type(imol,ibox)

	  ! Return if imol is external structure
	  If (initstyle(itype,ibox) .eq. 'coords') Return

	  ! Choose field type
	  ! Hard wall 
	  If (field_type .eq. HARD_WALL) Then

	  	! Set upper bound and lower bound for hard wall
	  	upbound = box(3,ibox) - wall_radius
	  	lobound = wall_radius
	  	
	  	! Loop over all sites on imol
	  	Do isite = 1, n_sites(itype)
	  		
	  		! Check if isite is interacting with hard wall
	  		If ((rz_s(isite,imol,ibox) .le. lobound) .or. (rz_s(isite,imol,ibox) .ge. upbound)) Then
	  			OVERLAP = .TRUE.
	  			RETURN
	  		End If

	  	End Do

	  	! If no overlap, set energy to zero
	  	engfield = 0.0d0

	  ! Hard slit wall compatible with Yun Long's model
	  Else If (field_type .eq. HARD_SLIT_FINITEX) Then

	  	! Set upper bound and lower bound for hard wall
	  	upbound = box(3,ibox) - wall_radius
	  	lobound = wall_radius
	  	
	  	! Loop over all sites on imol
	  	Do isite = 1, n_sites(itype)
	  		
	  		! Check if isite is interacting with hard wall
	  		If ((rz_s(isite,imol,ibox) .le. lobound) .or. (rz_s(isite,imol,ibox) .ge. upbound)) Then
	  			OVERLAP = .TRUE.
	  			RETURN
	  		End If

	  	End Do

	  	! If no overlap, set energy to zero
	  	engfield = 0.0d0

	  ! 10-4-3 Steele potential for only one side 
	  Else If (field_type .eq. STEELE) Then
	  	
	  	! Set upper bound for hard wall (This could be changed to upper Steele potential later)
	  	upbound = box(3,ibox) - steele_position(1)

	  	! Loop over all sites on imol
	  	Do isite = 1, n_sites(itype)

	  		! Check if isite is inside the wall
	  		If ((rz_s(isite,imol,ibox) .le. steele_position(1)) .or. (rz_s(isite,imol,ibox) .ge. upbound)) Then
	  			OVERLAP = .TRUE.
	  			RETURN
	  		End If

	  		! Get site type
	  		isitetype = site_type(isite,itype)

	  		! Get z distance
	  		zz = rz_s(isite,imol,ibox) - steele_position(1)
	  		! Double check
	  		If(zz .lt. 0.0d0) Then
	  			Write(*,*) 'ENG_FIELD: MOLECULE BELOW GRAPHENE SURFACE (Z < STEELE_POSITION)'
	  			STOP
	  		End if

	  		! Check with cutoff radius
	  		If (zz .gt. steele_cut) CYCLE

	  		! Calculate 10-4-3 Steele potential energy [K]
	  		engfield = engfield + & 
	  			& two_Pi*steele_epsilonsf(isitetype,itype)*steele_rhos*steele_sigmasfsq(isitetype,itype)*steele_delta*&
	  			& (0.4d0 * (steele_sigmasfsq(isitetype,itype)/zz**2)**5 - &
	  			&  (steele_sigmasfsq(isitetype,itype)/zz**2)**2 - &
	  			&  steele_sigmasfsq(isitetype,itype)**2/(3.0d0*steele_delta*(zz+0.61d0*steele_delta)**3))
	  			
	  		
	  	End Do


	  ! 10-4-3 Steele potential for slit pore geometry
	  Else If (field_type .eq. STEELE_SLIT_PORE) Then

	  	! Interactions with BOTTOM walls
	  	! Loop over all sites on imol
	  	Do isite = 1, n_sites(itype)

	  		! Check if isite is inside the wall
	  		If ((rz_s(isite,imol,ibox) .le. steele_position(1)) .or. (rz_s(isite,imol,ibox) .ge. steele_position(2))) Then
	  			OVERLAP = .TRUE.
	  			RETURN
	  		End If

	  		! Get site type
	  		isitetype = site_type(isite,itype)

	  		! Get z distance
	  		zz = rz_s(isite,imol,ibox) - steele_position(1)
	  		! Double check
	  		If(zz .lt. 0.0d0) Then
	  			Write(*,*) 'ENG_FIELD: MOLECULE BELOW GRAPHENE SURFACE (Z < STEELE_POSITION)'
	  			STOP
	  		End if

	  		! Check with cutoff radius
	  		If (zz .gt. steele_cut) CYCLE

	  		! Calculate 10-4-3 Steele potential energy [K]
	  		engfield = engfield + & 
	  			& two_Pi*steele_epsilonsf(isitetype,itype)*steele_rhos*steele_sigmasfsq(isitetype,itype)*steele_delta*&
	  			& (0.4d0 * (steele_sigmasfsq(isitetype,itype)/zz**2)**5 - &
	  			&  (steele_sigmasfsq(isitetype,itype)/zz**2)**2 - &
	  			&  steele_sigmasfsq(isitetype,itype)**2/(3.0d0*steele_delta*(zz+0.61d0*steele_delta)**3))
	  			
	  	! End loop over all sites on imol	
	  	End Do

	  	! Interactions with TOP walls
	  	! Loop over all sites on imol
	  	Do isite = 1, n_sites(itype)

	  		! Get site type
	  		isitetype = site_type(isite,itype)

	  		! Get z distance
	  		zz =  steele_position(2) - rz_s(isite,imol,ibox)
	  		! Double check
	  		If(zz .lt. 0.0d0) Then
	  			Write(*,*) 'ENG_FIELD: MOLECULE BELOW GRAPHENE SURFACE (Z < STEELE_POSITION)'
	  			STOP
	  		End if

	  		! Check with cutoff radius
	  		If (zz .gt. steele_cut) CYCLE

	  		! Calculate 10-4-3 Steele potential energy [K]
	  		engfield = engfield + & 
	  			& two_Pi*steele_epsilonsf(isitetype,itype)*steele_rhos*steele_sigmasfsq(isitetype,itype)*steele_delta*&
	  			& (0.4d0 * (steele_sigmasfsq(isitetype,itype)/zz**2)**5 - &
	  			&  (steele_sigmasfsq(isitetype,itype)/zz**2)**2 - &
	  			&  steele_sigmasfsq(isitetype,itype)**2/(3.0d0*steele_delta*(zz+0.61d0*steele_delta)**3))
	  			
	  	! End loop over all sites on imol	
	  	End Do


	  ! 10-4-3 Steele potential for finite pore model with finite surface in x-direction
	  ! Added on 3-24-2018
	  Else If (field_type .eq. STEELE_SLIT_FINITEX) Then

	  	! Interactions with BOTTOM walls
	  	! Loop over all sites on imol
	  	Do isite = 1, n_sites(itype)

	  		! Check if isite is inside the wall
	  		If ((rz_s(isite,imol,ibox) .le. steele_position(1)) .or. (rz_s(isite,imol,ibox) .ge. steele_position(2))) Then
	  			OVERLAP = .TRUE.
	  			RETURN
	  		End If

	  		! Check if outside the pore in x-direction
	  		If ((rx_s(isite,imol,ibox) .lt. steele_posx(1)) .or. (rx_s(isite,imol,ibox) .gt. steele_posx(2))) CYCLE

	  		! Get site type
	  		isitetype = site_type(isite,itype)

	  		! Get z distance
	  		zz = rz_s(isite,imol,ibox) - steele_position(1)
	  		! Double check
	  		If(zz .lt. 0.0d0) Then
	  			Write(*,*) 'ENG_FIELD: MOLECULE BELOW GRAPHENE SURFACE (Z < STEELE_POSITION)'
	  			STOP
	  		End if

	  		! Check with cutoff radius
	  		If (zz .gt. steele_cut) CYCLE

	  		! Calculate 10-4-3 Steele potential energy [K]
	  		engfield = engfield + & 
	  			& two_Pi*steele_epsilonsf(isitetype,itype)*steele_rhos*steele_sigmasfsq(isitetype,itype)*steele_delta*&
	  			& (0.4d0 * (steele_sigmasfsq(isitetype,itype)/zz**2)**5 - &
	  			&  (steele_sigmasfsq(isitetype,itype)/zz**2)**2 - &
	  			&  steele_sigmasfsq(isitetype,itype)**2/(3.0d0*steele_delta*(zz+0.61d0*steele_delta)**3))
	  			
	  	! End loop over all sites on imol	
	  	End Do

	  	! Interactions with TOP walls
	  	! Loop over all sites on imol
	  	Do isite = 1, n_sites(itype)

	  		! Check if outside the pore in x-direction
	  		If ((rx_s(isite,imol,ibox) .lt. steele_posx(1)) .or. (rx_s(isite,imol,ibox) .gt. steele_posx(2))) CYCLE

	  		! Get site type
	  		isitetype = site_type(isite,itype)

	  		! Get z distance
	  		zz =  steele_position(2) - rz_s(isite,imol,ibox)
	  		! Double check
	  		If(zz .lt. 0.0d0) Then
	  			Write(*,*) 'ENG_FIELD: MOLECULE BELOW GRAPHENE SURFACE (Z < STEELE_POSITION)'
	  			STOP
	  		End if

	  		! Check with cutoff radius
	  		If (zz .gt. steele_cut) CYCLE

	  		! Calculate 10-4-3 Steele potential energy [K]
	  		engfield = engfield + & 
	  			& two_Pi*steele_epsilonsf(isitetype,itype)*steele_rhos*steele_sigmasfsq(isitetype,itype)*steele_delta*&
	  			& (0.4d0 * (steele_sigmasfsq(isitetype,itype)/zz**2)**5 - &
	  			&  (steele_sigmasfsq(isitetype,itype)/zz**2)**2 - &
	  			&  steele_sigmasfsq(isitetype,itype)**2/(3.0d0*steele_delta*(zz+0.61d0*steele_delta)**3))
	  			
	  	! End loop over all sites on imol	
	  	End Do


	  ! Coarse-grained wall potential (from pre-calculated potential file)
	  Else If (field_type .eq. CG_WALL) Then

	  	! Interactions with BOTTOM walls
	  	! Loop over all sites on imol
	  	Do isite = 1, n_sites(itype)

	  		! Check if isite is inside the wall
	  		If ((rz_s(isite,imol,ibox) .le. cg_wall_position(1)) .or. (rz_s(isite,imol,ibox) .ge. cg_wall_position(2))) Then
	  			OVERLAP = .TRUE.
	  			RETURN
	  		End If

	  		! Double check
	  		If(rz_s(isite,imol,ibox) .lt. cg_wall_initz) Then
	  			engfield = engfield + 1.0d4
	  			CYCLE
	  		End if

	  		! Get distance to the surface
	  		zz = rz_s(isite,imol,ibox) - cg_wall_position(1)
	  		
	  		! Get bin number for energy grid
	  		ibin = CEILING((rz_s(isite,imol,ibox)-cg_wall_initz)/cg_wall_res)
	  		! case 1: ibin=0
	  		If (ibin .eq. 0) Then
	  			engfield = engfield + cg_wall_egrid(1)
	  			CYCLE
	  		End If
	  		! check if beyond the cutoff to avoid segmentation fault
	  		If ((ibin+1) .gt. cg_wall_maxnum) CYCLE
	  		
	  		! Do linear interpolation to calculate the energy [K]
	  		engfield = engfield + (zz-cg_wall_dgrid(ibin))*&
	  			&((cg_wall_egrid(ibin+1)-cg_wall_egrid(ibin))/(cg_wall_dgrid(ibin+1)-cg_wall_dgrid(ibin))) + &
	  			& cg_wall_egrid(ibin)
	  			
	  			
	  	! End loop over all sites on imol	
	  	End Do

	  	! Interactions with TOP walls
	  	! Loop over all sites on imol
	  	Do isite = 1, n_sites(itype)

	  		! Double check
	  		If((box(3,ibox)-rz_s(isite,imol,ibox)) .lt. cg_wall_initz) Then
	  			engfield = engfield + 1.0d4
	  			CYCLE
	  		End if

	  		! Get distance to the surface
	  		zz =  cg_wall_position(2) - rz_s(isite,imol,ibox) 
	  		 
	  		! Get bin number for energy grid
	  		ibin = CEILING((box(3,ibox)-rz_s(isite,imol,ibox)-cg_wall_initz)/cg_wall_res)
	  		! case 1: ibin=0
	  		If (ibin .eq. 0) Then
	  			engfield = engfield + cg_wall_egrid(1)
	  			CYCLE
	  		End If
	  		! check if beyond the cutoff to avoid segmentation fault
	  		If ((ibin+1) .gt. cg_wall_maxnum) CYCLE
	  		
	  		! Do linear interpolation to calculate the energy [K]
	  		engfield = engfield + (zz-cg_wall_dgrid(ibin))*&
	  			&((cg_wall_egrid(ibin+1)-cg_wall_egrid(ibin))/(cg_wall_dgrid(ibin+1)-cg_wall_dgrid(ibin))) + &
	  			& cg_wall_egrid(ibin)
	  			
	  			
	  	! End loop over all sites on imol	
	  	End Do


	  ! Coarse-grained wall potential with piecewise fluid-fluid interaction (from pre-calculated potential file)
	  Else If (field_type .eq. CG_WALL_FFPW) Then

	  	! Interactions with BOTTOM walls
	  	! Loop over all sites on imol
	  	Do isite = 1, n_sites(itype)

	  		! Check if isite is inside the wall
	  		If ((rz_s(isite,imol,ibox) .le. cg_wall_position(1)) .or. (rz_s(isite,imol,ibox) .ge. cg_wall_position(2))) Then
	  			OVERLAP = .TRUE.
	  			RETURN
	  		End If

	  		! Double check
	  		If(rz_s(isite,imol,ibox) .lt. cg_wall_initz) Then
	  			engfield = engfield + 1.0d4
	  			CYCLE
	  		End if

	  		! Get distance to the surface
	  		zz = rz_s(isite,imol,ibox) - cg_wall_position(1)
	  		
	  		! Get bin number for energy grid
	  		ibin = CEILING((rz_s(isite,imol,ibox)-cg_wall_initz)/cg_wall_res)
	  		! case 1: ibin=0
	  		If (ibin .eq. 0) Then
	  			engfield = engfield + cg_wall_egrid(1)
	  			CYCLE
	  		End If
	  		! check if beyond the cutoff to avoid segmentation fault
	  		If ((ibin+1) .gt. cg_wall_maxnum) CYCLE
	  		
	  		! Do linear interpolation to calculate the energy [K]
	  		engfield = engfield + (zz-cg_wall_dgrid(ibin))*&
	  			&((cg_wall_egrid(ibin+1)-cg_wall_egrid(ibin))/(cg_wall_dgrid(ibin+1)-cg_wall_dgrid(ibin))) + &
	  			& cg_wall_egrid(ibin)
	  			
	  			
	  	! End loop over all sites on imol	
	  	End Do

	  	! Interactions with TOP walls
	  	! Loop over all sites on imol
	  	Do isite = 1, n_sites(itype)

	  		! Double check
	  		If((box(3,ibox)-rz_s(isite,imol,ibox)) .lt. cg_wall_initz) Then
	  			engfield = engfield + 1.0d4
	  			CYCLE
	  		End if

	  		! Get distance to the surface
	  		zz =  cg_wall_position(2) - rz_s(isite,imol,ibox) 
	  		 
	  		! Get bin number for energy grid
	  		ibin = CEILING((box(3,ibox)-rz_s(isite,imol,ibox)-cg_wall_initz)/cg_wall_res)
	  		! case 1: ibin=0
	  		If (ibin .eq. 0) Then
	  			engfield = engfield + cg_wall_egrid(1)
	  			CYCLE
	  		End If
	  		! check if beyond the cutoff to avoid segmentation fault
	  		If ((ibin+1) .gt. cg_wall_maxnum) CYCLE
	  		
	  		! Do linear interpolation to calculate the energy [K]
	  		engfield = engfield + (zz-cg_wall_dgrid(ibin))*&
	  			&((cg_wall_egrid(ibin+1)-cg_wall_egrid(ibin))/(cg_wall_dgrid(ibin+1)-cg_wall_dgrid(ibin))) + &
	  			& cg_wall_egrid(ibin)
	  			
	  			
	  	! End loop over all sites on imol	
	  	End Do

	  ! Coarse-grained wall potential with cosine function mimicing the geometric roughness of the surface
	  Else If (field_type .eq. CG_WALL_COS) Then

	  	! Interactions with BOTTOM walls
	  	! Loop over all sites on imol
	  	Do isite = 1, n_sites(itype)

	  		! Get isite type 
	  		isitetype = site_type(isite,itype)

	  		! Check if isite is inside the wall
	  		If ((rz_s(isite,imol,ibox) .le. cg_wall_position(1)) .or. (rz_s(isite,imol,ibox) .ge. cg_wall_position(2))) Then
	  			OVERLAP = .TRUE.
	  			RETURN
	  		End If

	  		! Double check
	  		If(rz_s(isite,imol,ibox) .lt. cg_wall_initz) Then
	  			engfield = engfield + 1.0d4
	  			CYCLE
	  		End if

	  		! Get distance to the surface
	  		zz = rz_s(isite,imol,ibox) - cg_wall_position(1)

	  		! If larger than zlimit, no need to worry about cosine-like surface
	  		zlimit = 2*cg_wall_cosap + 0.5d0*sigma(isitetype,itype,isitetype,itype)
	  		If (zz .ge. zlimit) Then
	  			Continue
	  		Else if (zz .le. (cg_wall_cosap*COS(two_Pi*rx_s(isite,imol,ibox)/cg_wall_cospe)+cg_wall_cosap)) Then
	  			OVERLAP = .TRUE.
	  			RETURN
	  		Else 
	  			! Locate rough x-position of the site by dividing x-direction into bins of width 0.05 A
	  			xbin = CEILING(rx_s(isite,imol,ibox)/0.05d0)

	  			! Loop over neighbor points
	  			! 60 was chosen because to cover left/right range of 3 A, which is enough for Ar (diameter=3.405 A)
	  			Do cbin = 1, 60

	  				! Check left
	  				! x-position of the point on cosine function
	  				rxcos = DBLE((xbin-cbin))*0.05d0

	  				! x-distance between site and point on cosine function
	  				rxx = rx_s(isite,imol,ibox) - rxcos 

	  				! Check left part of the moelcular site
	  				rzz = zz - (cg_wall_cosap*COS(two_Pi*rxcos/cg_wall_cospe)+cg_wall_cosap)

	  				! Calculate distance square
	  				rsq = rxx**2+rzz**2

	  				If (rsq .ge. 0.25d0*sigmasq(isitetype,itype,isitetype,itype)) Then
	  					Continue
	  				Else
	  					OVERLAP = .TRUE.
	  					RETURN
	  				End If

	  				! Check right 
	  				! x-position of the point on cosine function
	  				rxcos = DBLE((xbin+cbin-1))*0.05d0

	  				! x-distance between site and point on cosine function
	  				rxx = rx_s(isite,imol,ibox) - rxcos 

	  				! Check left part of the moelcular site
	  				rzz = zz - (cg_wall_cosap*COS(two_Pi*rxcos/cg_wall_cospe)+cg_wall_cosap)

	  				! Calculate distance square
	  				rsq = rxx**2+rzz**2

	  				If (rsq .ge. 0.25d0*sigmasq(isitetype,itype,isitetype,itype)) Then
	  					Continue
	  				Else
	  					OVERLAP = .TRUE.
	  					RETURN
	  				End If
	  			End Do		

	  		! End checking contact with cosine-like surface
	  		End If

	  		! Get bin number for energy grid
	  		ibin = CEILING((rz_s(isite,imol,ibox)-cg_wall_initz)/cg_wall_res)
	  		! case 1: ibin=0
	  		If (ibin .eq. 0) Then
	  			engfield = engfield + cg_wall_egrid(1)
	  			CYCLE
	  		End If
	  		! check if beyond the cutoff to avoid segmentation fault
	  		If ((ibin+1) .gt. cg_wall_maxnum) CYCLE
	  		
	  		! Do linear interpolation to calculate the energy [K]
	  		engfield = engfield + (zz-cg_wall_dgrid(ibin))*&
	  			&((cg_wall_egrid(ibin+1)-cg_wall_egrid(ibin))/(cg_wall_dgrid(ibin+1)-cg_wall_dgrid(ibin))) + &
	  			& cg_wall_egrid(ibin)
	  			
	  			
	  	! End loop over all sites on imol	
	  	End Do

	  	! Interactions with TOP walls
	  	! Loop over all sites on imol
	  	Do isite = 1, n_sites(itype)

	  		! Double check
	  		If((box(3,ibox)-rz_s(isite,imol,ibox)) .lt. cg_wall_initz) Then
	  			engfield = engfield + 1.0d4
	  			CYCLE
	  		End if

	  		! Get distance to the surface
	  		zz =  cg_wall_position(2) - rz_s(isite,imol,ibox) 

	  		! If larger than zlimit, no need to worry about cosine-like surface
	  		If (zz .ge. zlimit) Then
	  			Continue
	  		Else if (zz .le. (cg_wall_cosap*COS(two_Pi*rx_s(isite,imol,ibox)/cg_wall_cospe)+cg_wall_cosap)) Then
	  			OVERLAP = .TRUE.
	  			RETURN
	  		Else 
	  			! Locate rough x-position of the site by dividing x-direction into bins of width 0.05 A
	  			xbin = CEILING(rx_s(isite,imol,ibox)/0.05d0)

	  			! Loop over neighbor points
	  			! 60 was chosen to cover left/right range of 3 A, which is enough for Ar (diameter=3.405 A)
	  			Do cbin = 1, 60

	  				! Check left
	  				! x-position of the point on cosine function
	  				rxcos = DBLE((xbin-cbin))*0.05d0

	  				! x-distance between site and point on cosine function
	  				rxx = rx_s(isite,imol,ibox) - rxcos 

	  				! Check left part of the moelcular site
	  				rzz = zz - (cg_wall_cosap*COS(two_Pi*rxcos/cg_wall_cospe)+cg_wall_cosap)

	  				! Calculate distance square
	  				rsq = rxx**2+rzz**2

	  				If (rsq .ge. 0.25d0*sigmasq(isitetype,itype,isitetype,itype)) Then
	  					Continue
	  				Else
	  					OVERLAP = .TRUE.
	  					RETURN
	  				End If

	  				! Check right 
	  				! x-position of the point on cosine function
	  				rxcos = DBLE((xbin+cbin-1))*0.05d0

	  				! x-distance between site and point on cosine function
	  				rxx = rx_s(isite,imol,ibox) - rxcos 

	  				! Check left part of the moelcular site
	  				rzz = zz - (cg_wall_cosap*COS(two_Pi*rxcos/cg_wall_cospe)+cg_wall_cosap)

	  				! Calculate distance square
	  				rsq = rxx**2+rzz**2

	  				If (rsq .ge. 0.25d0*sigmasq(isitetype,itype,isitetype,itype)) Then
	  					Continue
	  				Else
	  					OVERLAP = .TRUE.
	  					RETURN
	  				End If
	  			End Do		

	  		! End checking contact with cosine-like surface
	  		End If
	  		 
	  		! Get bin number for energy grid
	  		ibin = CEILING((box(3,ibox)-rz_s(isite,imol,ibox)-cg_wall_initz)/cg_wall_res)
	  		! case 1: ibin=0
	  		If (ibin .eq. 0) Then
	  			engfield = engfield + cg_wall_egrid(1)
	  			CYCLE
	  		End If
	  		! check if beyond the cutoff to avoid segmentation fault
	  		If ((ibin+1) .gt. cg_wall_maxnum) CYCLE
	  		
	  		! Do linear interpolation to calculate the energy [K]
	  		engfield = engfield + (zz-cg_wall_dgrid(ibin))*&
	  			&((cg_wall_egrid(ibin+1)-cg_wall_egrid(ibin))/(cg_wall_dgrid(ibin+1)-cg_wall_dgrid(ibin))) + &
	  			& cg_wall_egrid(ibin)
	  			
	  			
	  	! End loop over all sites on imol	
	  	End Do

	  ! Coarse-grained wall potential with accessible volume/length from full strcture
	  ! 3-4-2018, now the defects is at the middle of x-direction. and there is one defects
	  Else If (field_type .eq. CG_WALL_STRC) Then

	  	! Interactions with BOTTOM walls
	  	! Loop over all sites on imol
	  	Do isite = 1, n_sites(itype)

	  		! Get isite type 
	  		isitetype = site_type(isite,itype)

	  		! Check if isite is inside the wall
	  		If ((rz_s(isite,imol,ibox) .le. cg_wall_position(1)) .or. (rz_s(isite,imol,ibox) .ge. cg_wall_position(2))) Then
	  			OVERLAP = .TRUE.
	  			RETURN
	  		End If

	  		! Double check
	  		If(rz_s(isite,imol,ibox) .lt. cg_wall_initz) Then
	  			engfield = engfield + 1.0d4
	  			CYCLE
	  		End if

	  		! Get distance to the surface
	  		zz = rz_s(isite,imol,ibox) - cg_wall_position(1)

	  		! Check with lower bound of cg_wall_strc
	  		If(zz .lt. cg_strc_dgrid(1)) Then
	  			OVERLAP = .TRUE.
	  			RETURN
	  		End if

	  		! Check if overlap with strcture
	  		rxx = dABS(rx_s(isite,imol,ibox)-0.5d0*box(1,ibox))
	  		! Get bin number for accessible length
	  		zbin = CEILING((zz-cg_strc_dgrid(1))/cg_wall_res)
	  		! case 1: zbin=0
	  		If (zbin .eq. 0) Then
				If (rxx .le. cg_strc_dacc(1)) Then
					goto 100
				Else
					OVERLAP = .TRUE.
	  				RETURN
				End If
	  		End If
	  		! check if beyond the cutoff to avoid segmentation fault
	  		If ((zbin+1) .gt. cg_strc_maxnum) goto 100
	  		
	  		! Do linear interpolation to calculate the accessible volume/length [A]
	  		dacc = (zz-cg_strc_dgrid(zbin))*&
	  			&((cg_strc_dacc(zbin+1)-cg_strc_dacc(zbin))/(cg_strc_dgrid(zbin+1)-cg_strc_dgrid(zbin))) + &
	  			& cg_strc_dacc(zbin)

	  		If (rxx .le. dacc) Then
	  			goto 100
	  		Else 
	  			OVERLAP = .TRUE.
	  			RETURN
	  		End If

100			Continue 

	  		! Get bin number for energy grid
	  		ibin = CEILING((rz_s(isite,imol,ibox)-cg_wall_initz)/cg_wall_res)
	  		! case 1: ibin=0
	  		If (ibin .eq. 0) Then
	  			engfield = engfield + cg_wall_egrid(1)
	  			CYCLE
	  		End If
	  		! check if beyond the cutoff to avoid segmentation fault
	  		If ((ibin+1) .gt. cg_wall_maxnum) CYCLE
	  		
	  		! Do linear interpolation to calculate the energy [K]
	  		engfield = engfield + (zz-cg_wall_dgrid(ibin))*&
	  			&((cg_wall_egrid(ibin+1)-cg_wall_egrid(ibin))/(cg_wall_dgrid(ibin+1)-cg_wall_dgrid(ibin))) + &
	  			& cg_wall_egrid(ibin)
	  			
	  	
	  	End do		


	  	! Interactions with TOP walls
	  	! Loop over all sites on imol
	  	Do isite = 1, n_sites(itype)

	  		! Double check
	  		If((box(3,ibox)-rz_s(isite,imol,ibox)) .lt. cg_wall_initz) Then
	  			engfield = engfield + 1.0d4
	  			CYCLE
	  		End if

	  		! Get distance to the surface
	  		zz =  cg_wall_position(2) - rz_s(isite,imol,ibox) 

	  		! Check with lower bound of cg_wall_strc
	  		If(zz .lt. cg_strc_dgrid(1)) Then
	  			OVERLAP = .TRUE.
	  			RETURN
	  		End if

	  		! Check if overlap with strcture
	  		rxx = dABS(rx_s(isite,imol,ibox)-0.5d0*box(1,ibox))
	  		! Get bin number for accessible length
	  		zbin = CEILING((zz-cg_strc_dgrid(1))/cg_wall_res)
	  		! case 1: zbin=0
	  		If (zbin .eq. 0) Then
				If (rxx .le. cg_strc_dacc(1)) Then
					goto 200
				Else
					OVERLAP = .TRUE.
	  				RETURN
				End If
	  		End If
	  		! check if beyond the cutoff to avoid segmentation fault
	  		If ((zbin+1) .gt. cg_strc_maxnum) goto 200
	  		
	  		! Do linear interpolation to calculate the accessible volume/length [A]
	  		dacc = (zz-cg_strc_dgrid(zbin))*&
	  			&((cg_strc_dacc(zbin+1)-cg_strc_dacc(zbin))/(cg_strc_dgrid(zbin+1)-cg_strc_dgrid(zbin))) + &
	  			& cg_strc_dacc(zbin)

	  		If (rxx .le. dacc) Then
	  			goto 200
	  		Else 
	  			OVERLAP = .TRUE.
	  			RETURN
	  		End If

200	  		Continue
	  		 
	  		! Get bin number for energy grid
	  		ibin = CEILING((box(3,ibox)-rz_s(isite,imol,ibox)-cg_wall_initz)/cg_wall_res)
	  		! case 1: ibin=0
	  		If (ibin .eq. 0) Then
	  			engfield = engfield + cg_wall_egrid(1)
	  			CYCLE
	  		End If
	  		! check if beyond the cutoff to avoid segmentation fault
	  		If ((ibin+1) .gt. cg_wall_maxnum) CYCLE
	  		
	  		! Do linear interpolation to calculate the energy [K]
	  		engfield = engfield + (zz-cg_wall_dgrid(ibin))*&
	  			&((cg_wall_egrid(ibin+1)-cg_wall_egrid(ibin))/(cg_wall_dgrid(ibin+1)-cg_wall_dgrid(ibin))) + &
	  			& cg_wall_egrid(ibin)
	  			
	  			
	  	! End loop over all sites on imol	
	  	End Do


	  Else 

	  	Write(*,*) 'ENG_FIELD: INVALID FIELD_TYPE. ONLY HARD WALL TYPE HAS BEEN ADDED NOW'
	  	STOP
	  End If
	  
	  	Return
	  
	  End Subroutine 







	  	! 	!!!! 3-7-2018: New code for 'cg_wall_strc' with modfied potential for region 1
	  	! 	!! 3-7-2018, defect is 'region 1' with modified potential (cg_wall1.in),
	  	! 	! 		  perfect is 'region 2' with original potential (cg_wall2.in)
	  	! 	! Check if overlap with strcture
	  	! 	rxx = dABS(rx_s(isite,imol,ibox)-0.5d0*box(1,ibox))
	  	! 	! Get bin number for accessible length
	  	! 	zbin = CEILING((zz-cg_strc_dgrid(1))/cg_wall_res)
	  	! 	! check if beyond the cutoff to avoid segmentation fault
	  	! 	If ((zbin+1) .le. cg_strc_maxnum) Then
		  ! 		! case 1: zbin=0
		  ! 		If (zbin .eq. 0) Then
				! 	If (rxx .le. cg_strc_dacc(1)) Then
				! 		! Get bin number for energy grid
				!   		ibin = CEILING((rz_s(isite,imol,ibox)-cg_wall_initz)/cg_wall_res)
				!   		! case 1: ibin=0
				!   		If (ibin .eq. 0) Then
				!   			engfield = engfield + cg_wall1_egrid(1)
				!   			CYCLE
				!   		End If
				!   		! check if beyond the cutoff to avoid segmentation fault
				!   		If ((ibin+1) .gt. cg_wall_maxnum) CYCLE
				  		
				!   		! Do linear interpolation to calculate the energy [K]
				!   		engfield = engfield + (zz-cg_wall1_dgrid(ibin))*&
				!   			&((cg_wall1_egrid(ibin+1)-cg_wall1_egrid(ibin))/(cg_wall1_dgrid(ibin+1)-cg_wall1_dgrid(ibin)))&
				!   			& + cg_wall1_egrid(ibin)
				! 	Else
				! 		OVERLAP = .TRUE.
		  ! 				RETURN
				! 	End If
		  ! 		End If
	  	

		  ! 		If (rxx .le. cg_strc_dacc(1)) Then
		  ! 			! Get bin number for energy grid
			 !  		ibin = CEILING((rz_s(isite,imol,ibox)-cg_wall_initz)/cg_wall_res)
			 !  		! case 1: ibin=0
			 !  		If (ibin .eq. 0) Then
			 !  			engfield = engfield + cg_wall1_egrid(1)
			 !  			CYCLE
			 !  		End If
			 !  		! check if beyond the cutoff to avoid segmentation fault
			 !  		If ((ibin+1) .gt. cg_wall_maxnum) CYCLE
			  		
			 !  		! Do linear interpolation to calculate the energy [K]
			 !  		engfield = engfield + (zz-cg_wall1_dgrid(ibin))*&
			 !  			&((cg_wall1_egrid(ibin+1)-cg_wall1_egrid(ibin))/(cg_wall1_dgrid(ibin+1)-cg_wall1_dgrid(ibin)))&
			 !  			& + cg_wall1_egrid(ibin)
		  ! 		Else 
		  ! 			OVERLAP = .TRUE.
		  ! 			RETURN
		  ! 		End If

		  ! 	Else 

		  ! 		If (rxx .le. cg_strc_dacc(1)) Then
		  ! 			! Get bin number for energy grid
			 !  		ibin = CEILING((rz_s(isite,imol,ibox)-cg_wall_initz)/cg_wall_res)
			 !  		! case 1: ibin=0
			 !  		If (ibin .eq. 0) Then
			 !  			engfield = engfield + cg_wall1_egrid(1)
			 !  			CYCLE
			 !  		End If
			 !  		! check if beyond the cutoff to avoid segmentation fault
			 !  		If ((ibin+1) .gt. cg_wall_maxnum) CYCLE
			  		
			 !  		! Do linear interpolation to calculate the energy [K]
			 !  		engfield = engfield + (zz-cg_wall1_dgrid(ibin))*&
			 !  			&((cg_wall1_egrid(ibin+1)-cg_wall1_egrid(ibin))/(cg_wall1_dgrid(ibin+1)-cg_wall1_dgrid(ibin)))&
			 !  			& + cg_wall1_egrid(ibin)
			 !  	Else 
			 !  		! Get bin number for energy grid
			 !  		ibin = CEILING((rz_s(isite,imol,ibox)-cg_wall_initz)/cg_wall_res)
			 !  		! case 1: ibin=0
			 !  		If (ibin .eq. 0) Then
			 !  			engfield = engfield + cg_wall_egrid(1)
			 !  			CYCLE
			 !  		End If
			 !  		! check if beyond the cutoff to avoid segmentation fault
			 !  		If ((ibin+1) .gt. cg_wall_maxnum) CYCLE
			  		
			 !  		! Do linear interpolation to calculate the energy [K]
			 !  		engfield = engfield + (zz-cg_wall_dgrid(ibin))*&
			 !  			&((cg_wall_egrid(ibin+1)-cg_wall_egrid(ibin))/(cg_wall_dgrid(ibin+1)-cg_wall_dgrid(ibin))) + &
			 !  			& cg_wall_egrid(ibin)
				  			
		  ! 		End If
		  ! 	End if

	  		
	  	! ! End loop over all sites on imol	
	  	! End Do		
	  			


	  	! ! Interactions with TOP walls
	  	! ! Loop over all sites on imol
	  	! Do isite = 1, n_sites(itype)

	  	! 	! Double check
	  	! 	If((box(3,ibox)-rz_s(isite,imol,ibox)) .lt. cg_wall_initz) Then
	  	! 		engfield = engfield + 1.0d4
	  	! 		CYCLE
	  	! 	End if

	  	! 	! Get distance to the surface
	  	! 	zz =  cg_wall_position(2) - rz_s(isite,imol,ibox) 

	  	! 	! Check with lower bound of cg_wall_strc
	  	! 	If(zz .lt. cg_strc_dgrid(1)) Then
	  	! 		OVERLAP = .TRUE.
	  	! 		RETURN
	  	! 	End if

	  	! 	! Check if overlap with strcture
	  	! 	rxx = dABS(rx_s(isite,imol,ibox)-0.5d0*box(1,ibox))
	  	! 	! Get bin number for accessible length
	  	! 	zbin = CEILING((zz-cg_strc_dgrid(1))/cg_wall_res)
	  	! 	! check if beyond the cutoff to avoid segmentation fault
	  	! 	If ((zbin+1) .le. cg_strc_maxnum) Then
		  ! 		! case 1: zbin=0
		  ! 		If (zbin .eq. 0) Then
				! 	If (rxx .le. cg_strc_dacc(1)) Then
				! 		! Get bin number for energy grid
				!   		ibin = CEILING((box(3,ibox)-rz_s(isite,imol,ibox)-cg_wall_initz)/cg_wall_res)
				!   		! case 1: ibin=0
				!   		If (ibin .eq. 0) Then
				!   			engfield = engfield + cg_wall1_egrid(1)
				!   			CYCLE
				!   		End If
				!   		! check if beyond the cutoff to avoid segmentation fault
				!   		If ((ibin+1) .gt. cg_wall_maxnum) CYCLE
				  		
				!   		! Do linear interpolation to calculate the energy [K]
				!   		engfield = engfield + (zz-cg_wall1_dgrid(ibin))*&
				!   			&((cg_wall1_egrid(ibin+1)-cg_wall1_egrid(ibin))/(cg_wall1_dgrid(ibin+1)-cg_wall1_dgrid(ibin)))&
				!   			& + cg_wall1_egrid(ibin)
				! 	Else
				! 		OVERLAP = .TRUE.
		  ! 				RETURN
				! 	End If
		  ! 		End If
	  	

		  ! 		If (rxx .le. cg_strc_dacc(1)) Then
		  ! 			! Get bin number for energy grid
			 !  		ibin = CEILING((box(3,ibox)-rz_s(isite,imol,ibox)-cg_wall_initz)/cg_wall_res)
			 !  		! case 1: ibin=0
			 !  		If (ibin .eq. 0) Then
			 !  			engfield = engfield + cg_wall1_egrid(1)
			 !  			CYCLE
			 !  		End If
			 !  		! check if beyond the cutoff to avoid segmentation fault
			 !  		If ((ibin+1) .gt. cg_wall_maxnum) CYCLE
			  		
			 !  		! Do linear interpolation to calculate the energy [K]
			 !  		engfield = engfield + (zz-cg_wall1_dgrid(ibin))*&
			 !  			&((cg_wall1_egrid(ibin+1)-cg_wall1_egrid(ibin))/(cg_wall1_dgrid(ibin+1)-cg_wall1_dgrid(ibin)))&
			 !  			& + cg_wall1_egrid(ibin)
		  ! 		Else 
		  ! 			OVERLAP = .TRUE.
		  ! 			RETURN
		  ! 		End If

		  ! 	Else 

		  ! 		If (rxx .le. cg_strc_dacc(1)) Then
		  ! 			! Get bin number for energy grid
			 !  		ibin = CEILING((box(3,ibox)-rz_s(isite,imol,ibox)-cg_wall_initz)/cg_wall_res)
			 !  		! case 1: ibin=0
			 !  		If (ibin .eq. 0) Then
			 !  			engfield = engfield + cg_wall1_egrid(1)
			 !  			CYCLE
			 !  		End If
			 !  		! check if beyond the cutoff to avoid segmentation fault
			 !  		If ((ibin+1) .gt. cg_wall_maxnum) CYCLE
			  		
			 !  		! Do linear interpolation to calculate the energy [K]
			 !  		engfield = engfield + (zz-cg_wall1_dgrid(ibin))*&
			 !  			&((cg_wall1_egrid(ibin+1)-cg_wall1_egrid(ibin))/(cg_wall1_dgrid(ibin+1)-cg_wall1_dgrid(ibin)))&
			 !  			& + cg_wall1_egrid(ibin)
			 !  	Else 
			 !  		! Get bin number for energy grid
			 !  		ibin = CEILING((box(3,ibox)-rz_s(isite,imol,ibox)-cg_wall_initz)/cg_wall_res)
			 !  		! case 1: ibin=0
			 !  		If (ibin .eq. 0) Then
			 !  			engfield = engfield + cg_wall_egrid(1)
			 !  			CYCLE
			 !  		End If
			 !  		! check if beyond the cutoff to avoid segmentation fault
			 !  		If ((ibin+1) .gt. cg_wall_maxnum) CYCLE
			  		
			 !  		! Do linear interpolation to calculate the energy [K]
			 !  		engfield = engfield + (zz-cg_wall_dgrid(ibin))*&
			 !  			&((cg_wall_egrid(ibin+1)-cg_wall_egrid(ibin))/(cg_wall_dgrid(ibin+1)-cg_wall_dgrid(ibin))) + &
			 !  			& cg_wall_egrid(ibin)
				  			
		  ! 		End If
		  ! 	End if
	  			
	  			
	  	! ! End loop over all sites on imol	
	  	! End Do
	  	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!












