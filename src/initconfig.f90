! ==========================================================
! This subroutine is used to set initial configuration of each 
! molecule type in each box
! Created on 12-12-2016 by Kaihang Shi
! Last modified on 8-16-2017:
!   Fixed Read initial configuration from old files part.
! 	Add vol_pore(ibox) for slit pore geometry
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Molecule needs to be read from 'coords.in' always start from
! imol = 1
! ==========================================================


	  Subroutine initconfig(initlattice)

	  Use global

	  IMPLICIT NONE

	  ! Passed
	  Integer, Dimension(3,n_box_max) :: initlattice

	  ! Local
	  Integer :: ibox, itype, jtype, isite, jsite, imol, idirec, jdirec
	  Integer :: ninsert(n_mol_types), imolstart, ntotal, nwrap
	  Integer :: nx, ny, nz
	  Character (Len = 64) :: site_name
	  Double Precision :: space, z_increase, minbox, initx, inity, initz
	  Double Precision :: q1i, q2i, q3i, q4i

	  ! Variables used to read from initconfig.in file
	  Integer :: num, n_sites_tot



	  ! Set initial configuration for each box
	  Do ibox = 1, n_box


	  	! Check linit
	  	If (linit(ibox)) Then
	  		! linit = .true.
	  		goto 300
	  	Else
			! linit = .false.
	  		! Read initial configuration from old files
	  		!!!!! Needs modification every time when reading files
	  		!!!! Check carefully when use this 
	  		!---------------------------
	  		! 3-28-2018: Now read in 'initconfig.in' works for fully atomistic model
	  		! But ONLY works for spherical adsorbate (like Argon)
	  		!---------------------------

	  		Open(Unit=FILE_OLDCONFIG,FILE='initconfig.in',Status='Old',Access='Sequential',Action= 'Read')

	  		Read(FILE_OLDCONFIG,*) n_sites_tot
	  		Read(FILE_OLDCONFIG,*) n_mol_tot(ibox), (box(idirec,ibox), idirec=1,3)

	  		! Check the validity of cut-off radius
      		If (r_cut .gt. MINVAL(box(:,:)/2.0d0)) Then
      			Write(*,*) 'INITCONFIG: R_CUT LARGER THAN HALF OF MINBOX LENGTH'
      			STOP
      		End If

      		! Calculate volume of the system
	  		If (lfield .and. ((field_type .eq. STEELE) .or. (field_type .eq. STEELE_SLIT_PORE) .OR. &
      				& (field_type .eq. STEELE_SLIT_FINITEX))) Then
	  			vol_pore(ibox) = box(1,ibox)*box(2,ibox)*(box(3,ibox)-2*steele_position(1))
	  			vol(ibox) = box(1,ibox)*box(2,ibox)*box(3,ibox)

	  		Else if (lfield .and. ((field_type .eq. HARD_WALL) .OR. (field_type .eq. HARD_SLIT_FINITEX))) Then
	  			vol_pore(ibox) = box(1,ibox)*box(2,ibox)*(box(3,ibox)-2*wall_radius)
	  			vol(ibox) = box(1,ibox)*box(2,ibox)*box(3,ibox)

	  		Else if (lfield .and. ((field_type .eq. CG_WALL) .or. (field_type .eq. CG_WALL_FFPW) &
	  				& .or. (field_type .eq. CG_WALL_COS) .or. (field_type .eq. CG_WALL_STRC))) Then
	  			vol_pore(ibox) = box(1,ibox)*box(2,ibox)*(box(3,ibox)-2*cg_wall_position(1))
	  			vol(ibox) = box(1,ibox)*box(2,ibox)*box(3,ibox)

	  		Else
	  			vol(ibox) = box(1,ibox)*box(2,ibox)*box(3,ibox)
	  			vol_pore(ibox) = vol(ibox)

	  		End if
	  		
	  		! Print initial volume of the system
  			Write(*,'(A,I0,A,F20.4)') 'Initial volume of the box ',ibox, ' is: ', vol(ibox)
  			Write(*,'(A,I0,A,F20.4)') 'Initial pore volume of the box ',ibox, ' is: ', vol_pore(ibox)
  			Write(*,*) '!--- Warning: pore volume should be equal to the total system volume in NON-slit-like pore. ---!'


  			If (initstyle(1,ibox) .eq. 'coords') Then

  				! Assume external structure is type 1
  				itype = 1
  				! Read in external structure
				Do isite = 1, n_sites(itype)

					Read(FILE_OLDCONFIG,*) site_name, &
		  					&rx_s(isite,1,ibox), ry_s(isite,1,ibox), rz_s(isite,1,ibox)

		  			! Set internal coordinates
					rx_i(isite,itype) = rx_s(isite,1,ibox)
					ry_i(isite,itype) = ry_s(isite,1,ibox)
					rz_i(isite,itype) = rz_s(isite,1,ibox)

					! Initialize site_type
					site_type(isite,itype) = 1
						
				End Do

				n_mol(itype,ibox) = 1

				! Set molecule type. Assume imol=1 is external structure
	  			mol_type(1,ibox) = 1

	  			! set adsorbate as type 2
	  			itype = 2
	  			! Read in follow-up spherical adsorbate moelcules
	  			Do imol = 2, n_mol_tot(ibox)

	  				! Set imol type (only argon molecule in the box)
		  			mol_type(imol,ibox) = 2

		  			Read(FILE_OLDCONFIG,*) site_name, &
			  					&rx_s(1,imol,ibox), ry_s(1,imol,ibox), rz_s(1,imol,ibox)

			  		rx(imol,ibox) = rx_s(1,imol,ibox)
			  		ry(imol,ibox) = ry_s(1,imol,ibox)
			  		rz(imol,ibox) = rz_s(1,imol,ibox)

			  		! Initialize site_type
				  	site_type(1,itype) = 1

	  			End Do

	  			n_mol(itype,ibox) = n_mol_tot(ibox)-1

	  		! No external structure
	  		Else 
	  			! Read in sites information (1 atoms per argon molecule)
		  		Do imol = 1, n_mol_tot(ibox)

		  			! Set imol type (only argon molecule in the box)
		  			mol_type(imol,ibox) = 1

		  			Read(FILE_OLDCONFIG,*) site_name, &
			  					&rx_s(1,imol,ibox), ry_s(1,imol,ibox), rz_s(1,imol,ibox)

			  		rx(imol,ibox) = rx_s(1,imol,ibox)
			  		ry(imol,ibox) = ry_s(1,imol,ibox)
			  		rz(imol,ibox) = rz_s(1,imol,ibox)

			  		! Initialize site_type
			  		site_type(1,1) = 1

		  		End Do

		  		n_mol(1,ibox) = n_mol_tot(ibox)

  				
  			End If


	  		Close(FILE_OLDCONFIG)

	  		Write(*,'(A,I0,A)') 'Box ', ibox, ' has successfully read from old configuration'
	  		CYCLE

	  	End If



	  		!-------------------------
	  		! Read in NIST water config
	  		! ------------------------
	  		! Read in box dimension
	  		! Open(Unit=FILE_OLDCONFIG,FILE='initconfig.in',Status='Old',Access='Sequential',Action= 'Read')
	  		! Read(FILE_OLDCONFIG,*) (box(idirec,ibox), idirec=1,3)
	  		! Read(FILE_OLDCONFIG,*) n_mol_tot(ibox)
	  		! vol(ibox) = box(1,ibox)*box(2,ibox)*box(3,ibox)

	  		! Write(*,'(A,I0,A,F8.3,3X,F8.3,3X,F8.3)') 'Box ', ibox, &
     !  		& ' Dimensions: ', box(1,ibox), box(2,ibox), box(3,ibox)

	  		! ! Read in sites information (3 atoms per water molecule)
 	  	! 	Do imol = 1, n_mol_tot(ibox)
 	  	! 		! Set imol type (only water molecule in the box)
 	  	! 		mol_type(imol,ibox) = 1
				
 	  	! 		rx(imol,ibox) = 0.0d0
 	  	! 		ry(imol,ibox) = 0.0d0
 	  	! 		rz(imol,ibox) = 0.0d0

 	  	! 		Do isite = 1, n_sites(1)

 				! 	! NIST WATER
 				! 	Read(FILE_OLDCONFIG,*) num, rx_s(isite,imol,ibox), ry_s(isite,imol,ibox), &
 	  	! 			& rz_s(isite,imol,ibox), site_name

 	  	! 			! Initialize site_type
 		  ! 			site_type(isite,1) = -1

 		  ! 			! Set site type
   		! 			Do jsite = 1, n_site_types(1)
   		! 				If (Trim(site_name) .eq. site_type_name(jsite,1)) site_type(isite,1)=jsite
   		! 			End Do

  			! 		! Throw an error if no site type is found
   		! 			If (site_type(isite,1) .eq. -1) Then
   		! 				Write(*,'(2A)') 'FATAL ERROR: FAILED TO READ FROM OLD_CONFIG FOR ', site_name
   		! 				STOP
   		! 			End If

   		! 			rx(imol,ibox) = rx(imol,ibox) + rx_s(isite,imol,ibox)
   		! 			ry(imol,ibox) = ry(imol,ibox) + ry_s(isite,imol,ibox)
   		! 			rz(imol,ibox) = rz(imol,ibox) + rz_s(isite,imol,ibox)

 	  	! 		End Do

 	  	! 		! Calculate COM position
 	  	! 		rx(imol,ibox) = rx(imol,ibox)/3.0d0
 	  	! 		ry(imol,ibox) = ry(imol,ibox)/3.0d0
 	  	! 		rz(imol,ibox) = rz(imol,ibox)/3.0d0

	  			
 	  	! 	End Do
 	  	! 	Close(FILE_OLDCONFIG)

 	  	! 	Write(*,'(A,I0,A)') 'Box ', ibox, ' has successfully read from old configuration'
 	  	! 	CYCLE

 	  	! End If
	  	


300		Continue

	  	! Initialize local variables for each box
	  	imolstart = 1
	  	ntotal = 0
	  	ninsert(:) = 1
	  	space = 0.0d0
	  	z_increase = 0.0d0
	  	initx = 0.0d0
	  	inity = 0.0d0
	  	initz = 0.0d0

	  	! Initialize external structure flag (assume no external structure initially)
	  	ext_struc(ibox) = .false.

	  	! Check the validity of cut-off radius
      	If (r_cut .gt. MINVAL(box(:,:)/2.0d0)) Then
      		Write(*,*) 'INITCONFIG: R_CUT LARGER THAN HALF OF MINBOX LENGTH'
      		STOP
      	End If

	  	! Calculate box volume
  	    If (lfield .and. ((field_type .eq. STEELE) .or. (field_type .eq. STEELE_SLIT_PORE) .OR. &
  	    		& (field_type .eq. STEELE_SLIT_FINITEX))) Then
  			vol_pore(ibox) = box(1,ibox)*box(2,ibox)*(box(3,ibox)-2*steele_position(1))
  			vol(ibox) = box(1,ibox)*box(2,ibox)*box(3,ibox)
  		Else if (lfield .and. ((field_type .eq. HARD_WALL) .OR. (field_type .eq. HARD_SLIT_FINITEX))) Then
  			vol_pore(ibox) = box(1,ibox)*box(2,ibox)*(box(3,ibox)-2*wall_radius)
  			vol(ibox) = box(1,ibox)*box(2,ibox)*box(3,ibox)
  		Else if (lfield .and. ((field_type .eq. CG_WALL) .or. (field_type .eq. CG_WALL_FFPW) &
  				& .or. (field_type .eq. CG_WALL_COS) .or. (field_type .eq. CG_WALL_STRC))) Then
  			vol_pore(ibox) = box(1,ibox)*box(2,ibox)*(box(3,ibox)-2*cg_wall_position(1))
  			vol(ibox) = box(1,ibox)*box(2,ibox)*box(3,ibox)
  		Else
  			vol(ibox) = box(1,ibox)*box(2,ibox)*box(3,ibox)
  			vol_pore(ibox) = vol(ibox)
  		End if

  		! Print to screen
  		Write(*,'(A,I0,A,F20.4)') 'Initial volume of the box ',ibox, ' is: ', vol(ibox)
  		Write(*,'(A,I0,A,F20.4)') 'Initial pore volume of the box ',ibox, ' is: ', vol_pore(ibox)
  		Write(*,*) '!--- Warning: pore volume should be equal to the total system volume in NON-slit-like pore. ---!'

	  	! Determine the simple cubic lattice parameter
	  	! Initialize minbox value
	  	minbox = box(1,ibox)
	  	jdirec = 1
	  	! Find the minimum value of box length
	  	! For the moment, we assume cubic box, this step is trivial
	  	Do idirec = 2,3
	  		If (minbox .gt. box(idirec,ibox)) Then
	  			minbox = box(idirec,ibox)
	  			jdirec = idirec
	  		End If		
	  	End Do

	  	! Check validity of cutoff radius
	  	If (r_cut .gt. 0.5d0*minbox) Then
	  		Write(*,*) 'FATAL ERROR: CUTOFF RAIDUS IS TOO LARGE (r_cut > 1/2*box_length)'
	  		STOP
	  	End If

	  	! Calculate space parameter, we leave 0.1A space on each boundary
!	  	space = (box(jdirec,ibox)-0.2d0) / (initlattice(jdirec,ibox)-1)
		space = 2.0d0

	  	! Initialize the position of the first molecule of this type
	  	! x-coordinate
	  	If (Mod(initlattice(1,ibox),2) .eq. 0) Then
	  		initx = box(1,ibox)/2 - space/2 - (initlattice(1,ibox)/2 - 1.0d0)*space 
	  	Else
	  		initx = box(1,ibox)/2 - ((initlattice(1,ibox)-1)/2)*space 
	  	End If
	  	If (initx .le. 0.0d0) Then
	  		Write(*,*) 'FATAL ERROR: BAD SETTING OF INIX OR BOX LENGTH'
	  		Write(*,*) 'SHOULD DECREASE INIX VALUE'
	  		STOP
	  	End If
		! y-coordinate
		If (Mod(initlattice(2,ibox),2) .eq. 0) Then
			inity = box(2,ibox)/2 - space/2 - (initlattice(2,ibox)/2 - 1.0d0)*space 
		Else
			inity = box(2,ibox)/2 - ((initlattice(2,ibox)-1)/2)*space 
		End If
		If (inity .le. 0.0d0) Then
			Write(*,*) 'FATAL ERROR: BAD SETTING OF INIY OR BOX LENGTH'
			Write(*,*) 'SHOULD DECREASE INIY VALUE'
			STOP
		End If
		! z-coordinate
		If (Mod(initlattice(3,ibox),2) .eq. 0) Then
			initz = box(3,ibox)/2 - space/2 - (initlattice(3,ibox)/2 - 1.0d0)*space 
		Else
			initz = box(3,ibox)/2 - ((initlattice(3,ibox)-1)/2)*space 
		End If
		If (initz .le. 0.0d0) Then
			Write(*,*) 'FATAL ERROR: BAD SETTING OF INIZ OR BOX LENGTH' 
			Write(*,*) 'SHOULD DECREASE INIZ VALUE'
			STOP
		End If

		! Choose the imol start number for 'simple_cubic' style
  		Do jtype = 1, n_mol_types
  			if(initstyle(jtype,ibox) .eq. 'coords') Then
  				imolstart = 2
  				ext_struc(ibox) = .true.
  			End if
  		End Do

	  	! Loop over each molecule type
	  	Do itype =1, n_mol_types 

	  		! Choose different initial styles
	  		Select Case (initstyle(itype,ibox))
	  			
		  		! Read initial configration from coords.in file
		  		Case('coords')

		  			Write(*,'(A,A,A,I0,A)') 'Reading molecule ', mol_type_name(itype), &
		  			&' in box ', ibox, ' from coords.in file'

		  			! Open coords.in file to read
		  			Open(Unit=FILE_COORDS,File='coords.in',Status='Old',Access='Sequential',Action= 'Read')

		  			! Loop over each molecule of itype 
		  			! 'coords' molecule always start from imol=1 and we assume ONLY ONE molecule needs to 
		  			! read from 'coords.in'
		  			Do imol = 1, n_mol(itype,ibox)

		  				! Loop over site and read site info
		  				Do isite = 1, n_sites(itype)
		  					Read(FILE_COORDS,*) site_name, &
		  					&rx_s(isite,imol,ibox), ry_s(isite,imol,ibox), rz_s(isite,imol,ibox)

		  					! Set internal coordinates
		  					rx_i(isite,itype) = rx_s(isite,imol,ibox)
		  					ry_i(isite,itype) = ry_s(isite,imol,ibox)
		  					rz_i(isite,itype) = rz_s(isite,imol,ibox)

		  					! Check site_name length
		  					If (Len(Trim(site_name)) .gt. 6) Then
		  						Write(*,*) &
		  						&'FATAL ERROR: SITE NAME IN coords.in EXCEEDS THE LENGTH OF 6 CHARACTERS'
		  						Write(*,*) &
		  						&'CHECK global.f90 TO CHANGE THE LENGTH OF sites_name OR USE SIMPLER SITE NAME'
		  						STOP
		  					End If

		  					! Initialize site_type
		  					site_type(isite,itype) = -1

		  					! Set site type
		  					Do jsite = 1, n_site_types(itype)
		  						If (Trim(site_name) .eq. site_type_name(jsite,itype)) site_type(isite,itype)=jsite
		  					End Do

		  					! Throw an error if no site type is found
		  					If (site_type(isite,itype) .eq. -1) Then
		  						Write(*,'(2A)') 'FATAL ERROR: NO SITE TYPE IS FOUND FOR ', site_name
		  						STOP
		  					End If
		  				End Do

		  				! Throw finish reading notice
		  				Write(*,'(3A)') 'Finish reading initial structure of molecule ', &
		  				&mol_type_name(itype), ' from coords.in file'

		  				! Set molecule type
		  				mol_type(imol,ibox) = itype
		  			End Do

		  			Close(FILE_COORDS)
		  			


		  		Case('simple_cubic')

		  			Write(*,'(A,A,A,I0,A)') 'Initializing molecule ', mol_type_name(itype), &
		  			&' in box ', ibox, ' using simple cubic lattice structure'

		  			! Generate random quaternion (orientation) for all molecules
!		  			Call rand_quat(q1i,q2i,q3i,q4i)
					q1i = 0
					q2i = -0.2d0
					q3i = 0
					q4i = 0.9798d0

		  			! Loop over molecules of itype
		  			Do nz = 1, initlattice(3,ibox)
		  				Do ny = 1, initlattice(2,ibox)
		  					Do nx = 1, initlattice(1,ibox)

		  						! See if enough number of itype molecule have been inserted
		  						If (ninsert(itype) .gt. n_mol(itype,ibox)) Then
		  							Goto 100
		  						End If

		  						! Set COM coordinates
		  						rx(imolstart+ntotal,ibox) = initx + space*(nx - 1)
		  						ry(imolstart+ntotal,ibox) = inity + space*(ny - 1)
		  						rz(imolstart+ntotal,ibox) = initz + space*(nz - 1) + z_increase

		  						! Set molecule type
		  						mol_type(imolstart+ntotal,ibox) = itype

		  						! Set the molecule's orientation
		  						q1(imolstart+ntotal,ibox) = q1i
		  						q2(imolstart+ntotal,ibox) = q2i
		  						q3(imolstart+ntotal,ibox) = q3i
		  						q4(imolstart+ntotal,ibox) = q4i

		  						! Compute the coordinates of the sites on the molecule
		  						! & set site types
		  						Call site_coords(imolstart+ntotal,ibox)

		  						! Increase the counter
		  						ntotal = ntotal + 1
		  						ninsert(itype) = ninsert(itype) + 1

		  					End Do		
		  				End Do		
		  			End Do

100		  			Continue
		  			
		  			! In order to avoid overlap of different types of molecules 
		  			z_increase = z_increase + 1.0d0


		  		! Random initial configuration (Added on Oct 14, 2019)
		  		Case('random')

		  			Write(*,'(A,A,A,I0,A)') 'Initializing molecule ', mol_type_name(itype), &
		  			&' in box ', ibox, ' using random configuration'

		  			! Generate random quaternion (orientation) for all molecules
!		  			Call rand_quat(q1i,q2i,q3i,q4i)
					q1i = 0
					q2i = -0.2d0
					q3i = 0
					q4i = 0.9798d0

		  			! Loop over molecules of itype
		  			Do nz = 1, initlattice(3,ibox)
		  				Do ny = 1, initlattice(2,ibox)
		  					Do nx = 1, initlattice(1,ibox)

		  						! See if enough number of itype molecule have been inserted
		  						If (ninsert(itype) .gt. n_mol(itype,ibox)) Then
		  							Goto 200
		  						End If

		  						! Set COM coordinates
		  						rx(imolstart+ntotal,ibox) = random(idum)*box(1,ibox)
		  						ry(imolstart+ntotal,ibox) = random(idum)*box(2,ibox)
		  						rz(imolstart+ntotal,ibox) = random(idum)*box(3,ibox)

		  						! Set molecule type
		  						mol_type(imolstart+ntotal,ibox) = itype

		  						! Set the molecule's orientation
		  						q1(imolstart+ntotal,ibox) = q1i
		  						q2(imolstart+ntotal,ibox) = q2i
		  						q3(imolstart+ntotal,ibox) = q3i
		  						q4(imolstart+ntotal,ibox) = q4i

		  						! Compute the coordinates of the sites on the molecule
		  						! & set site types
		  						Call site_coords(imolstart+ntotal,ibox)

		  						! Increase the counter
		  						ntotal = ntotal + 1
		  						ninsert(itype) = ninsert(itype) + 1

		  					End Do		
		  				End Do		
		  			End Do

200		  			Continue
		  			


		  		Case('none')
		  			Goto 500

		  		
		  		
		  		Case default
		  			Write(*,*) 'FATAL ERROR: INVALID INITIAL STYLE IN input.in'
		  			STOP
	  		
	  
	  		End Select

500	  		Continue

	  	End do 

	  	! Double check with n_mol_tot
	  	If ((imolstart+ntotal-1) .ne. n_mol_tot(ibox)) Then
	  		Write(*,'(A,I0)') 'FATAL ERROR: INSERTED MOLECULES ARE NOT EQUAL TO TOTOAL NUMBER OF MOLS IN BOX ', ibox
	  		Write(*,'(A,I0)') 'INSERTED NUMBER: ', (imolstart+ntotal-1)
	  		Write(*,'(A,I0)') 'TOTAL NUMBER OF MOLS: ', n_mol_tot(ibox)
	  		STOP
	  	End If


	  End do 

	  Write(*,*) 'Initialization of configurations is completed.'	


	  ! Added on 2-19-2019
	  ! Wrap all atoms into the central box 
	  ! This is mainly for reading output structures generated by other programs (NAMD, Amber etc.)
	  ! Loop over box
	  Do ibox = 1, n_box

	  	! initialize number of atoms already wrapped
	  	nwrap = 0

	  	! Loop over each molecule type
	  	Do itype =1, n_mol_types 

	  		! Choose different initial styles
	  		Select Case (initstyle(itype,ibox))
	  			
		  		! for external structure, directly wrap site coordinate
		  		Case('coords')

		  			! 'coords' molecule always start from imol=1 and we assume ONLY ONE molecule (n_mol =1)
		  			Do imol = 1, n_mol(itype,ibox)

		  				! Loop over site 
		  				Do isite = 1, n_sites(itype)

		  					! wrap into central box
		  					rx_s(isite,imol,ibox) = &
		  						& rx_s(isite,imol,ibox) - FLOOR(rx_s(isite,imol,ibox)/box(1,ibox))*box(1,ibox)
		  					ry_s(isite,imol,ibox) = &
		  						& ry_s(isite,imol,ibox) - FLOOR(ry_s(isite,imol,ibox)/box(2,ibox))*box(2,ibox)
		  					rz_s(isite,imol,ibox) = &
		  						& rz_s(isite,imol,ibox) - FLOOR(rz_s(isite,imol,ibox)/box(3,ibox))*box(3,ibox)

		  					! Set internal coordinates
		  					rx_i(isite,itype) = rx_s(isite,imol,ibox)
		  					ry_i(isite,itype) = ry_s(isite,imol,ibox)
		  					rz_i(isite,itype) = rz_s(isite,imol,ibox)

		  				End Do
		  			End Do

		  		! wrap center-of-mass for other types of molecules
		  		! for now, its simple molecules, like argon, tip3p water...
		  		Case default
 				
		  			! Read in molecules of itype
		  			Do imol = nwrap+1, nwrap+n_mol(itype,ibox)

		  				! Put coordinates into central box 
						rx(imol,ibox) = rx(imol,ibox) - FLOOR(rx(imol,ibox)/box(1,ibox))*box(1,ibox) 
						ry(imol,ibox) = ry(imol,ibox) - FLOOR(ry(imol,ibox)/box(2,ibox))*box(2,ibox) 
						rz(imol,ibox) = rz(imol,ibox) - FLOOR(rz(imol,ibox)/box(3,ibox))*box(3,ibox) 

					  	! Compute the new sites coordinates
					  	Call site_coords(imol,ibox)

		  			End Do

	  		End Select

	  		nwrap = nwrap + n_mol(itype,ibox)

	    End Do

	  	Write(*,'(A,I0,A)') 'Box ', ibox, ' has been successfully wrapped'

	  ! End loop over boxes
	  ! End wrapping
	  End Do


	  Return
	  
	  End Subroutine initconfig




















