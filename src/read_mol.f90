! ==========================================================
! This subroutine is used to read molecular information (mol.in) 
! for each type of molecule
! Created on 12-5-2016 by Kaihang Shi
! Last modified 12-10-2016
! ==========================================================


	  Subroutine read_mol

	  Use global

	  IMPLICIT NONE

	  ! Local
	  Integer :: itype, isite, isitetype, jtype, jsitetype, inonbond
	  Integer :: n_nonbond
	  Integer :: nsite(n_nonbond_max)
!	  Integer :: ierr
	  Character(Len=64) :: temp_name
	  ! Counter 
	  Integer :: nonbond_start


	  ! Initialize the number and list of dispersion sites and electrostatic sites
	  n_sites(:) = 0

	  ! Open mol.in file to read 
	  Open(Unit=FILE_MOL,File='mol.in',Status='Old',Access='Sequential',Action= 'Read')

	  ! Read in number of nonbonded types for the use of 'Explicit' mixing rules
	  ! format like that set in towhee forcefield file (towhee_ff)
	  ! Seems redundant
	  If (label_check(FILE_MOL,'Number of Nonbonded Types')) Then
	  	Read(FILE_MOL,*) n_nonbond

	  	! Check with n_nonbond_types set in input.in
	  	If (n_nonbond .ne. n_nonbond_types) Then
	  		Write(*,*) 'READ_MOL: NUMBER OF NONBONDED TYPES IS NOT EQUAL TO THAT SET IN input.in'
	  		Write(*,*) ' CHECK input.in TO CHANGE n_nonbond_types VALUE ACCORDINGLY'
	  		STOP
	  	Else
	  		Write(*,'(A,I0)') 'Number of nonbonded types: ', n_nonbond 
	  	End If 

	  Else
	  	Write(*,*) 'FATAL ERROR: FAILED TO READ NUMBER OF NONBONDED TYPES LABEL IN mol.in'
	  	STOP
	  End If

	  ! Initialize counter 
	  nonbond_start = 0
	  

	  ! Loop for each molecule type
	  Do itype = 1,n_mol_types

	  	! Read the name of the molecule
	  	If (label_check(FILE_MOL,'Molecule Name')) Then
	  		Read(FILE_MOL,*) temp_name

	  		! Check the name length
	  		If (Len(Trim(temp_name)) .gt. 6) Then
	  			Write(*,'(3A)') 'FATAL ERROR: THE NAME LENGTH OF MOLECULE ', Trim(temp_name), &
	  			&' in mol.in IS MORE THAN 6 CHARACTER'
	  			Write(*,*) 'CHECK global.f90 TO CHANGE THE LENGTH OF mol_type_name'
	  			STOP
	  		End If
	  	Else 
	  		Write(*,*) 'FATAL ERROR: FAILED TO READ MOLECULE NAME LABEL IN mol.in'
	  		STOP
	  	End if 

	  	! Set the name of the molecule 
	  	mol_type_name(itype) = Trim(temp_name)

	  	! Write molecule type No.
	  	Write(*,'(A)') 'No.          Molecule'
	  	Write(*,'(I0,12X,A)') itype, mol_type_name(itype)
	  	Write(*,'(A)') 'No.          Site'

	  	! Read the mass of the molecule
	  	If (label_check(FILE_MOL,'Mass [g/mol]')) Then
	  		Read(FILE_MOL,*) mol_mass(itype)
	  	Else 
	  		Write(*,*) 'FATAL ERROR: FAILED TO READ MOLECULE MASS LABEL IN mol.in'
	  		STOP
	  	End if

	  	! Read number of sites on the molecule
	  	If (label_check(FILE_MOL,'Number of Site Types')) Then
	  		Read(FILE_MOL,*) n_site_types(itype)
	  	Else
	  		Write(*,*) 'FATAL ERROR: FAILED TO READ NUMBER OF SITES LABEL IN mol.in'
	  		STOP
	  	End if

	  	! Check n_site_types should be smaller than n_nonbond_max
	  	If (n_site_types(itype) .gt. n_nonbond) Then
	  		Write(*,*) 'FATAL ERROR: n_site_types SHOULD BE SMALLER THAN n_nonbond_max'
	  		STOP
	  	End if 

	  	! Read number for each site types to get total number of sites on each molecule
	  	If (label_check(FILE_MOL,'Number for Each Site Type')) Then

	  		Read(FILE_MOL,*,END=100) (nsite(isitetype), isitetype = 1,n_site_types(itype))

100			Continue
	  		
	  		! Calculate total number of sites on each molecule type
	  		Do isitetype = 1, n_site_types(itype)
	  			n_sites(itype) = n_sites(itype) + nsite(isitetype)
	  		End Do

	  	Else
	  		Write(*,*) 'FATAL ERROR: FAILED TO READ NUMBER FOR EACH SITE TYPES LABEL'
	  		STOP
	  	End if


	  	! Check with max number of sites on each molecule set in input.in
	  	If (n_sites(itype) .gt. n_sites_max) Then
	  		Write(*,'(3A)') 'FATAL ERROR: MOLECULE ', Trim(mol_type_name(itype)), &
	  		&' HAS MORE THAN THE MAXIMUM NUMBER OF SITES SET IN THE input.in'
	  		STOP
	  	End if


	  	! Read the information for each molecular site
	  	Do isitetype = 1, n_site_types(itype)

	  		! Increment counter
	  		nonbond_start = nonbond_start + 1

	  		! Read in site name
	  		If (label_check(FILE_MOL,'Site Name')) Then
	  			Read(FILE_MOL,*) temp_name

	  			! Check site name length 
	  			If (Len(Trim(temp_name)) .gt. 6) Then
	  				Write(*,'(5A)') 'FATAL ERROR: THE NAME LENGTH OF SITE ', Trim(temp_name), ' IN MOLECULE ', &
	  				& mol_type_name(itype), ' in mol.in IS MORE THAN 6 CHARACTER'
	  				Write(*,*) 'CHECK global.f90 TO CHANGE THE LENGTH OF sites_name'
	  				STOP
	  			End If

	  			! Set site name
	  			site_type_name(isitetype,itype) = Trim(temp_name)

	  			! Write site name
	  			Write(*,'(I0,12X,A)') isitetype, site_type_name(isitetype,itype)

	  		Else
	  			Write(*,*) 'FATAL ERROR: FAILED TO READ SITE NAME LABEL IN mol.in'
	  			STOP
	  		End if

	  		! Read in potential parameters. 
	  		! sigma in the unit of [Angstrom], epsilon in the unit of [K]
	  		! point charge in the unit of [e]
	  		! Lorentz-Berthelot
	  		If (mix_rule .eq. LORENTZ_BERTHELOT) Then

		  		! Read in point charge in the unit of [e]
		  		If (label_check(FILE_MOL,'Charge [e]')) Then
		  			Read(FILE_MOL,*) q(isitetype,itype,isitetype,itype)
		  		Else
		  			Write(*,'(2A)') 'FATAL ERROR: FAILED TO READ POINT CHARGE FOR SITE: ', site_type_name(isitetype,itype)
		  			STOP
		  		End if

		  		! Read in vdW potential parameters 
	  			If (label_check(FILE_MOL,'Nonbonded Coefficients')) Then
	  				! Read in sigma value in Angstrom
	  				Read(FILE_MOL,*) sigma(isitetype,itype,isitetype,itype)
	  				! Read in epsilon value in Kelvin
	  				Read(FILE_MOL,*) epsilon(isitetype,itype,isitetype,itype)
	  			Else 
	  				Write(*,*) 'FATAL ERROR: FAILED TO READ NONBONDED COEFFICIENTS LABEL IN mol.in'
	  				STOP
	  			End If

	  		! Explicit
	  		Else if (mix_rule .eq. EXPLICIT) Then

	  			! Read in point charge in the unit of [e]
		  		If (label_check(FILE_MOL,'Charge [e]')) Then
		  			! Read in original charge on site
		  			Read(FILE_MOL,*) q(isitetype,itype,isitetype,itype)
		  			! Read in perturbated charge on site (for now is Ad-Surf)
		  			Read(FILE_MOL,*) qex(AS,isitetype,itype,isitetype,itype)
		  		Else
		  			Write(*,'(2A)') 'FATAL ERROR: FAILED TO READ POINT CHARGE FOR SITE: ', site_type_name(isitetype,itype)
		  			STOP
		  		End if

		  		! Read in vdW potential parameters
		  		! Will be processed for simulation use in initialize subroutine
				Do inonbond = 1, n_nonbond_types

				 	If (label_check(FILE_MOL,'Nonbonded Coefficients')) Then
	  					! Read in sigma value in Angstrom
	  					Read(FILE_MOL,*) sigma(isitetype,itype,inonbond,inonbond)
	  					! Read in epsilon value in Kelvin
	  					Read(FILE_MOL,*) epsilon(isitetype,itype,inonbond,inonbond)
	  				Else 
	  					Write(*,*) 'FATAL ERROR: FAILED TO READ NONBONDED COEFFICIENTS LABEL IN mol.in'
	  					STOP
	  				End If

				End Do 


			! Other mixing rules
	  		Else 
	  			Write(*,*) 'FATAL ERROR: INVALID MIX_RULE VALUE'
	  			Write(*,*) "ONLY 'Lorentz-Berthelot' or 'Explicit' IS VALID"
	  			STOP
	  		End If

	  	! End looping over site types
	  	End Do

	  	! 'Coords' molecule's internal coordinates will be set in initconfig.f90 subroutine
 	  	If (initstyle(itype,1) .ne. 'coords') Then
	  		! Read in internal coordiantes of each site on the molecule
		  	! Coordinates are those with respect to Center of Mass
		  	If (label_check(FILE_MOL,'Internal Coordinates (wrt COM)')) Then
		  		
		  		Do isite = 1, n_sites(itype)

		  			! Read in site name and internal coordinates
		  			Read(FILE_MOL,*) sites_name(isite,itype), rx_i(isite,itype), ry_i(isite,itype), rz_i(isite,itype)

		  		End Do

		  	Else 
		  		Write(*,*) 'FATAL ERROR: FAILED TO READ INTERNAL COORDINATES IN mol.in FILE'
		  		STOP
		  	End If
	  	End If
	  	

	  ! End loop over molecule types
	  End do

	  Close(FILE_MOL)


	  Return

	  End Subroutine read_mol















