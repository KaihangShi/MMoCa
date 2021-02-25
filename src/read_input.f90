! ==========================================================
! This subroutine is used to read input file input.in and other
! related files
! Created on Dec. 3rd, 2016 by Kaihang Shi
! Last modified on June, 2020
! ==========================================================



      Subroutine read_input

      Use global

      IMPLICIT NONE

      ! Local 
      Character(Len=50) :: ensemblstr, pot_name, mix_rule_name, ewldtype, fieldname, optionname, pressdef, calcdef
      Integer :: ibox, itype, jtype, idirec, imove, inonbond
      Integer :: initlattice(3,n_box_max), ntot(n_box_max)

      ! Initialize local variables
      ntot(:) = 1


      ! Open input.in to read input parameters
      Open(Unit=FILE_INPUT,File='input.in',Status='Old',Access='Sequential',Action= 'Read')
      Rewind(FILE_INPUT)

      ! Read the emsemble type (string)
      If (label_check(FILE_INPUT,'Ensemble')) Then
      	Read(FILE_INPUT,*) ensemblstr

      	! Set the ensemble type
      	If (Trim(ensemblstr) .eq. 'NVT') Then
      		ensmbl = ENS_NVT
      	Else if (Trim(ensemblstr) .eq. 'NPT') Then
      		ensmbl = ENS_NPT
            Else if (Trim(ensemblstr) .eq. 'uVT') Then
                  ensmbl = ENS_uVT
      	Else 
      		Write(*,*) 'READ_INPUT: INVALID ENSEMBLE TYPE SPECIFIED IN input.in'
      		STOP
      	End If
      	Write(*,'(2A)') 'Ensemble: ', Trim(ensemblstr)

      Else 
      	Write(*,*) 'READ_INPUT: LABEL_CHECK FAILED TO READ ENSEMBLE LABEL IN input.in'
      	STOP
      End if

      ! Read the temperature in the unit of K
      If (label_check(FILE_INPUT,'Temperature [K]')) Then
      	Read(FILE_INPUT,*) temp
      	Write(*,'(A,F7.2)') 'Temperature [K]: ', temp
      Else
      	Write(*,*) 'READ_INPUT: LABEL_CHECK FAILED TO READ TEMPERATURE LABEL IN input.in'
      	STOP
      End if

      ! Read the pressure for NPT ensemble in the unit of bar
      If (ensmbl .eq. ENS_NPT) Then
      	If (label_check(FILE_INPUT,'Pressure [bar]')) Then
      		Read(FILE_INPUT,*) press 
      		Write(*,'(A,F15.7)') 'Pressure [bar]: ', press
      	Else
      		Write(*,*) 'READ_INPUT: FAILED TO READ PRESSURE LABEL FOR NPT ENSEMBLE IN input.in'
      	    STOP
      	End if 
      Else
      	Call skip_lines(FILE_INPUT,2)
      End if  


      ! Read number of box
      If(label_check(FILE_INPUT,'Number of Boxes')) Then
      	Read(FILE_INPUT,*) n_box

      	! Check with max number of boxes pre-set in global.f90
      	If (n_box .gt. n_box_max) Then
      		Write(*,*) 'READ_INPUT: NUMBER OF BOXES EXCEEDS THE PRESET NUMBER'
      		STOP
      	End If

      	Write(*,'(A,I0)') 'Number of boxes: ', n_box
      Else 
      	Write(*,*) 'READ_INPUT: LABEL_CHECK FAILED TO READ NUMBER OF BOXES LABEL IN input.in'
      	STOP
      End if 

      ! Read number of molecular types
      If(label_check(FILE_INPUT,'Number of Molecule Types')) Then
      	Read(FILE_INPUT,*) n_mol_types

      	! Check with n_mol_types_max set in global.f90
      	If (n_mol_types .gt. n_mol_types_max) Then
      		Write(*,*) 'READ_INPUT: NUMBER OF MOLECULE TYPES SHOULD BE SMALLER THAN PRE-SET ONE'
      		Write(*,*) 'CHECK global.f90 TO RESET THE VALUE OF n_mol_types_max'
      		STOP
      	End If

      	Write(*,'(A,I0)') 'Number of molecule types: ', n_mol_types
      Else 
      	Write(*,*) 'READ_INPUT: LABEL_CHECK FAILED TO READ NUMBER OF MOLECULAR TYPES LABEL IN input.in'
      	STOP
      End if


      ! Read maximum number of sites on each molecule
      If(label_check(FILE_INPUT,'Maximum Number of Sites on the Molecule')) Then
      	Read(FILE_INPUT,*) n_sites_max
      	Write(*,'(A,I0)') 'Maximum number of sites on a molecule: ', n_sites_max
      Else 
      	Write(*,*) 'READ_INPUT: LABEL_CHECK FAILED TO READ MAX NUMBER OF SITES LABEL IN input.in'
      	STOP
      End if

      ! Read number of nonbonded types in the system
      If (label_check(FILE_INPUT,'Number of Nonbonded Types')) Then
            Read(FILE_INPUT,*) n_nonbond_types
      Else
            Write(*,*) 'READ_INPUT: FAILED TO READ IN NUMBER OF NONBONDED TYPES LABEL'
            STOP
      End if


      ! Allocate the box dimension 
      Call global_allocate('Box')

      ! Allocate system arrays (molecular topology, potential parameters)
      Call global_allocate('System')


      ! Read in predefined chemical potential for uVT ensemble in the unit of K
      If (ensmbl .eq. ENS_uVT) Then
            If (label_check(FILE_INPUT,'Chemical Potential [K]')) Then
                  Read(FILE_INPUT,*) (mu(itype), itype = 1, n_mol_types)
                  ! Dump to screen
                  Do itype = 1, n_mol_types
                        Write(*,'(A,I0,A,F15.7)') &
                              & 'Molecule Type: ',itype,'  Chemical Potential [K]: ',mu(itype)        
                  End Do 
            Else
                  Write(*,*) 'READ_INPUT: FAILED TO READ PRESSURE LABEL FOR NPT ENSEMBLE IN input.in'
                STOP
            End if 
      Else
            Call skip_lines(FILE_INPUT,2)
      End if  


      ! Initialize total number of molecules in each box
      n_mol_tot(:) = 0

      ! Read number of molecules of each type in each box
      ! Reading sequence of molecules should be in accordance with that in mol.in
      ! For N-fixed ensemble, this is N;  
      ! For uVT, this is the initial number of molecule in the box
      If (label_check(FILE_INPUT,'Number of Molecules')) Then

      	Do ibox = 1, n_box
      		Read(FILE_INPUT,*) (n_mol(jtype,ibox), jtype = 1, n_mol_types)
      		! Calculate total number of molecules in each box
      		Do jtype = 1, n_mol_types
      			n_mol_tot(ibox) = n_mol_tot(ibox) + n_mol(jtype,ibox)
      		End Do

      		! Print total number of molecules to the screen
      		Write(*,'(A,I0,A,I0)') 'Total number of molecules in box ', ibox, ': ',n_mol_tot(ibox)	 
        End Do
     
      Else
      	Write(*,*) 'READ_INPUT: FAILED TO READ NUMBER OF MOLECULES LABEL IN input.in'
      	STOP
      End if

      

      ! Read in potential form
      ! Other types of potential form could be added here
      If (label_check(FILE_INPUT, 'Potential Form')) Then
      	Read(FILE_INPUT,*) pot_name

      	! Set potential form
      	If (Trim(pot_name) .eq. 'Lennard-Jones') Then
      		potential = LENNARD_JONES
      	Else 
      		Write(*,*) 'READ_INPUT: INVALID POTENTIAL FORM IN input.in'
      		Write(*,*) "Only 'Lennard-Jones' is valid so far"
      		STOP
      	End if
!      	Write(*,'(2A)') 'Potential form: ', Trim(pot_name)

      Else 
      	Write(*,*) 'READ_INPUT: LABEL_CHECK FAILED TO READ POTENTIAL FORM LABEL IN input.in'
      	STOP
      End if

      ! Read in mixing rule for sigma and epsilon 
      ! EXPLICIT could be used to test combining rule in CS model
      If(label_check(FILE_INPUT,'Mixing Rules')) Then
      	Read(FILE_INPUT,*) mix_rule_name

      	! Set mixing rule
      	If (Trim(mix_rule_name) .eq. 'Lorentz-Berthelot') Then
      		mix_rule = LORENTZ_BERTHELOT
      	Else if (Trim(mix_rule_name) .eq. 'Explicit') Then
      		mix_rule = EXPLICIT 
      	Else 
      		Write(*,*) 'READ_INPUT: INVALID MIXING RULES IN input.in'
      		Write(*,*) "Either 'Lorentz-Berthelot' or 'Explicit'"
      		STOP
      	End if
      	Write(*,'(2A)') 'Mixing rules: ', Trim(mix_rule_name)

      Else 
      	Write(*,*) 'READ_INPUT: LABEL_CHECK FAILED TO READ MIXING RULES LABEL IN input.in'
      	STOP
      End if


      ! Set up tail correction flag on/off
      If (label_check(FILE_INPUT,'Tail Correction')) Then
      	Read(FILE_INPUT,*) ltailc
      Else
      	Write(*,*) 'READ_INPUT: FAILED TO READ TAIL CORRECTION LABEL IN input.in'
      	STOP
      End if

      
      ! Read in cut-off radius for all types of molecule
      If (label_check(FILE_INPUT,'r_cutoff [A]')) Then
      	Read(FILE_INPUT,*) r_cut
      	
      	r_cutsq = r_cut**2

      	Write(*,'(A,F7.3)') 'r_cutoff: ', r_cut
      Else
      	Write(*,*) 'READ_INPUT: FAILED TO READ CUT-OFF RADIUS LABEL IN input.in file'
      	STOP
      End if

      ! Read in inner cut-off radius to speed up LJ system
      If (label_check(FILE_INPUT,'r_min [A]')) Then
      	Read(FILE_INPUT,*) r_min
      	r_minsq = r_min**2

      	Write(*,'(A,F7.3)') 'r_min: ', r_min
      Else
      	Write(*,*) 'READ_INPUT: FAILED TO READ INNER CUT-OFF RAIDUS LABEL in input.in file'
      	STOP
      End if

      ! Read in External Field flag
      If (label_check(FILE_INPUT,'External Field')) Then
            Read(FILE_INPUT,*) lfield
      Else
            Write(*,*) 'READ_INPUT: FAILED TO READ IN EXTERNAL FIELD LABEL'
            STOP
      End if

            !Read in External Field type
            If (label_check(FILE_INPUT,'Field Type')) Then
                  Read(FILE_INPUT,*) fieldname

                  ! Set field type
                  If (Trim(fieldname) .eq. 'hard_wall') Then
                        field_type = HARD_WALL
                  Else if (Trim(fieldname) .eq. 'hard_slit_finitex') Then
                        field_type = HARD_SLIT_FINITEX          
                  Else if (Trim(fieldname) .eq. 'steele_wall') Then
                        field_type = STEELE
                  Else if (Trim(fieldname) .eq. 'steele_slit_pore') Then
                        field_type = STEELE_SLIT_PORE
                  Else if (Trim(fieldname) .eq. 'steele_slit_finitex') Then
                        field_type = STEELE_SLIT_FINITEX
                  Else if (Trim(fieldname) .eq. 'cg_wall') Then
                        field_type = CG_WALL
                  Else if (Trim(fieldname) .eq. 'cg_wall_ffpw') Then
                        field_type = CG_WALL_FFPW
                  Else if (Trim(fieldname) .eq. 'cg_wall_cos') Then
                        field_type = CG_WALL_COS
                  Else if (Trim(fieldname) .eq. 'cg_wall_strc') Then
                        field_type = CG_WALL_STRC                  
                  Else
                        Write(*,*) 'READ_INPUT: ONLY HARD_WALL FIELD TYPE HAS BEEN CODED'
                        STOP
                  End If
            Else
                  Write(*,*) 'READ_INPUT: FAILED TO READ IN FIELD TYPE LABEL'
                  STOP
            End if

            ! Read in external field parameters
            If (lfield) Then
                  ! Hard wall
                  If (field_type .eq. HARD_WALL) Then

                        Write(*,*) 'Hard wall is used...'

                        If (label_check(FILE_INPUT,'Half Thickness')) Then
                              ! Read in radius (half of thickness of hard wall)
                              Read(FILE_INPUT,*) wall_radius

                              Write(*,'(A,F6.2)') 'Hard wall half thickness of ', wall_radius
                              Write(*,*) 'Note: Center of hard wall is put on z=0'
                        Else
                              Write(*,*) 'READ_INPUT: FAILED TO READ IN HALF THICKNESS LABEL FOR HARD WALL'
                              STOP
                        End if

                  ! Hard slit wall with averaging region in the middle (compatible with Yun's model)
                  Else if (field_type .eq. HARD_SLIT_FINITEX) Then                       

                        If (label_check(FILE_INPUT,'Half Thickness')) Then
                              ! Read in radius (half of thickness of hard wall)
                              Read(FILE_INPUT,*) wall_radius

                              Write(*,'(A,F6.2)') 'Hard wall half thickness of ', wall_radius
                              Write(*,*) 'Note: Center of hard wall is put on z=0'
                        Else
                              Write(*,*) 'READ_INPUT: FAILED TO READ IN HALF THICKNESS LABEL FOR HARD WALL'
                              STOP
                        End if

                        ! Averaing region
                        If (label_check(FILE_INPUT,'Averaging Region')) Then
                              ! Read in lower and upper bound of averaging region
                              Read(FILE_INPUT,*) steele_avgx(1), steele_avgx(2)
                              ! Length of the averaging region in x-direction (added on 6-9-2020)
                              l_avgx = steele_avgx(2) - steele_avgx(1)
                              If (l_avgx .LT. r_cut) Then
                                    Write(*,*) 'read_input: Averaging region should be > = r_cut (HARD_SLIT_FINITEX)'
                                    STOP
                              End If
                        Else
                              Write(*,*) 'READ_INPUT: FAILED TO READ IN AVERAGING REGION LABEL '
                              STOP
                        End if

                        Write(*,*) 'Hard slit wall is used...'
                        Write(*,'(A,F10.4)') 'steele_avgx(1) = ', steele_avgx(1)
                        Write(*,'(A,F10.4)') 'steele_avgx(2) = ', steele_avgx(2)

                  ! 10-4-3 Steele wall
                  Else if (field_type .eq. STEELE) Then

                        Write(*,*) '10-4-3 Steele potential is used...'

                        ! Position of Steele wall 
                        ! There will also be a corresponding hard wall on the top of the box (z direction)
                        If (label_check(FILE_INPUT,'Position [A]')) Then
                              Read(FILE_INPUT,*) steele_position(1)
                              Write(*,'(A,F7.3)') 'Position: ', steele_position(1)

                        Else
                              Write(*,*) 'READ_INPUT: FAILED TO READ IN POSITION LABEL FOR STEELE WALL'
                              STOP
                        End if
                        ! Cutoff raidus for Steele wall
                        If (label_check(FILE_INPUT,'Cutoff [A]')) Then
                              Read(FILE_INPUT,*) steele_cut
                              Write(*,'(A,F7.3)') 'Cutoff: ', steele_cut
                        Else
                              Write(*,*) 'READ_INPUT: FAILED TO READ IN CUTOFF LABEL FOR STEELE WALL'
                              STOP
                        End if
                        ! Delta value 
                        If (label_check(FILE_INPUT,'Delta [A]')) Then
                              Read(FILE_INPUT,*) steele_delta
                              Write(*,'(A,F7.3)') 'Delta: ', steele_delta
                        Else
                              Write(*,*) 'READ_INPUT: FAILED TO READ IN DELTA LABEL FOR STEELE WALL'
                              STOP
                        End if
                        ! Solid density
                        If (label_check(FILE_INPUT,'Rho_s [A^-3]')) Then
                              Read(FILE_INPUT,*) steele_rhos
                              Write(*,'(A,F7.3)') 'Rho_s: ', steele_rhos
                        Else
                              Write(*,*) 'READ_INPUT: FAILED TO READ IN RHO_S LABEL FOR STEELE WALL'
                              STOP
                        End if
                        ! Steele interaction parameters
                        If (label_check(FILE_INPUT,'Interaction Parameters')) Then
                              
                              ! Loop all nonbonded types
                              Do inonbond = 1, n_nonbond_types

                                    ! Read in values 
                                    ! We use n_mol_types_max rank here only for storage of value
                                    ! Array will be processed later in initialize subroutine
                                    Read(FILE_INPUT,*) steele_site_name(inonbond), &
                                          & steele_sigmasf(inonbond,n_mol_types_max), &
                                          & steele_epsilonsf(inonbond,n_mol_types_max)  

!                                    Write(*,*) 'Interaction parameters (for test only): '
!                                    Write(*,'(A,2F7.3)') steele_site_name(inonbond), &
!                                          & steele_sigmasf(inonbond,n_mol_types_max), &
!                                          & steele_epsilonsf(inonbond,n_mol_types_max)  
                                    
                              End Do
                        Else
                              Write(*,*) 'READ_INPUT: FAILED TO READ IN INTERACTION PARAMETERS FOR STEELE WALL'
                              STOP
                        End if
                  
                  ! 10-4-3 Steele wall for slit pore geometry
                  Else if ((field_type .eq. STEELE_SLIT_PORE) .OR. (field_type .eq. STEELE_SLIT_FINITEX)) Then

                        Write(*,*) '10-4-3 Steele potential for slit pore geometry is used...'

                        ! Position of Steele wall (downside)
                        ! There will also be a corresponding steele wall on the top of the box (z direction)
                        If (label_check(FILE_INPUT,'Position [A]')) Then
                              Read(FILE_INPUT,*) steele_position(1)
!                              Write(*,'(A,F7.3)') 'Position for Downside Wall: ', steele_position(1)
                        Else
                              Write(*,*) 'READ_INPUT: FAILED TO READ IN POSITION LABEL FOR STEELE WALL'
                              STOP
                        End if
                        ! Cutoff raidus for Steele wall
                        If (label_check(FILE_INPUT,'Cutoff [A]')) Then
                              Read(FILE_INPUT,*) steele_cut
                              Write(*,'(A,F7.3)') 'Cutoff: ', steele_cut
                        Else
                              Write(*,*) 'READ_INPUT: FAILED TO READ IN CUTOFF LABEL FOR STEELE WALL'
                              STOP
                        End if
                        ! Delta value 
                        If (label_check(FILE_INPUT,'Delta [A]')) Then
                              Read(FILE_INPUT,*) steele_delta
                              Write(*,'(A,F7.3)') 'Delta: ', steele_delta
                        Else
                              Write(*,*) 'READ_INPUT: FAILED TO READ IN DELTA LABEL FOR STEELE WALL'
                              STOP
                        End if
                        ! Solid density
                        If (label_check(FILE_INPUT,'Rho_s [A^-3]')) Then
                              Read(FILE_INPUT,*) steele_rhos
                              Write(*,'(A,F7.3)') 'Rho_s: ', steele_rhos
                        Else
                              Write(*,*) 'READ_INPUT: FAILED TO READ IN RHO_S LABEL FOR STEELE WALL'
                              STOP
                        End if
                        ! Steele interaction parameters
                        If (label_check(FILE_INPUT,'Interaction Parameters')) Then
                              
                              ! Loop all nonbonded types
                              Do inonbond = 1, n_nonbond_types

                                    ! Read in values 
                                    ! We use n_mol_types_max rank here only for storage of value
                                    ! Array will be processed later in initialize subroutine
                                    Read(FILE_INPUT,*) steele_site_name(inonbond), &
                                          & steele_sigmasf(inonbond,n_mol_types_max), &
                                          & steele_epsilonsf(inonbond,n_mol_types_max)  

!                                    Write(*,*) 'Interaction parameters (for test only): '
!                                    Write(*,'(A,2F7.3)') steele_site_name(inonbond), &
!                                          & steele_sigmasf(inonbond,n_mol_types_max), &
!                                          & steele_epsilonsf(inonbond,n_mol_types_max)  
                                    
                              End Do
                        Else
                              Write(*,*) 'READ_INPUT: FAILED TO READ IN INTERACTION PARAMETERS FOR STEELE WALL'
                              STOP
                        End if

                  ! Coarse-grained fluid-solid potential (w(D)) from free energy averaged method
                  Else If (field_type .eq. CG_WALL) Then

                        Write(*,*) 'Coarse-grained wall is used...'
                        Write(*,*) 'Potential energy w(D) is pre-calculated using free energy averaged method (sf)'

                        If (label_check(FILE_INPUT,'Position [A]')) Then
                              ! Read in bottom position of cg wall
                              Read(FILE_INPUT,*) cg_wall_position(1)
                        Else
                              Write(*,*) 'READ_INPUT: FAILED TO READ IN POSITION LABEL FOR CG_WALL'
                              STOP
                        End if

                        ! Starting point for fluid-solid potential calculation
                        If (label_check(FILE_INPUT,'Initz [A]')) Then
                              Read(FILE_INPUT,*) cg_wall_initz
                        Else
                              Write(*,*) 'FATAL ERROR: FAILED TO READ IN INITZ LABEL FOR CG_WALL'
                              STOP
                        End if

                        ! Resolution of fluid-solid potential
                        If (label_check(FILE_INPUT,'Resolution [A]')) Then
                              Read(FILE_INPUT,*) cg_wall_res
                        Else
                              Write(*,*) 'READ_INPUT: FAILED TO READ IN RESOLUTION LABEL FOR CG_WALL '
                              STOP
                        End if


                  ! Coarse-grained fluid-solid potential with piecewise effective fluid-fluid interactions
                  Else If (field_type .eq. CG_WALL_FFPW) Then

                        Write(*,*) 'Coarse-grained wall with piecewise fluid-fluid interaction is used...'

                        If (label_check(FILE_INPUT,'Position [A]')) Then
                              ! Read in bottom position of cg wall
                              Read(FILE_INPUT,*) cg_wall_position(1)
                        Else
                              Write(*,*) 'READ_INPUT: FAILED TO READ IN POSITION LABEL FOR CG_WALL'
                              STOP
                        End if

                        ! Starting point for fluid-solid potential calculation
                        If (label_check(FILE_INPUT,'Initz [A]')) Then
                              Read(FILE_INPUT,*) cg_wall_initz
                        Else
                              Write(*,*) 'FATAL ERROR: FAILED TO READ IN INITZ LABEL FOR CG_WALL'
                              STOP
                        End if

                        ! Resolution of fluid-solid potential
                        If (label_check(FILE_INPUT,'Resolution [A]')) Then
                              Read(FILE_INPUT,*) cg_wall_res !, cg_ff_res
                        Else
                              Write(*,*) 'READ_INPUT: FAILED TO READ IN RESOLUTION LABEL FOR CG_WALL '
                              STOP
                        End if

                        ! Cutoff distance for the piecewise interaction
                        If (label_check(FILE_INPUT,'Range [A]')) Then
                              Read(FILE_INPUT,*) cg_ff_lob, cg_ff_upb
                        Else
                              Write(*,*) 'READ_INPUT: FAILED TO READ IN RANGE LABEL FOR CG_WALL_FFPW'
                              STOP
                        End if


                  ! Coarse-grained fluid-solid potential (w(D)) from free energy averaged method and use 
                  ! Cosine function to mimic the geometric roughness of the surface
                  Else If (field_type .eq. CG_WALL_COS) Then

                        Write(*,*) 'Coarse-grained wall with cosine function to mimic geometric roughness is used...'
                        Write(*,*) 'Potential energy w(D) is pre-calculated using free energy averaged method (sf)'

                        If (label_check(FILE_INPUT,'Position [A]')) Then
                              ! Read in bottom position of cg wall
                              Read(FILE_INPUT,*) cg_wall_position(1)
                        Else
                              Write(*,*) 'READ_INPUT: FAILED TO READ IN POSITION LABEL FOR CG_WALL_COS'
                              STOP
                        End if

                        ! Starting point for fluid-solid potential calculation
                        If (label_check(FILE_INPUT,'Initz [A]')) Then
                              Read(FILE_INPUT,*) cg_wall_initz
                        Else
                              Write(*,*) 'FATAL ERROR: FAILED TO READ IN INITZ LABEL FOR CG_WALL_COS'
                              STOP
                        End if

                        ! Resolution of fluid-solid potential
                        If (label_check(FILE_INPUT,'Resolution [A]')) Then
                              Read(FILE_INPUT,*) cg_wall_res
                        Else
                              Write(*,*) 'READ_INPUT: FAILED TO READ IN RESOLUTION LABEL FOR CG_WALL_COS '
                              STOP
                        End if

                        ! Parameters for Cosine function
                        If (label_check(FILE_INPUT,'Parameters')) Then
                              ! Read in number of periodicity and amplitude
                              Read(FILE_INPUT,*) cg_wall_cosnp, cg_wall_cosap
                        Else
                              Write(*,*) 'READ_INPUT: FAILED TO READ IN PARAMETERS LABEL FOR CG_WALL_COS '
                              STOP
                        End if
                  ! Read in accessible volume (length) from all atomistic structure
                  Else If (field_type .eq. CG_WALL_STRC) Then

                        Write(*,*) 'Coarse-grained wall is used...'
                        Write(*,*) 'Potential energy w(D) is pre-calculated using free energy averaged method (sf)'

                        If (label_check(FILE_INPUT,'Position [A]')) Then
                              ! Read in bottom position of cg wall
                              Read(FILE_INPUT,*) cg_wall_position(1)
                        Else
                              Write(*,*) 'READ_INPUT: FAILED TO READ IN POSITION LABEL FOR CG_WALL'
                              STOP
                        End if

                        ! Starting point for fluid-solid potential calculation
                        If (label_check(FILE_INPUT,'Initz [A]')) Then
                              Read(FILE_INPUT,*) cg_wall_initz
                        Else
                              Write(*,*) 'FATAL ERROR: FAILED TO READ IN INITZ LABEL FOR CG_WALL'
                              STOP
                        End if

                        ! Resolution of fluid-solid potential
                        If (label_check(FILE_INPUT,'Resolution [A]')) Then
                              Read(FILE_INPUT,*) cg_wall_res
                        Else
                              Write(*,*) 'READ_INPUT: FAILED TO READ IN RESOLUTION LABEL FOR CG_WALL '
                              STOP
                        End if

                  End If

            Else 
                  ! Skip to next Label line
                  If (field_type .eq. HARD_WALL) Then
                        ! skip next 2 lines
                        Call skip_lines(FILE_INPUT,2)
                  Else If (field_type .eq. HARD_SLIT_FINITEX) Then
                        ! skip next 4 lines
                        Call skip_lines(FILE_INPUT,4)
                  Else If ((field_type .eq. STEELE) .or. (field_type .eq. STEELE_SLIT_PORE) .OR. &
                        & (field_type .eq. STEELE_SLIT_FINITEX)) Then
                        ! skip next (9+n_nonbond_types) lines
                        Call skip_lines(FILE_INPUT,9+n_nonbond_types)
                  Else If ((field_type .eq. CG_WALL) .or. (field_type .eq. CG_WALL_STRC)) Then
                        ! skip next 6 lines
                        Call skip_lines(FILE_INPUT,6)
                  Else If (field_type .eq. CG_WALL_FFPW) Then
                        ! skip next 8 lines
                        Call skip_lines(FILE_INPUT,8)
                  Else If (field_type .eq. CG_WALL_COS) Then
                        ! skip next 8 lines
                        Call skip_lines(FILE_INPUT,8)
                  End If
                  
            End If



      ! Check if needs to initialize the configuration from the scractch
      ! if .false., read the initial config from external file for ibox
      ! if .true., create initial configurations by initstyle
      If (label_check(FILE_INPUT,'linit')) Then
      	Do ibox = 1, n_box
      		Read(FILE_INPUT,*) linit(ibox)
      	End Do
      	
      Else
      	Write(*,*) 'READ_INPUT: FAILED TO READ linit flag in input.in'
      	STOP
      End if

      ! Set up initial style for each type of molecule
      ! 'coords' - initial config read from coords.in
      ! 'simple_cubic' - initial config by simple cubic lattice
      ! 'random' - initial random config
      ! 'none' - no such type of molecule in ibox
      If (label_check(FILE_INPUT,'Initial Style')) Then

      	! Loop over simulation boxes
      	Do ibox = 1, n_box
      		Read(FILE_INPUT,*) (initstyle(itype,ibox), itype = 1,n_mol_types)
      	End Do
      
      Else
      	Write(*,*) 'READ_INPUT: FAILED TO READ INITIAL STYLE LABEL IN input.in'
      	STOP
      End if

      ! Read in initial number of molecules in each direction in each box
      ! The product of inix*iniy*iniz should be .gt. or .eq. totoal number of molecules in each box
      If (label_check(FILE_INPUT,'Inix Iniy Iniz')) Then
      	
      	! Loop over simulation boxes
      	Do ibox = 1, n_box
      		Read(FILE_INPUT,*) (initlattice(idirec,ibox), idirec = 1, 3)

      		! Check if the product of inix*iniy*iniz .gt. or .eq. n_mol_tot
      		Do idirec = 1,3
      			ntot(ibox) = ntot(ibox)*initlattice(idirec,ibox)
      		End Do
      		If(ntot(ibox) .lt. n_mol_tot(ibox)) Then
      			Write(*,*) 'READ_INPUT: INSUFFICIENT MOLECULES TO CREATE INITIAL CONFIGURATION'
      			Write(*,*) 'INIX*INIY*INIZ SHOULD BE .GT. OR .EQ. TOTAL NUMBER OF MOLECULES IN EACH BOX'
      			STOP
      		End if 
      	End Do

      Else
      	Write(*,*) 'READ_INPUT: FAILED TO READ INIX INIY INIZ LABEL IN input.in'
      	STOP
      End if


      ! Read box dimensions
      If (label_check(FILE_INPUT,'Box Dimension')) Then

      	! Loop over each box
      	Do ibox = 1, n_box
      		Read(FILE_INPUT,*) (box(idirec,ibox), idirec = 1,3)

      		if(linit(ibox)) &
      		& Write(*,'(A,I0,A,F8.3,3X,F8.3,3X,F8.3)') 'Box ', ibox, &
      		& ' Dimensions: ', box(1,ibox), box(2,ibox), box(3,ibox)
      	End do 

      	! Check if cubic box in NPT ensemble
      	If (ensmbl .eq. ENS_NPT) Then
      		If ((box(1,1) .ne. box(2,1)) .or. (box(1,1) .ne. box(3,1)) .or. (box(2,1) .ne. box(3,1))) Then
      			Write(*,*) 'READ_INPUT: NPT SIMULATION REQURES CUBIC BOX'
      			STOP
      		End If
      	End If

      Else
      	Write(*,*) 'READ_INPUT: FAILED TO READ BOX DIMENSION LABEL IN input.in'
      	STOP
      End if


      ! Read in the information for each molecule type (see read_mol.f90)
      ! Also write molecular and site No. to screen
      Call read_mol

      ! Initialize potential parameters 
      Call initialize('Potential')

      ! Initialize Steele potential parameters
      if(lfield .and. ((field_type .eq. STEELE) .or. (field_type .eq. STEELE_SLIT_PORE) .or. &
            & (field_type .eq. STEELE_SLIT_FINITEX))) Call initialize('Steele')

      ! Initialize coarse-grained wall potentials
      if(lfield .and. (field_type .eq. CG_WALL)) Call initialize('CG_WALL')
      ! Initialize coarse-grained wall potential with cosine function mimicing the geometric roughness
      if(lfield .and. (field_type .eq. CG_WALL_COS)) Call initialize('CG_WALL')

      ! Initialize coarse-grained wall potentials with piecewise fluid-fluid interactions
      if(lfield .and. (field_type .eq. CG_WALL_FFPW)) Call initialize('CG_WALL_FFPW')

      ! Initialize coarse-grained wall potentials with accessible volume/length read from full structure
      if(lfield .and. (field_type .eq. CG_WALL_STRC)) Call initialize('CG_WALL_STRC')

      ! Construct initial configurations for each box
      Call initconfig(initlattice)


      ! Read in Ewald Summation flag 
      If (label_check(FILE_INPUT,'Ewald Sum')) Then
      	Read(FILE_INPUT,*) lewld
      Else
      	Write(*,*) 'READ_INPUT: FAILED TO READ EWALD SUM LABEL IN input.in file'
      	STOP
      End If

      ! Read in Ewald sum parameters
        ! Read in Ewald style
        If (label_check(FILE_INPUT,'Ewald Style')) Then
              Read(FILE_INPUT,*) ewldtype

              ! Set Ewald Style
              If (Trim(ewldtype) .eq. 'ewald_fix_kmax') Then
                  ewld_style = ewald_fix_kmax
              Else
                  Write(*,*) 'READ_INPUT: INVALID EWALD STYLE'
                  STOP
              End If
        Else
              Write(*,*) 'READ_INPUT: FAILED TO READ IN EWALD STYLE LABEL'
              STOP
        End if

        ! Kmax-fixed scheme
        If (ewld_style .eq. ewald_fix_kmax) Then
              ! Read in kalp value (alpha = kalp/MINVAL(box(:,ibox)))
              If (label_check(FILE_INPUT,'kalp')) Then
                  Read(FILE_INPUT,*) kalp
              Else
                  Write(*,*) 'READ_INPUT: FAILED TO READ IN kalp VALUE FOR EWALD SUM USE'
                  STOP
              End If

              ! Read in k_max value
              If (label_check(FILE_INPUT,'kmax')) Then
              	! Loop over boxes
              	Do ibox = 1, n_box
              		Read(FILE_INPUT,*) (k_max(idirec,ibox), idirec = 1, 3)

              		! Print to screen
              		Write(*,'(A,I0,A,I0,3x,I0,3x,I0)') 'kmax for box ',ibox, ': ', &
              			& k_max(1,ibox), k_max(2,ibox), k_max(3,ibox)
              	End Do
                  
              Else
              	Write(*,*) 'READ_INPUT: FAILED TO READ IN KMAX VALUE'
              	STOP
              End if
        End If
	
      ! Read in Slab correction flag for slit pore model
      If (label_check(FILE_INPUT,'Slab Correction')) Then
            Read(FILE_INPUT,*) lslabc
      Else
            Write(*,*) 'READ_INPUT: FAILED TO READ IN SLAB CORRECTION LABEL'
            STOP
      End if  

      ! Setup and allocate ewald sum parameters
      if(lewld) Call set_ewld(0,1)

      ! Set basic parameters for Ewald sum
      if(lewld) Call set_ewld(2,1)
      WRITE(*,'(A,F10.3)') 'Initial Alpha in Ewald Sum: ',alpha
      WRITE(*,'(A,F10.3)') 'Initial rcelect: ',rcelect


      ! Setup the k-vectors for each box
      If (lewld) Then
      	Do ibox = 1,n_box
      		Call set_ewld(1,ibox)
      		Write(*,'(A,I0)') 'Successfully set up k-vectors for box ',ibox
      	End Do
      End If


      ! Read in spherical cut-off for adsorbate-surface interaction feature
      If (label_check(FILE_INPUT,'Spherical cut-off for ad-surf interaction')) Then
      	READ(FILE_INPUT,*) lcoulsc
      	! Read in spherical cut-off raidus [A]
      	Read(FILE_INPUT,*) rscelect

      	! Check the validity of cut-off radius
      	If (rscelect .gt. MINVAL(box(:,:)/2.0d0)) Then
      		Write(*,*) 'READ_INPUT: rscelect LARGER THAN HALF OF MINBOX LENGTH'
      		STOP
      	End If

      	rscelectsq = rscelect**2

      	! Check consistency 
      	If (lcoulsc) Then
      		If (.not. lewld) Then
      			Write(*,*) 'FATAL ERROR: SPHERICAL CUT-OFF FOR AD-SURF INTERACTION &
      						& AND EWALD SUM SHOULD BOTH BE TURNED ON'
      			STOP
      		End If
      	End If

      	! Print to Screen
      	If (lcoulsc) Then
      		Write(*,*) 'Spherical cut-off for ad-surf interaction is applied'
      		Write(*,'(A,F5.2)') 'cut-off radius:', rscelect
      	End If

      Else
      	Write(*,*) 'READ_INPUT: FAILED TO READ IN SPHERICAL CUT-OFF FOR AD-SURF INTERACTION LAEBL'
      	STOP
      End if


      ! Read in MC run parameters
      ! Read in number of relax blocks
      If (label_check(FILE_INPUT,'Relax Blocks')) Then
            Read(FILE_INPUT,*) n_blocks_relax
            Write(*,'(A,I0)') 'Number of relax blocks: ', n_blocks_relax
      Else
            Write(*,*) 'READ_INPUT: FAILED TO READ RELAX BLOCKS LABEL FROM input.in'
            STOP
      End if

      ! Read in total number of simulation blocks
      If (label_check(FILE_INPUT,'Total Number of Blocks')) Then
            Read(FILE_INPUT,*) n_blocks_tot
        	Write(*,'(A,I0)') 'Total number of simulation blocks: ', n_blocks_tot
      Else
            Write(*,*) 'READ_INPUT: FAILED TO READ TOTAL NUMBER OF BLOCKS LABEL FROM input.in'
        	STOP
      End if

      ! Read in equilibrium blocks number
      If (label_check(FILE_INPUT,'Equilibrium Blocks')) Then
        	Read(FILE_INPUT,*) n_blocks_equil
        	Write(*,'(A,I0)') 'Number of equilibrium blocks: ', n_blocks_equil
      Else
        	Write(*,*) 'READ_INPUT: FAILED TO READ EQUILIBRIUM BLOCKS LABEL IN input.in'
        	STOP
      End if
        
      ! Read in block size
      If (label_check(FILE_INPUT,'Block Size')) Then
        	Read(FILE_INPUT,*) block_size
        	Write(*,'(A,I0)') 'Block size: ', block_size
      Else
            Write(*,*) 'READ_INPUT: FAILED TO READ BLOCK SIZE LABEL IN input.in'
        	STOP
      End if

!	Write(*,*) 'Code is sucessful so far'
!	STOP

      ! Allocate the block average array
      Call global_allocate('Block')



      ! Read in probabilities of selecting different MC move types
      If (label_check(FILE_INPUT,&
      	&'Probability of Selecting MC Move Types (Translation, Rotation, Transfer, Volume)')) Then

            ! For the first half part of equilibrium
            Read(FILE_INPUT,*) (move_type_prob(imove), imove=1,4)
            ! For the later half part of equilibrium
!            Read(FILE_INPUT,*) (move_type_prob_update(1,imove), imove=1,4)
            ! For production stage
!      	Read(FILE_INPUT,*) (move_type_prob_update(2,imove), imove=1,4)


		! Check for bad setting of probabilities
		! 1 for translational move, 2 for Rotational move, 3 for transfer move, 4 for volume change
		Do imove = 2,4 
			If (move_type_prob(imove) .lt. move_type_prob(imove-1)) Then
				Write(*,*) 'READ_INPUT: BAD SETTING OF move_type_prob'
				Write(*,*) '0 <= prob(translate) <= prob(rotate) <= prob(transfer) <=prob(volchange) = 1'
				STOP
			Else if (move_type_prob(4) .ne. 1.0d0) Then
				Write(*,*) 'READ_INPUT: move_type_prob for transfer move must be 1.0d0'
				STOP
			End If
		End Do
      	
      Else
      	Write(*,*) 'READ_INPUT: FAILED TO READ PROBABILITIES OF SELECTING MC MOVE TYPE LABEL'
      	STOP
      End if


      ! Read in probabilities of translational move for each molecule type
      If (label_check(FILE_INPUT,'Probability of Translational Move for Each Molecule Type')) Then
      	Read(FILE_INPUT,*) (trans_prob(itype), itype = 1, n_mol_types)

      	! Check for bad setting of probabilities
      	If (trans_prob(n_mol_types) .ne. 1.0d0) Then
      		Write(*,*) 'READ_INPUT: PROBABILITIES FOR LAST MOL TYPE MUST BE 1.0d0'
      		STOP
      	End If

      Else
      	Write(*,*) 'READ_INPUT: FAILED TO READ PROBABILITIES OF Translational MOVE FOR EACH MOLECULE TYPE LABEL'
      	STOP
      End if

      ! Read in probabilities of rotational move for each molecule type
      If (label_check(FILE_INPUT,'Probability of Rotational Move for Each Molecule Type')) Then
      	Read(FILE_INPUT,*) (rotat_prob(itype), itype = 1, n_mol_types)

      	! Check for bad setting of probabilities
      	If (rotat_prob(n_mol_types) .ne. 1.0d0) Then
      		Write(*,*) 'READ_INPUT: PROBABILITIES FOR LAST MOL TYPE MUST BE 1.0d0'
      		STOP
      	End If

      Else
      	Write(*,*) 'READ_INPUT: FAILED TO READ PROBABILITIES OF ROTATIONAL MOVE FOR EACH MOLECULE TYPE LABEL'
      	STOP
      End if

      ! Read in probabilities of transferral move for each molecule type
      If (label_check(FILE_INPUT,'Probability of Transferral Move for Each Molecule Type')) Then
            Read(FILE_INPUT,*) (transfer_prob(itype), itype = 1, n_mol_types)

            ! Check for bad setting of probabilities
            If (transfer_prob(n_mol_types) .ne. 1.0d0) Then
                  Write(*,*) 'READ_INPUT: PROBABILITIES FOR LAST MOL TYPE MUST BE 1.0d0'
                  STOP
            End If

      Else
            Write(*,*) 'READ_INPUT: FAILED TO READ PROBABILITIES OF TRANSFERRAL MOVE FOR EACH MOLECULE TYPE LABEL'
            STOP
      End if


      ! Read in Widom Insertion Method info
      If (label_check(FILE_INPUT,'Widom Insertion')) Then
      	Read(FILE_INPUT,*) lwdm
      Else
      	Write(*,*) 'READ_INPUT: FAILED TO READ IN WIDOM INSERTION LABEL IN input.in'
      	STOP
      End if

      ! Check ensemble
      If (lwdm) Then
      	if(ensmbl .ne. ENS_NPT) Then
      		Write(*,*) 'READ_INPUT: WIDOM METHOD ONLY FOR NPT ENSEMBLE NOW'
      		STOP
      	End If
      End If

      ! Read in Widom insertion frequency
      If (label_check(FILE_INPUT,'Insertion Frequency for Each Molecule Type')) Then
      	
      	! Loop over each box
      	Do ibox = 1, n_box
        	   Read(FILE_INPUT,*) (widom_freq(itype,ibox), itype = 1, n_mol_types)
      	End Do

      Else
      	Write(*,*) 'READ_INPUT: FAILED TO READ INSERTION FREQUENCY IN input.in'
      	STOP
      End if


      ! Read in other sampling options
      If (label_check(FILE_INPUT,'Sampling Options')) Then

      	! Initialize all the sampling options flag
      	lsampling = .false.
            lzdensity = .false.
      	lrdensity = .false.
      	lsurfex = .false.
            lthermopress_slit = .false.
            lvirialpress_slit = .false.
            lvirialpress = .false.
            lvirialpress_cylin = .false.
            llattconst = .false.
            lqst = .false.
            ldumpxyz = .false.
            dump_moltype = 0
            ldumpdensity = .false.
            ldumpenergy = .false.
            lwriterestart = .false. 
            lno_pbc = .false.
            lclist = .false.
            ldumpvir = .false.

      	! Continue reading sampling options until EOF
      	Do

      		! Read in options
      		Read(FILE_INPUT,*) optionname

      		! Set up z-density flag 
      		! Density profile of each molecule (COM) along z-axis 
      		If (Trim(optionname) .eq. 'z_density') Then

      			lzdensity = .true.
      			lsampling = .true.

      			! Read in bins to produce density profile 
      			Read(FILE_INPUT,*) zden_bins

                        ! Read in calculation frequency
                        Read(FILE_INPUT,*) zden_freq

      			! Print
      			Write(*,*) 'Z-Density will be calculated...'


                  ! Set up r-density flag 
                  ! Density profile of each molecule (COM) along contant-R surface in cylinder
                  Else If (Trim(optionname) .eq. 'r_density') Then

                        lrdensity = .true.
                        lsampling = .true.

                        ! Read in bins to produce density profile 
                        Read(FILE_INPUT,*) rden_bins

                        !Read in calculation frequency
                        Read(FILE_INPUT,*) rden_freq

                        ! Read in cutoff radius for calculation and force cutoff for cylindrical pressure tensor calculation
                        Read(FILE_INPUT,*) rden_cut

                        ! Print
                        Write(*,'(A,F6.3,A)') 'R-Density will be calculated within ',rden_cut,' Angstrom...' 


      		! Set up surface excess flag 
      		! Excess adsorption number of adsorbates
      		Else If (Trim(optionname) .eq. 'surface_excess') Then

      			lsurfex = .true.
      			lsampling = .true.

      			! Read in preset bulk density for each molecule type
      			! surfex_bulk for external structure is redundant, but necessary to input.
      			Read(FILE_INPUT,*) (surfex_bulk(itype), itype = 1, n_mol_types)

      			! Read in surface area and accessible volume
      			Read(FILE_INPUT,*) surfex_area, surfex_vol

      			! Print
      			Write(*,*) 'Surface excess adsorption number will be calculated...'


                  ! Calculate the pressure tensor for slit geometry from thermo route
                  Else If (Trim(optionname) .eq. 'thermo_press_slit') Then

                        lthermopress_slit = .true. 
                        lsampling = .true.

                        ! Read in volume change frequency
                        Read(FILE_INPUT,*) thermopress_slit_freq

                        ! Read in length of box change factor 
                        Read(FILE_INPUT,*) thermopress_slit_ratio

                        ! Read in the bins for producing pressure tensor
                        Read(FILE_INPUT,*) thermopress_slit_bins

                        ! Print
                        Write(*,*) 'Pressure tensor for the slit geometry will be calculated by thermodynamic route'
                        Write(*,'(A,I0)') 'Volume change frequency:', thermopress_slit_freq
                        Write(*,'(A,F5.3)') 'Length (dLx/Lx) change ratio:', thermopress_slit_ratio

                  ! Calculate the pressure tensor for slit geometry from virial route
                  Else If (Trim(optionname) .eq. 'virial_press_slit') Then

                        If (.not. lzdensity) Then
                              Write(*,*) 'READ_INPUT: lzdensity should be turned on to calculate kinetic part of the pressure'
                              STOP
                        End If

                        lvirialpress_slit = .true. 
                        lsampling = .true.

                        ! Read in the definition of integral contour for pressure tensor calculation and the contribution to be accounted
                        Read(FILE_INPUT,*) pressdef, calcdef

                        ! Set & Print type
                        If ((pressdef .eq. 'IK') .or. (pressdef .eq. 'Irving-Kirkwood')) Then
                              virialpress_ctype = 1
                              Write(*,*) 'Pressure tensor calculation ONLY by Irving-Kirkwood definition'
                        Else if ((pressdef .eq. 'H') .or. (pressdef .eq. 'Harasima')) Then
                              virialpress_ctype = 2
                              Write(*,*) 'Pressure tensor calculation ONLY by Harasima definition'
                        Else if (pressdef .eq. 'IK&H') Then
                              virialpress_ctype = 3
                              Write(*,*) 'Pressure tensor calculation by both IK and Harasima definition'
                        ! Added on May 29, 2019
                        Else if (pressdef .eq. 'H-VR') Then 
                              virialpress_ctype = 4
                              Write(*,*) 'Pressure tensor by the variation of the Harasima type (H-VR)'
                        ! Added on May 30, 2019
                        Else if (pressdef .eq. 'IK-VR') Then 
                              virialpress_ctype = 5
                              Write(*,*) 'Pressure tensor by the variation of the IK type (IK-VR)'
                        ! Added on May 30, 2019
                        Else if (pressdef .eq. 'ALL') Then 
                              virialpress_ctype = 6
                              Write(*,*) 'Pressure tensor by all types of the contour definition: IK, H, H-VR, IK-VR'
                        Else 
                              Write(*,*) 'READ_INPUT: INVALID VIRIAL CONTOUR TYPE FOR PRESSURE TENSOR CALCULATION'
                              STOP
                        End If

                        ! Set & Print contribution type
                        If (calcdef .eq. 'fw') Then
                              calc_type = 1
                              Write(*,*) 'Fluid-Wall contribution to the pressure tensor will be calculated...'
                        Else if (calcdef .eq. 'ff') Then
                              calc_type = 2
                              Write(*,*) 'Fluid-Fluid contribution to the pressure tensor will be calculated...'
                        Else if (calcdef .eq. 'fw&ff') Then
                              calc_type = 3
                              Write(*,*) 'Both Fluid-Fluid and Fluid-Wall contribution to the pressure tensor will be calculated...'
                        Else 
                              Write(*,*) 'READ_INPUT: INVALID calcdef input'
                              STOP
                        End If


                        ! Read in the bins for producing pressure tensor
                        Read(FILE_INPUT,*) virialpress_slit_bins

                        ! Read in calculation frequency
                        Read(FILE_INPUT,*) virialpress_slit_freq

                        ! Check the consistency of number of bins chosen
                        If (zden_bins .ne. virialpress_slit_bins) Then
                              Write(*,*) &
                              &'READ_INPUT: Number of bins should be consistent for zdensity and virial_press_slit'
                              STOP
                        End If

                        ! Print
                        Write(*,'(A,I0)') 'Calculation frequency:', virialpress_slit_freq


                  ! Calculate the pressure tensor for planar surface (perpendicular to z-direction) from virial route
                  Else If (Trim(optionname) .eq. 'virial_press') Then

                        If (.not. lzdensity) Then
                              Write(*,*) 'READ_INPUT: lzdensity should be turned on to calculate kinetic part of the pressure'
                              STOP
                        End If

                        lvirialpress = .true. 
                        lsampling = .true.

                        ! Read in the definition of integral contour for pressure tensor calculation
                        Read(FILE_INPUT,*) pressdef

                        ! Set & Print type
                        If ((pressdef .eq. 'IK') .or. (pressdef .eq. 'Irving-Kirkwood')) Then
                              virialpress_ctype = 1
                              Write(*,*) 'Pressure tensor calculation ONLY by Irving-Kirkwood definition'
                        Else if ((pressdef .eq. 'H') .or. (pressdef .eq. 'Harasima')) Then
                              virialpress_ctype = 2
                              Write(*,*) 'Pressure tensor calculation ONLY by Harasima definition'
                        Else if (pressdef .eq. 'IK&H') Then
                              virialpress_ctype = 3
                              Write(*,*) 'Pressure tensor calculation by both IK and Harasima definition'
                        Else 
                              Write(*,*) 'READ_INPUT: INVALID VIRIAL CONTOUR TYPE FOR PRESSURE TENSOR CALCULATION'
                              STOP
                        End If

                        ! Read in calculation frequency
                        Read(FILE_INPUT,*) virialpress_freq

                        ! Print
                        Write(*,'(A,I0)') 'Calculation frequency:', virialpress_freq


                  ! Calculate the pressure tensor for cylindrical geometry from virial route
                  Else If (Trim(optionname) .eq. 'virial_press_cylin') Then

                        If (.not. lrdensity) Then
                              Write(*,*) 'READ_INPUT: lrdensity should be turned on to calculate kinetic part of the pressure'
                              STOP
                        End If

                        lvirialpress_cylin = .true. 
                        lsampling = .true.

                        ! Read in the definition of integral contour for pressure tensor calculation
                        Read(FILE_INPUT,*) pressdef

                        ! Set & Print type
                        If ((pressdef .eq. 'IK') .or. (pressdef .eq. 'Irving-Kirkwood')) Then
                              virialpress_ctype = 1
                              Write(*,'(A)') &
                                    &'Pressure tensor calculation ONLY by Irving-Kirkwood definition HAS NOT BEEN INCLUDED YET'
                              STOP
                        Else if ((pressdef .eq. 'H') .or. (pressdef .eq. 'Harasima')) Then
                              virialpress_ctype = 2
                              Write(*,'(A)') 'Pressure tensor calculation ONLY by Harasima definition'
                        Else if (pressdef .eq. 'IK&H') Then
                              virialpress_ctype = 3
                              Write(*,'(A)') &
                                    &'Pressure tensor calculation by both IK and Harasima definition HAS NOT BEEN INCLUDED YET'
                              STOP
                        Else 
                              Write(*,'(A)') 'READ_INPUT: INVALID VIRIAL CONTOUR TYPE FOR PRESSURE TENSOR CALCULATION'
                              STOP
                        End If

                        ! Read in calculation frequency
                        Read(FILE_INPUT,*) virialpress_cylin_freq

                        ! Print
                        Write(*,'(A)') 'Pressure tensor in the cylindrical sytem will be calculated!'


                  ! 2D lattice constant 
                  Else if (Trim(optionname) .eq. 'lattice_const') Then

                        llattconst = .true.
                        lsampling = .true.

                        ! Read in number of confined layers to be considered
                        Read(FILE_INPUT,*) lattconst_n

                        ! Read in cutoff boundary for each layer
                        ! Lower boundary for each layer
                        Read(FILE_INPUT,*) (lattconst_cut(1,itype), itype = 1, lattconst_n)
                        ! Upper boundary for each layer
                        Read(FILE_INPUT,*) (lattconst_cut(2,itype), itype = 1, lattconst_n)

                        ! Read in calculation frequency
                        Read(FILE_INPUT,*) lattconst_freq

                        ! False check
                        If((steele_avgx(1) .eq. 0.0d0) .AND. (steele_avgx(2) .eq. 0.0d0))  Then
                              Write(*,*) "READ_INPUT: Averaing region needs to be specified for lattice_const"
                              STOP
                        End if

                        ! Print
                        Write(*,*) '2D lattice constant will be calculated...'


                  ! Isosteric heat of adsorption (added on Aug 8, 2018)
                  Else if (Trim(optionname) .eq. 'isosteric_heat') Then

                        lsampling = .true.
                        lqst = .true.

                        ! Make sure only single adsorbate component is present
                        If (initstyle(1,1) .eq. 'coords') Then
                              If (n_mol_types .gt. 2) Then
                                    Write(*,*) 'ERROR: isosteric_heat only works for single adsorbate component now!'
                                    STOP
                              End If
                        Else 
                              If (n_mol_types .gt. 1) Then
                                    Write(*,*) 'ERROR: isosteric_heat only works for single component now!'
                                    STOP
                              End If 
                        End If

                  ! Dump xyz coordinates of the whole system to the file (added on Nov 21, 2018)
                  Else if (Trim(optionname) .eq. 'dump_xyz_all') Then

                        lsampling = .true.
                        ldumpxyz = .true.

                        ! Read in output frequency 
                        Read(FILE_INPUT,*) dumpxyz_freq

                  ! Dump xyz coordinates of the specified molecule type to the file (added on Jan 22, 2019)
                  Else if (Trim(optionname) .eq. 'dump_xyz_mol') Then

                        lsampling = .true.
                        ldumpxyz = .true.
                        
                        ! Read in dumping modes
                        ! 'center' - center-of-mass only; 'sites' - all sites 
                        Read(FILE_INPUT,*) dumpxyz_mode

                        ! Read in output frequency and moltype 
                        Read(FILE_INPUT,*) dumpxyz_freq, dump_moltype


                  ! Dump density of the whole system to the file (added on April 25, 2019)
                  Else if (Trim(optionname) .eq. 'dump_density') Then

                        lsampling = .true.
                        ldumpdensity = .true.
                        
                        ! Read in output frequency 
                        Read(FILE_INPUT,*) dumpdensity_freq


                  ! Option to turn off/on periodic boundary conditions/hard-wall boundary condition
                  ! Added on July 14, 2019
                  Else If (Trim(optionname) .eq. 'no_pbc') Then

                        ! This option not working for Ewald method
                        If (lewld) Then
                              Write(*,*) 'READ_INPUT: no_pbc is not compatiable with Ewald'
                              STOP
                        End If
                        ! not working for tail correction which assumes g(r) =1 in long range
                        If (ltailc) Then
                              Write(*,*) 'READ_INPUT: no_pbc is not compatiable with tail correction'
                              STOP
                        End If

                        lno_pbc = .true. 

                        ! Print
                        Write(*,'(A)') 'WARNING: Hard wall boundary will be applied (instead of PBC)!'  


                  ! Using cell list and linked-list algorithm to accelerate simulation
                  ! Added on July 15, 2019
                  Else If (Trim(optionname) .eq. 'cell_list') Then

                        ! check compatibility
                        If (lewld) Then
                              Write(*,*) 'READ_INPUT: Cell list is not compatible with current Ewald method (ewald_fix_kmax)'
                              STOP
                        End If

                        lclist = .true. 

                        ! Allocate cell list variables
                        Call global_allocate('Clist')


                  ! Dump total energy of the whole system to the file (added on October 13, 2019)
                  Else if (Trim(optionname) .eq. 'dump_energy') Then

                        lsampling = .true.
                        ldumpenergy = .true.
                        
                        ! Read in output frequency 
                        Read(FILE_INPUT,*) dumpenergy_freq

                  ! Write instantaneous configuraitons (in XYZ format) for restart (added on October 13, 2019)
                  Else if (Trim(optionname) .eq. 'write_restart') Then

                        lsampling = .true.
                        lwriterestart = .true.
                        
                        ! Read in output frequency 
                        Read(FILE_INPUT,*) rst_freq

                  ! Dump instant virial, hypervirial and hydrostatic pressure to file (added on Jan, 2020)
                  Else if (Trim(optionname) .eq. 'dump_virial') Then

                        ! check compatibility
                        If (lewld) Then
                              Write(*,*) 'READ_INPUT: Virial computation is ONLY compatible with pairwise forces for now!'
                              Write(*,*) 'READ_INPUT: Virial computation is NOT compatible with the Ewald method. Revise the code!'
                              STOP
                        End If

                        If (n_sites_max .GT. 1) Then
                              Write(*,*) 'READ_INPUT: Virial computation is ONLY valid for atomic fluids now!'
                              STOP
                        End If

                        lsampling = .true.
                        ldumpvir = .true.
                        
                        ! Read in output frequency 
                        Read(FILE_INPUT,*) dumpvir_freq

                  ! Recalculate and check total energy to mitigate the error in accumulated energies
                  ! Added on June 9, 2020, following the "check_energy" subroutine in the DL_MOTE-2 program
                  ! Write instantaneous configuraitons (in XYZ format) for restart (added on October 13, 2019)
                  Else if (Trim(optionname) .eq. 'check_energy') Then
                       
                        ! Read in check frequency and overwrite the default value
                        Read(FILE_INPUT,*) check_freq

      		! EOF
      		Else if((Trim(optionname) .eq. 'EOF') .or. (Trim(optionname) .eq. 'eof')) Then
      			! exit this loop and continue
      			EXIT

      		! Throw and error
      		Else
      			Write(*,*) 'READ_INPUT: INVALID SAMPLING OPTIONS COMMAND'
      			STOP
      		End If

      	! End reading optionname
      	End do
      	
      Else
      	Write(*,*) 'READ_INPUT: FAILED TO READ IN SAMPLING OPTIONS LABEL'
      	STOP
      End if

      ! Allocate sampling variables
      if(lsampling) Call global_allocate('Sampling')


      ! Initialize system parameters
      Call initialize('System')



      Close(FILE_INPUT)

      Return 

      End Subroutine read_input


















      