! ==========================================================
! This subroutine is used to initialize the system parameters
! Created on Dec. 10th, 2016 by Kaihang Shi
! Last modified
! ==========================================================

      Subroutine initialize(selector)

      Use global

      IMPLICIT NONE

      ! Passed
      Character(Len=*) :: selector

      ! Local
      Integer :: ierr, icount
      Integer :: isite, jsite, itype, jtype, ibox, isitetype, inonbond
      Double Precision :: prob
      


      Select Case (selector)

        ! Initialize coarse-grained wall potential with accessible volume/length read from full structure
        Case('CG_WALL_STRC')

            ! initialize i/o variables
            ierr = 0
            ! Initialize point counter (defined in global.f90)
            cg_wall_maxnum = 0
            cg_strc_maxnum = 0

            ! Calculate & Print wall position
            Write(*,'(A,F7.3)') 'Position for the Bottom CG Wall: ', cg_wall_position(1)
            cg_wall_position(2) = box(3,1) - cg_wall_position(1)
            Write(*,'(A,F7.3)') 'Position for the Top CG Wall: ', cg_wall_position(2)

            ! Open pre-calculated potential file
            !!!!!!!!!!! Check out FILE NAME
            Open(Unit=FILE_CG_WALL,File='cg_wall.in',Status='Old',Access='Sequential',Action= 'Read')

            ! Keep reading data
            Do
                cg_wall_maxnum = cg_wall_maxnum + 1
                Read(FILE_CG_WALL,*,IOSTAT=ierr) cg_wall_dgrid(cg_wall_maxnum), cg_wall_egrid(cg_wall_maxnum)
                If (ierr .eq. 0) Then
                    Continue
                Else
                    cg_wall_maxnum = cg_wall_maxnum - 1
                    EXIT
                End If
            End Do

            Write(*,*) "Successfully read data from cg_wall.in"
            Close(FILE_CG_WALL)

            !!!!! Code block for modified version of cg_wall_strc !!!!!
            ! ! Open pre-calculated potential file
            ! !!!!!!!!!!! Check out FILE NAME
            ! cg_wall_maxnum = 0
            ! Open(Unit=FILE_CG_WALL1,File='cg_wall1.in',Status='Old',Access='Sequential',Action= 'Read')

            ! ! Keep reading data
            ! Do
            !   cg_wall_maxnum = cg_wall_maxnum + 1
            !   Read(FILE_CG_WALL1,*,IOSTAT=ierr) cg_wall1_dgrid(cg_wall_maxnum), cg_wall1_egrid(cg_wall_maxnum)
            !   If (ierr .eq. 0) Then
            !       Continue
            !   Else
            !       cg_wall_maxnum = cg_wall_maxnum - 1
            !       EXIT
            !   End If
            ! End Do

            ! Write(*,*) "Successfully read data from cg_wall.in"
            ! Close(FILE_CG_WALL1)
            !!!!!!!!!!!!!!!!!!!!!!!!!

            ! Do icount = 1, cg_wall_maxnum
            !   Write(*,*) cg_wall_dgrid(icount), cg_wall_egrid(icount)
            ! End Do
            ! STOP

            ! Initialize I/O variable
            ierr = 0

            ! Open pre-calculated accessible volume/length file from full strcture
            Open(Unit=FILE_CG_WALL_STRC,File='cg_wall_strc.in',Status='Old',Access='Sequential',Action= 'Read')

            ! Keep reading data
            Do
                cg_strc_maxnum = cg_strc_maxnum + 1
                Read(FILE_CG_WALL_STRC,*,IOSTAT=ierr) cg_strc_dgrid(cg_strc_maxnum), cg_strc_dacc(cg_strc_maxnum)
                If (ierr .eq. 0) Then
                    Continue
                Else
                    cg_strc_maxnum = cg_strc_maxnum - 1
                    EXIT
                End If
            End Do
            Write(*,*) "Successfully read data from cg_wall_strc.in"
            Close(FILE_CG_WALL_STRC)

        ! Initialize coarse-grained wall potential with piecewise fluid-fluid interactions
        Case('CG_WALL_FFPW')

            ! initialize i/o variables
            ierr = 0
            ! Initialize point counter (defined in global.f90)
            cg_wall_maxnum = 0
            cg_ff_maxnum = 0

            ! Calculate & Print wall position
            Write(*,'(A,F7.3)') 'Position for the Bottom CG Wall: ', cg_wall_position(1)
            cg_wall_position(2) = box(3,1) - cg_wall_position(1)
            Write(*,'(A,F7.3)') 'Position for the Top CG Wall: ', cg_wall_position(2)

            ! Open pre-calculated potential file
            Open(Unit=FILE_CG_WALL,File='cg_wall.in',Status='Old',Access='Sequential',Action= 'Read')

            ! Keep reading data
            Do
                cg_wall_maxnum = cg_wall_maxnum + 1
                Read(FILE_CG_WALL,*,IOSTAT=ierr) cg_wall_dgrid(cg_wall_maxnum), cg_wall_egrid(cg_wall_maxnum)
                If (ierr .eq. 0) Then
                    Continue
                Else
                    cg_wall_maxnum = cg_wall_maxnum - 1
                    EXIT
                End If
            End Do

            Write(*,*) "Successfully read data from cg_wall.in"
            Close(FILE_CG_WALL)

            ! Do icount = 1, cg_wall_maxnum
            !   Write(*,*) cg_wall_dgrid(icount), cg_wall_egrid(icount)
            ! End Do
            ! STOP

            ! Initialize I/O variable
            ierr = 0

            ! Open pre-calculated effective fluid-fluid potential file
            Open(Unit=FILE_EFF_SIG,File='eff_sig.in',Status='Old',Access='Sequential',Action= 'Read')

            ! Keep reading data
            Do
                cg_ff_maxnum = cg_ff_maxnum + 1
                Read(FILE_EFF_SIG,*,IOSTAT=ierr) eff_sig_z(cg_ff_maxnum), eff_sig(cg_ff_maxnum)
                If (ierr .eq. 0) Then
                    Continue
                Else
                    cg_ff_maxnum = cg_ff_maxnum - 1
                    EXIT
                End If
            End Do
            Write(*,*) "Successfully read data from eff_sig.in"
            Close(FILE_EFF_SIG)

        ! Initialize coarse-grained wall potential
        Case('CG_WALL')

            ! initialize i/o variables
            ierr = 0
            ! Initialize point counter (defined in global.f90)
            cg_wall_maxnum = 0

            ! Calculate & Print wall position
            Write(*,'(A,F7.3)') 'Position for the Bottom CG Wall: ', cg_wall_position(1)
            cg_wall_position(2) = box(3,1) - cg_wall_position(1)
            Write(*,'(A,F7.3)') 'Position for the Top CG Wall: ', cg_wall_position(2)

            ! Open pre-calculated potential file
            Open(Unit=FILE_CG_WALL,File='cg_wall.in',Status='Old',Access='Sequential',Action= 'Read')

            ! Keep reading data
            Do
                cg_wall_maxnum = cg_wall_maxnum + 1
                Read(FILE_CG_WALL,*,IOSTAT=ierr) cg_wall_dgrid(cg_wall_maxnum), cg_wall_egrid(cg_wall_maxnum)
                If (ierr .eq. 0) Then
                    Continue
                Else
                    cg_wall_maxnum = cg_wall_maxnum - 1
                    EXIT
                End If
            End Do

            ! Do icount = 1, cg_wall_maxnum
            !   Write(*,*) cg_wall_dgrid(icount), cg_wall_egrid(icount)
            ! End Do
            ! STOP

            ! Initialize external field
            If (field_type .eq. CG_WALL_COS) Then
                ! Calculate periodicity for the cosine function [A]
                cg_wall_cospe = box(1,1)/cg_wall_cosnp
                ! Print
                Write(*,'(A,F8.4)') 'Periodicity of the Cosine function for cg_wall_cos is: ', cg_wall_cospe
                Write(*,'(A,F8.4)') 'Amplitude of the Cosine function for cg_wall_cos is: ', cg_wall_cosap
            End If


        ! Initialize Steele potential interaction parameters
        Case('Steele')

            ! False check
            If ((steele_position(1) .gt. box(3,1)) .or. (steele_position(1) .lt. 0.0d0)) Then
                Write(*,*) 'INITIALIZE.F90: INVALID steele_position. RESET IT!'
                STOP
            End If

            ! Print wall position
            Write(*,'(A,F7.3)') 'Position for Downside Wall: ', steele_position(1)

            If ((field_type .eq. STEELE_SLIT_PORE) .OR. (field_type .eq. STEELE_SLIT_FINITEX)) Then
                ! Upper wall position
                steele_position(2) = box(3,1) - steele_position(1)
                Write(*,'(A,F7.3)') 'Position for Upper Wall: ', steele_position(2)
            End If
            

            Write(*,*) 'Initialize 10-4-3 Steele potential interaction parameters:'

            ! Loop over molecule types
            Do itype = 1, n_mol_types

                ! Loop over all the sites
                Do isitetype = 1, n_site_types(itype)

                    ! Process according to site name
                    Do inonbond = 1, n_nonbond_types
                        If (steele_site_name(inonbond) .eq. site_type_name(isitetype,itype)) Then
                            ! Reset array for later more convenient computational use in the program
                            steele_sigmasf(isitetype,itype) = steele_sigmasf(inonbond,n_mol_types_max)
                            steele_sigmasfsq(isitetype,itype) = steele_sigmasf(isitetype,itype)**2
                            steele_epsilonsf(isitetype,itype) = steele_epsilonsf(inonbond,n_mol_types_max)

                            ! Print to screen
                            Write(*,'(A,F15.3,2X,F15.3)') site_type_name(isitetype,itype), &
                                          & steele_sigmasf(isitetype,itype), &
                                          & steele_epsilonsf(isitetype,itype)

                            EXIT
                        End If
                    End Do

                    
                End Do
                
            End Do

            ! Initialize parameters for steele finite slit pore model
            If (field_type .eq. STEELE_SLIT_FINITEX) Then

                ! Set the pore length in x-direction which is now automatically set to 1/3 of the Lx
                ! 1 - lower boundary, 2- upper boundary
                steele_posx(1) = 17.0d0
                steele_posx(2) = 97.0d0

                ! Set the averaging region boundary which is now automatically set to 6sig_ff
                ! False check, 
                ! If (steele_posx(1) .lt. 68.866d0) Then
                !   Write(*,*) 'INITIALIZE: Lx is too small for steele_slit_finitex potential'
                !   STOP
                ! End If

                ! away from pore mouth of at least 5sigma_ff. Now set to 30 A (larger than Yun's paper).
                steele_avgx(1) = steele_posx(1) + 29.785d0
                steele_avgx(2) = steele_avgx(1) + 20.43d0

                ! Write to file
                Write(*,'(A,F10.4)') 'steele_posx(1) = ', steele_posx(1)
                Write(*,'(A,F10.4)') 'steele_posx(2) = ', steele_posx(2)
                Write(*,'(A,F10.4)') 'steele_avgx(1) = ', steele_avgx(1)
                Write(*,'(A,F10.4)') 'steele_avgx(2) = ', steele_avgx(2)


                
            End If
        
        ! Initialize potential parameters
        Case('Potential')

            ! Output potential form to screen
            If (potential .eq. LENNARD_JONES) Then
                Write(*,*) ' '
                Write(*,*) 'Potential form: Lennard-Jones'
            Else
                Write(*,*) 'INITIALIZE: INVALID potential form'
                STOP
            End If 
                
            ! Calculate interaction parameters using different mixing rules set in input.in
            If (mix_rule .eq. LORENTZ_BERTHELOT) Then

                ! Write header
                Write(*,*) 'Mixing rules: Lorentz-Berthelot'
                Write(*,*) 'No.  Atom(i)   No.  Atom(j)         sigma[A]           epsilon[K]         charge[e]'

                ! Calculate interaction parameters 
                Do itype = 1, n_mol_types
                    Do jtype = 1, n_mol_types
                        Do isite = 1, n_site_types(itype)
                            Do jsite = 1, n_site_types(jtype)

                                epsilon(isite,itype,jsite,jtype) = &
                                &Dsqrt(epsilon(isite,itype,isite,itype)*epsilon(jsite,jtype,jsite,jtype))
                                sigma(isite,itype,jsite,jtype) = &
                                &0.5d0 *(sigma(isite,itype,isite,itype)+sigma(jsite,jtype,jsite,jtype))
                                sigmasq(isite,itype,jsite,jtype) = sigma(isite,itype,jsite,jtype)**2    

                            End do 
                        End Do 
                    End do  
                End Do  

                ! Calculate point charges [e]
                Do itype = 1, n_mol_types
                    Do jtype = 1, n_mol_types
                        Do isite = 1, n_site_types(itype)
                            Do jsite = 1, n_site_types(jtype)

                                qsq(isite,itype,jsite,jtype) = &
                                &q(isite,itype,isite,itype)*q(jsite,jtype,jsite,jtype)

                            End do 
                        End Do 
                    End do  
                End Do  



                ! Print vdW potential parameters to the screen
                ! Print intramolecular cross term
                Do itype = 1, n_mol_types 
                    Do isite = 1, n_site_types(itype)
                        Do jsite = isite, n_site_types(itype)

                            If (site_type_name(isite,itype) .eq. site_type_name(jsite,itype)) Then
                                ! Write sigma, epsilon and charge value of each pair to screen
                                Write(*,'(I2,3X,A6,4X,I2,3X,A6,10X,F9.4,10X,F9.4,10X,F9.4)') &
                                &isite, site_type_name(isite,itype), jsite, site_type_name(jsite,itype),&
                                &sigma(isite,itype,jsite,itype), epsilon(isite,itype,jsite,itype), &
                                &q(isite,itype,jsite,itype)
                            Else
                                ! Only Write sigma and epsilon value of each pair to screen
                                Write(*,'(I2,3X,A6,4X,I2,3X,A6,10X,F9.4,10X,F9.4)') &
                                &isite, site_type_name(isite,itype), jsite, site_type_name(jsite,itype),&
                                &sigma(isite,itype,jsite,itype), epsilon(isite,itype,jsite,itype)
                            End If  
                            
                        End Do
                    End Do                  
                End Do

                ! Print intermolecular cross term
                Do itype = 1, n_mol_types - 1
                    Do jtype = itype + 1, n_mol_types
                        Do isite = 1, n_site_types(itype)
                            Do jsite = 1, n_site_types(jtype)

                                ! Only Write sigma and epsilon value of each pair to screen
                                Write(*,'(I2,3X,A6,4X,I2,3X,A6,10X,F9.4,10X,F9.4)') &
                                &isite, site_type_name(isite,itype), jsite, site_type_name(jsite,jtype),&
                                &sigma(isite,itype,jsite,jtype), epsilon(isite,itype,jsite,jtype)
                                
                            End Do
                        End Do                      
                    End Do                  
                End Do





            Else if (mix_rule .eq. EXPLICIT) Then


                Write(*,*) 'No.  Atom(i)   No.  Atom(j)         sigma[A]           epsilon[K]         charge[e]'

                ! Calculate interaction parameters 
                Do itype = 1, n_mol_types
                    Do isite = 1, n_site_types(itype)

                        ! Initialize 
                        inonbond = 1 

                        Do jtype = 1, n_mol_types
                            Do jsite = 1, n_site_types(jtype)

                                epsilon(isite,itype,jsite,jtype) = epsilon(isite,itype,inonbond,inonbond)
                                sigma(isite,itype,jsite,jtype) = sigma(isite,itype,inonbond,inonbond)
                                sigmasq(isite,itype,jsite,jtype) = sigma(isite,itype,jsite,jtype)**2    

                                inonbond = inonbond + 1

                            End do 
                        End Do 

                    End do  
                End Do  


                ! Calculate point charges [e]
                Do itype = 1, n_mol_types
                    Do jtype = 1, n_mol_types
                        Do isite = 1, n_site_types(itype)
                            Do jsite = 1, n_site_types(jtype)

                                ! Original charge 
                                qsq(isite,itype,jsite,jtype) = &
                                &q(isite,itype,isite,itype)*q(jsite,jtype,jsite,jtype)

                                ! Perturbated charge (for Ad-Surf)
                                qsqex(AS,isite,itype,jsite,jtype) = &
                                &qex(AS,isite,itype,isite,itype)*qex(AS,jsite,jtype,jsite,jtype)

                            End do 
                        End Do 
                    End do  
                End Do  



                ! Print vdW potential parameters to the screen
                ! Print intramolecular cross term
                Do itype = 1, n_mol_types 
                    Do isite = 1, n_site_types(itype)
                        Do jsite = isite, n_site_types(itype)

                            If (site_type_name(isite,itype) .eq. site_type_name(jsite,itype)) Then
                                ! Write sigma, epsilon and charge value of each pair to screen
                                Write(*,'(I2,3X,A6,4X,I2,3X,A6,10X,F9.4,10X,F9.4,10X,F9.4,10X,F9.4)') &
                                &isite, site_type_name(isite,itype), jsite, site_type_name(jsite,itype),&
                                &sigma(isite,itype,jsite,itype), epsilon(isite,itype,jsite,itype), &
                                &q(isite,itype,jsite,itype), qex(AS,isite,itype,jsite,itype)
                            Else
                                ! Only Write sigma and epsilon value of each pair to screen
                                Write(*,'(I2,3X,A6,4X,I2,3X,A6,10X,F9.4,10X,F9.4)') &
                                &isite, site_type_name(isite,itype), jsite, site_type_name(jsite,itype),&
                                &sigma(isite,itype,jsite,itype), epsilon(isite,itype,jsite,itype)
                            End If  
                            
                        End Do
                    End Do                  
                End Do

                ! Print intermolecular cross term
                Do itype = 1, n_mol_types - 1
                    Do jtype = itype + 1, n_mol_types
                        Do isite = 1, n_site_types(itype)
                            Do jsite = 1, n_site_types(jtype)

                                ! Only Write sigma and epsilon value of each pair to screen
                                Write(*,'(I2,3X,A6,4X,I2,3X,A6,10X,F9.4,10X,F9.4)') &
                                &isite, site_type_name(isite,itype), jsite, site_type_name(jsite,jtype),&
                                &sigma(isite,itype,jsite,jtype), epsilon(isite,itype,jsite,jtype)
                                
                            End Do
                        End Do                      
                    End Do                  
                End Do
                

            End If

        ! Initialize system parameters
        Case('System')


            ! Initialize thermodynamic statistics
            energy = 0.0d0
            energy_disp = 0.0d0
            energy_tail = 0.0d0
            energy_field = 0.0d0
            ewld_tot = 0.0d0
            ewld_real = 0.0d0
            ewld_fourier = 0.0d0
            ewld_slab = 0.0d0

            
            ! Initialize average block statistics
            avg_nmol = 0.0d0
            avg_rho = 0.0d0
            avg_eng = 0.0d0
            avg_eng_disp = 0.0d0
            avg_vol = 0.0d0

            ! Initialize MC move statistics 
            in_stat = 0
            rem_stat = 0

            ! Initialize widom insertion variables
            widom_sample = 0.0d0

            ! Initialize sampling variables
            If (lsampling) Then

                ! z_density
                If (lzdensity) Then

                    ! Initialize statistics
                    avg_zden = 0.0d0

                    ! If graphene slit pore, bins are set for accessible volume
                    IF (lfield .and. (field_type .eq. STEELE_SLIT_PORE)) THEN
                        ! Initialize related parameters
                        dz = (box(3,1)-2*steele_position(1))/zden_bins
                        dvol = box(1,1)*box(2,1)*dz
                    Else if (lfield .and. (field_type .eq. STEELE_SLIT_FINITEX)) Then
                        dz = (box(3,1)-2*steele_position(1))/zden_bins
                        dvol = (steele_avgx(2)-steele_avgx(1))*box(2,1)*dz
                    ! Modified on April 29, 2020
                    ! Now we want to sample the whole box instead of the pore volume
                    Else if (lfield .and. (field_type .eq. HARD_SLIT_FINITEX)) Then
                        ! Assume three graphene layers with interlayer spacing 3.35A
                        dz = (box(3,1)-2*wall_radius)/zden_bins
                        dvol = (steele_avgx(2)-steele_avgx(1))*box(2,1)*dz
                    End if

                End If


                ! r_density
                If (lrdensity) Then

                    ! Initialize statistics
                    avg_rden = 0.0d0

                    ! If graphene slit pore, bins are set for accessible volume
                    IF (lfield) Then
                        Write(*,*) 'Double check r_density setting with the current external field_type'
                        STOP

                    Else 
                        ! Initialize related parameters
                        ! Unlike in z-density case, we set up a cutoff radius for r-density calculation here.
                        rden_dr = rden_cut/DBLE(rden_bins)
                        rden_drsq = rden_dr**2
                        rden_lim = rden_cut+r_cut
                        ! delta value used in cylindrical pressure tensor of Harasima contour
                        delrr = rden_dr/2.0d0
                    End if

                End If


                ! surface_excess
                If (lsurfex) Then
                    
                    ! Initialize statistics
                    avg_surfex = 0.0d0

                    ! Loop over molecule types
                    Do itype = 1, n_mol_types

                        ! Skip external structure
                        If (initstyle(itype,1) .eq. 'coords') CYCLE

                        ! Calculate number of molecules in accessible volume
                        surfex_n_mol(itype) = surfex_vol * surfex_bulk(itype)

                    End Do

                End If

                ! Pressure tensor for slit pore from thermodynamic route
                If (lthermopress_slit) Then

                    ! Initialize statistics
                    thermopress_slit_sample = 0.0d0
                    thermopress_slit_pt = 0.0d0
                    thermopress_slit_pn = 0.0d0

                    ! Initialize related parameters for bins
                    thermopress_slit_dz = box(3,1)/thermopress_slit_bins
                    thermopress_slit_dvol = box(1,1)*box(2,1)*thermopress_slit_dz
                    
                End If

                ! Pressure tensor for slit pore from virial route
                If (lvirialpress_slit) Then

                    ! Initialize statistics
                    virialpress_slit_pt = 0.0d0
                    virialpress_slit_pn = 0.0d0

                    ! If graphene slit pore, bins are set for accessible volume
                    IF (lfield .and. (field_type .eq. STEELE_SLIT_PORE)) THEN
                        ! Initialize related parameters for bins
                        virialpress_slit_dz = (box(3,1)-2*steele_position(1))/virialpress_slit_bins

                    Else if (lfield .and. (field_type .eq. STEELE_SLIT_FINITEX)) Then
                        virialpress_slit_dz = (box(3,1)-2*steele_position(1))/virialpress_slit_bins

                    Else if (lfield .and. (field_type .eq. HARD_SLIT_FINITEX)) Then
                        ! Assume three graphene layers with interlayer spacing 3.35A
                        ! Revised on 4/29/2020 to be consistent with zdensity
                        virialpress_slit_dz = (box(3,1)-2*wall_radius)/virialpress_slit_bins    
                                            
                    Else 
                        ! Initialize related parameters for bins
                        virialpress_slit_dz = box(3,1)/virialpress_slit_bins
                    End if
                    
                End If

                ! Pressure tensor for planar surface from virial route
                If (lvirialpress) Then

                    ! Initialize statistics
                    virialpress_pt = 0.0d0
                    virialpress_pn = 0.0d0
                    
                End If

                ! Pressure tensor for cylindrical system from virial route
                If (lvirialpress_cylin) Then

                    ! Initialize statistics
                    virialpress_cylin_ptt = 0.0d0
                    virialpress_cylin_ptz = 0.0d0
                    virialpress_cylin_pnr = 0.0d0
                    
                End If

                ! Lattice constat
                If (llattconst) Then
                    ! Initialize
                    avg_lattconst = 0.0d0
                End If

                ! Isosteric heat of adsorption
                If (lqst) Then
                    ! Initialize
                    avg_qst = 0.0d0
                    qst = 0.0d0
                End If

                ! Initialize virial
                If (ldumpvir) Then
                    vir = 0.0d0
                    vir_hyp = 0.0d0
                End If

            End if


            ! Initialize cell list
            ! Added on July 15, 2019
            If (lclist) Then
                ! Construct the cell and assign neighbors
                Call set_clist

                ! Print out initial cell list information
                Write(*,'(A)') 'Cell list will be applied to accelerate simulations!'  
                Write(*,'(A,I3,3X,I3,3X,I3)') 'Initial number of cells in x, y, z-direction: ', clist_nx, clist_ny, clist_nz
                Write(*,'(A,F5.2,3X,F5.2,3X,F5.2)') 'Initial size of cells in x, y, z-direction: ', clist_dx, clist_dy, clist_dz
                
                ! Construct a new cell list
                Call update_clist
            End If
            

            ! Initialize maximum translational move distance [A]
            max_trans = 1.0d0

            ! Initialize the maximum rotation
            max_rotat = 5.0d0

            ! Initalize the maximum volume change
            max_vol = 0.01 


            ! Calculate MC run parameters
            ! Total number of simulation steps
            n_steps_tot = n_blocks_tot * block_size
            ! Number of production blocks
            n_blocks_prod = n_blocks_tot - n_blocks_equil

            ! Loop over each molecular type
            Do itype = 1, n_mol_types

                ! Successively change chempot 10 times (add 8-13-2017 for update_chempot subroutine)
!               dmu(itype) = (mu(itype) + 2084.716617d0)/10.0d0
                
                ! Calculate de Broglie wavelength per molecule [Angstrom] 
                lambda(itype) = (h/Dsqrt(2.0d0*Pi*temp*Kb*mol_mass(itype)/(Na*1000.0d0)))*1.0d10

                ! Calculate activity constant * volume (GCMC)
                if(ensmbl .eq. ENS_uVT) mol_act(itype) = vol_pore(1)*dEXP(mu(itype)/temp)/(lambda(itype)**3)
                
            End Do

            ! Calculate ideal part of chemical potential [K]
            If (lwdm) Then
                
                ! Loop over each molecule type
                Do itype = 1, n_mol_types

                    ! Choose different ensemble
                    ! Double check the second term when simulation with binary/ternary system...
                    If (ensmbl .eq. ENS_NPT) muid(itype) = -temp*dLOG(temp/(PVTOK*press*lambda(itype)**3))&
                                                         & -temp*dLOG(DBLE(n_mol_tot(1)+1)/DBLE(n_mol(itype,1)+1))

                End Do
                
            End If


            ! Print pre-set MC move probabilities to the screen
            ! Probabilities of selecting MC move types
            Write(*,'(A,F5.3)') 'Probabilities of selecting translational move: ',&
            &move_type_prob(TRANSLATION)
            Write(*,'(A,F5.3)') 'Probabilities of selecting rotational move: ', &
            &move_type_prob(ROTATION) - move_type_prob(TRANSLATION)
            Write(*,'(A,F5.3)') 'Probabilities of selecting transfer move: ',&
            &move_type_prob(TRANSFER) - move_type_prob(ROTATION)
            Write(*,'(A,F5.3)') 'Probabilities of selecting volume move: ',&
            &move_type_prob(VOLCHANGE) - move_type_prob(TRANSFER)





            
        Case default
            Write(*,*) 'FATAL ERROR: INVALID SELECTOR IN INITIALIZE SUBROUTINE'
            STOP
            
      End Select



      Return

      End Subroutine initialize


















