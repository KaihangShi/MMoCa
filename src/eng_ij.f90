! ==========================================================
! This subroutine is used to calculate dispersion energy, virial 
! & Real space Ewald sum between imol and jmol 
! Created on 12-26-2016 by Kaihang Shi
!
! Modified on 1-26-2017 notes: 
! We directly use site-site distance rather than
! COM-COM distance to compare with r_cut, this move might be more
! expansive, but more accurate. (results consistent with NIST Water Ref)
!
! Note: ewld_real is truncated at rcelect, different from vdW cut off
!       r_cut. rcelect is changing along with the volume for 'ewald_fix_kmax' mode, 
!       r_cut is a constant throughout the simulation.
! ==========================================================

      Subroutine eng_ij(ibox,imol,jmol,engdisp,engelec,vir_update,virhyp_update)
      

      use global

      IMPLICIT NONE

      ! Passed 
      Integer :: ibox, imol, jmol
      Double Precision :: engdisp, engelec, vir_update, virhyp_update

      ! Local
      Integer :: itype, jtype, isite, jsite, isitetype, jsitetype
      Double Precision :: rxijc, ryijc, rzijc, rxijs, ryijs, rzijs
      Double Precision :: rxis, ryis, rzis
      Double Precision :: rijssq, rijcsq, rijs
      Double Precision :: lj_factor, rdudr
      Double Precision :: engcgff
!        Logical :: WITHIN




      ! Initialize overlap flag
      ! Assume no overlap at the beginning
      OVERLAP = .False.
!        WITHIN = .FALSE.

      ! Initialize energy contribution
      engdisp = 0.0d0
      engelec = 0.0d0
      vir_update = 0.0d0
      virhyp_update = 0.0d0

      ! Get the molecule types 
      itype = mol_type(imol,ibox)
      jtype = mol_type(jmol,ibox)

      ! Check initial style
      ! If initstyle(itype,ibox) = 'coords', we assume this molecule type is external structure:
      ! open surface or pore structure. No calculation of COM separation is needed
      ! If 'coords' is specified, we assume imol = 1 corresponds to the ONLY external structure molecule 
      If (initstyle(itype,ibox) .eq. 'coords') Then

        ! Directly loop over imol's sites
        Do isite = 1, n_sites(itype)

            ! Get isite type 
            isitetype = site_type(isite,itype)

            rxis = rx_s(isite,imol,ibox)
            ryis = ry_s(isite,imol,ibox)
            rzis = rz_s(isite,imol,ibox)

            ! loop over jmol's sites (adsorbates)
            Do jsite = 1, n_sites(jtype)

                ! Get jsite type 
                jsitetype = site_type(jsite,jtype)

                ! Calculate the site separation vector 
                rxijs = rx_s(jsite,jmol,ibox) - rxis
                ryijs = ry_s(jsite,jmol,ibox) - ryis
                rzijs = rz_s(jsite,jmol,ibox) - rzis

                ! 
                If (lno_pbc) Then
                    Continue
                Else
                    ! Apply minimum image convension
                    rxijs = rxijs - dNINT(rxijs/box(1,ibox))*box(1,ibox)
                    ryijs = ryijs - dNINT(ryijs/box(2,ibox))*box(2,ibox)
                    rzijs = rzijs - dNINT(rzijs/box(3,ibox))*box(3,ibox)
                      
                End If


                ! Calculate the square of the site scalar separation distance
                rijssq = rxijs**2 + ryijs**2 + rzijs**2

                ! Check with rcelect first, because rcelect >= r_cut always
                ! Calculate Real space Ewald sum contribution
                If (lewld) Then

                    ! Newly added spherical cut-off method to calculate Coulomb interactions
                    ! between adsorbates and surface sites
                    If (lcoulsc) Then

                        ! Check with rscelect
                        if(rijssq .gt. rscelectsq) GOTO 100

                        If (qsq(isitetype,itype,jsitetype,jtype) .ne. 0.0d0) Then
                            rijs = dSQRT(rijssq)
                            engelec = engelec + qsq(isitetype,itype,jsitetype,jtype)/rijs
                        End If

                        GOTO 100

                    End If

                    ! Double check rcelect
                    If (rcelect .ne. MINVAL(box(:,ibox)/2.0d0)) Then
                        Write(*,*) 'ENG_IJ: RCELECT NOT EQUAL TO HALF OF MIN BOX LENGTH'
                        STOP
                    End If

                    ! Check with rcelect
                    if(rijssq .gt. rcelectsq) CYCLE

                    ! Mixing rule
                    If (mix_rule .eq. EXPLICIT) Then

                        ! If don't want perturbate charge, just set 'qex' to original charge
                        If (qsqex(AS,isitetype,itype,jsitetype,jtype) .ne. 0.0d0) Then
                              rijs = dSQRT(rijssq)
                              engelec = engelec + &
                                    & qsqex(AS,isitetype,itype,jsitetype,jtype)*dERFC(alpha*rijs)/rijs
                        End If

                    ! Redundant at this point
                    Else if (mix_rule .eq. LORENTZ_BERTHELOT) Then
                        If (qsq(isitetype,itype,jsitetype,jtype) .ne. 0.0d0) Then
                              rijs = dSQRT(rijssq)
                              engelec = engelec + &
                                    & qsq(isitetype,itype,jsitetype,jtype)*dERFC(alpha*rijs)/rijs
                        End If
                    End If

                ! End Coulombic calculation 
                End If

                ! Continue here when direct Coulombic calculation finishes
100                   Continue

                ! Check with r_cutoff
                If (rijssq .gt. r_cutsq) CYCLE

                ! Check with r_min
                If (rijssq .lt. r_minsq) Then
                    OVERLAP = .TRUE.
                    Return
                End If  

                ! Calculate dispersion energy
                If ((potential .eq. LENNARD_JONES) .and. (epsilon(isitetype,itype,jsitetype,jtype) .ne. 0.0d0)) Then

                    ! Calculate lJ factor
                    lj_factor = (sigmasq(isitetype,itype,jsitetype,jtype)/rijssq)**3
                    
                    ! Calculate dispersion energy using 12-6 Lennard-Jones potential
                    engdisp = engdisp + &
                    & 4.0d0*epsilon(isitetype,itype,jsitetype,jtype)*lj_factor*(lj_factor-1.0d0)


                End If

    
            End Do
            
        End Do


      Else if ((initstyle(itype,ibox) .eq. 'simple_cubic') .OR. (initstyle(itype,ibox) .eq. 'random')) Then


        ! Loop over sites on imol
        Do isite = 1, n_sites(itype)

            ! Get isite type 
            isitetype = site_type(isite,itype)

            ! loop over jmol's sites 
            Do jsite = 1, n_sites(jtype)

                ! Get jsite type 
                jsitetype = site_type(jsite,jtype)

                ! Calculate the site separation vector 
                rxijs = rx_s(jsite,jmol,ibox) - rx_s(isite,imol,ibox)
                ryijs = ry_s(jsite,jmol,ibox) - ry_s(isite,imol,ibox)
                rzijs = rz_s(jsite,jmol,ibox) - rz_s(isite,imol,ibox)

                ! Check if hard wall boundary will be applied
                If (lno_pbc) Then
                    Continue
                ! with PBC      
                Else
                    ! Apply minimum image convension
                    rxijs = rxijs - dNINT(rxijs/box(1,ibox))*box(1,ibox)
                    ryijs = ryijs - dNINT(ryijs/box(2,ibox))*box(2,ibox)
                    rzijs = rzijs - dNINT(rzijs/box(3,ibox))*box(3,ibox)

                End If
                        

                ! Calculate the square of the site scalar separation distance
                rijssq = rxijs**2 + ryijs**2 + rzijs**2

                ! Calculate Real space Ewald sum contribution
                If (lewld) Then

                    ! Double check rcelect
                    If (rcelect .ne. MINVAL(box(:,ibox)/2.0d0)) Then
                        Write(*,*) 'ENG_IJ: RCELECT NOT EQUAL TO HALF OF MIN BOX LENGTH'
                        STOP
                    End If

                    ! Check with rcelect
                    if(rijssq .gt. rcelectsq) CYCLE

                              ! For now, unperturbed for adsorbate-adsorbate interaction
                    If (qsq(isitetype,itype,jsitetype,jtype) .ne. 0.0d0) Then
                        rijs = dSQRT(rijssq)
                        engelec = engelec + qsq(isitetype,itype,jsitetype,jtype)*dERFC(alpha*rijs)/rijs
                    End If
                End If

                ! Check with r_cutoff
                If (rijssq .gt. r_cutsq) CYCLE

                ! Check with r_min
                If (rijssq .lt. r_minsq) Then
                    OVERLAP = .TRUE.
                    Return
                End If  

!               Call range_check(WITHIN,rz_s(isite,imol,ibox),rz_s(jsite,jmol,ibox))

                ! Check if to use effective fluid-fluid interactions
                If (lfield .and. (field_type .eq. CG_WALL_FFPW)) Then
                
                    Call eng_cgff(ibox,engcgff,rijssq,isite,imol,jsite,jmol)
                    engdisp = engdisp + engcgff
                    CYCLE 
                               
                End if

                ! Calculate dispersion energy
                If ((potential .eq. LENNARD_JONES) .and. (epsilon(isitetype,itype,jsitetype,jtype) .ne. 0.0d0)) Then

                    ! Calculate lJ factor
                    lj_factor = (sigmasq(isitetype,itype,jsitetype,jtype)/rijssq)**3
                    
                    ! Calculate dispersion energy using 12-6 Lennard-Jones potential
                    engdisp = engdisp + &
                    & 4.0d0*epsilon(isitetype,itype,jsitetype,jtype)*lj_factor*(lj_factor-1.0d0)

                    !Calculate virial and hypervirial
                    If(ldumpvir) Then

                        rdudr = 24.0*epsilon(isitetype,itype,jsitetype,jtype)* &
                                & lj_factor * (1.0 - 2.0*lj_factor)

                        ! ! Molecular Virial in units of [K] (A&T 2.63)
                        ! vir_update = vir_update + 24.0*epsilon(isitetype,itype,jsitetype,jtype)* &
                        !     & lj_factor * (1.0 - 2.0*lj_factor)/rijssq * &
                        !     & (rxij*rxijs + ryij*ryijs + rzij*rzijs)

                        ! Atomic virial in units of [K]
                        vir_update = vir_update + rdudr

                        ! Atomic Hypervirial in units of [K] (A&T 2.78)
                        virhyp_update = virhyp_update + rdudr + 24.0*epsilon(isitetype,itype,jsitetype,jtype)* &
                                & lj_factor * (26.0*lj_factor - 7.0)

                    End If

                End If

                
            ! End loop jmol's sites
            End Do
        ! End loop imol's sites 
        End Do

        
      End If

      
      Return
      
      End Subroutine 
        
      

