! ==========================================================
! This subroutine is used to calculate energy and virial 
! of imol using the cell list 
! Created on 7-15-2019 by Kaihang Shi
! For now, cell list is not compatible with Ewald (fix_kmax)
! ==========================================================


      Subroutine eng_ij_clist(ibox,imol,eng_update,vir_update,virhyp_update)

      Use global

      IMPLICIT NONE

      ! Passed
      Integer :: ibox, imol
      Double Precision :: eng_update, vir_update, virhyp_update

      ! Local
      Integer :: jmol, ibinx, ibiny, ibinz, icel, jcel, inei, itype, jtype, isite, jsite, isitetype, jsitetype
      Double Precision :: engdisp, engelec, lj_factor, rdudr
      Double Precision :: rxijs, ryijs, rzijs, rijssq
      Double Precision :: rxij, ryij, rzij

      ! Assume no overlap at the beginning
      OVERLAP = .False.

      ! Initialize variables
      eng_update = 0.0d0
      engdisp = 0.0d0
      engelec = 0.0d0
      vir_update = 0.0d0
      virhyp_update = 0.0d0

      ! Get the molecule types 
      itype = mol_type(imol,ibox)

      ! Check initial style
      ! If initstyle(itype,ibox) = 'coords', we assume this molecule type is external structure:
      ! open surface or pore structure. No calculation of COM separation is needed
      ! If 'coords' is specified, we assume imol = 1 corresponds to the ONLY external structure molecule 
      If (initstyle(itype,ibox) .eq. 'coords') Then

        ! Directly loop over subsrate sites
        Do isite = 1, n_sites(itype)

            ! Get isite type 
            isitetype = site_type(isite,itype)

            ! Determine the location of isite
            ibinx = FLOOR(rx_s(isite,imol,ibox)/clist_dx) + 1
            ibiny = FLOOR(ry_s(isite,imol,ibox)/clist_dy) + 1
            ibinz = FLOOR(rz_s(isite,imol,ibox)/clist_dz) + 1

            ! Locate cell number
            icel = clist_loca(ibinx,ibiny,ibinz)

            ! Loop over neighboring cells (including icel itself)
            Do inei = 1, 27

                jcel = clist_neigh(icel,inei)
                ! only substrate-adsorbate interaction is included
                ! No substrate self interactions
                jmol = clist_hoc(jcel)

                Do while (jmol .ne. 0)  

                    ! Get jmol's type
                    jtype = mol_type(jmol,ibox)

                    ! loop over jmol's sites (adsorbates)
                    Do jsite = 1, n_sites(jtype)

                        ! Get jsite type 
                        jsitetype = site_type(jsite,jtype)

                        ! Calculate the site separation vector 
                        rxijs = rx_s(jsite,jmol,ibox) - rx_s(isite,imol,ibox)
                        ryijs = ry_s(jsite,jmol,ibox) - ry_s(isite,imol,ibox)
                        rzijs = rz_s(jsite,jmol,ibox) - rz_s(isite,imol,ibox)

                        ! Solid boundary, no periodic boundary conditions
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
                    Enddo

                    ! Update jmol using linked-list
                    jmol = clist_llist(jmol)

                EndDo
            ! end inei cell
            EndDo
        ! end looping over substrate sites
        EndDo

        ! Update energy
        eng_update = eng_update + engdisp


      ! Adsorbate
      Else If ((initstyle(itype,ibox) .eq. 'simple_cubic') .OR. (initstyle(itype,ibox) .eq. 'random')) Then

        ! Determine the location of imol
        ibinx = FLOOR(rx(imol,ibox)/clist_dx) + 1
        ibiny = FLOOR(ry(imol,ibox)/clist_dy) + 1
        ibinz = FLOOR(rz(imol,ibox)/clist_dz) + 1

        ! Locate cell number
        icel = clist_loca(ibinx,ibiny,ibinz)
      
        ! Loop over neighboring cells (including icel itself)
        Do inei = 1, 27

            jcel = clist_neigh(icel,inei)
            
            ! -----------------------------------------
            ! Evaluate adsorbate-adsorbate interactions
            jmol = clist_hoc(jcel)

            Do while (jmol .ne. 0)  

                If (imol .eq. jmol) goto 100

                ! Get jmol's type
                jtype = mol_type(jmol,ibox)

                ! If (ldumpvir) Then

                !   ! Calculate the COM distance
                !   rxij = rx(jmol,ibox) - rx(imol,ibox)
                !   ryij = ry(jmol,ibox) - ry(imol,ibox)
                !   rzij = rz(jmol,ibox) - rz(imol,ibox)

      !               If (lno_pbc) Then
      !                   Continue
      !               Else
      !                   ! Apply minimum image convension
      !                   rxij = rxij - dNINT(rxij/box(1,ibox))*box(1,ibox)
      !                   ryij = ryij - dNINT(ryij/box(2,ibox))*box(2,ibox)
      !                   rzij = rzij - dNINT(rzij/box(3,ibox))*box(3,ibox)
                              
      !               End If
      !           End If


                ! Directly loop over imol's sites
                Do isite = 1, n_sites(itype)

                    ! Get isite type 
                    isitetype = site_type(isite,itype)

                    ! loop over jmol's sites (adsorbates)
                    Do jsite = 1, n_sites(jtype)

                        ! Get jsite type 
                        jsitetype = site_type(jsite,jtype)

                        ! Calculate the site separation vector 
                        rxijs = rx_s(jsite,jmol,ibox) - rx_s(isite,imol,ibox)
                        ryijs = ry_s(jsite,jmol,ibox) - ry_s(isite,imol,ibox)
                        rzijs = rz_s(jsite,jmol,ibox) - rz_s(isite,imol,ibox)

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

                    Enddo
                EndDo

                ! Update jmol using linked-list
100             jmol = clist_llist(jmol)

            EndDo

            ! ----------------------------------------
            ! Evaluate subsrate-adsorbate interactions
            jsite = clist_hoc_sub(jcel)

            Do while (jsite .ne. 0) 

                ! Get jsite type 
                ! jtype = 1 for external structure (hard-coded)
                jtype = 1
                jmol = 1
                jsitetype = site_type(jsite,jtype)

                ! Directly loop over imol's sites
                Do isite = 1, n_sites(itype)

                    ! Get isite type 
                    isitetype = site_type(isite,itype)
                
                    ! Calculate the site separation vector 
                    rxijs = rx_s(jsite,jmol,ibox) - rx_s(isite,imol,ibox)
                    ryijs = ry_s(jsite,jmol,ibox) - ry_s(isite,imol,ibox)
                    rzijs = rz_s(jsite,jmol,ibox) - rz_s(isite,imol,ibox)

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
                    
                EndDo

                ! Update jsite using linked-list
                jsite = clist_llist_sub(jsite)

            EndDo

        ! End looping over neighboring cells
        EndDo


        ! Update energy
        eng_update = eng_update + engdisp


      EndIF

        
      
      Return
      
      End Subroutine 











