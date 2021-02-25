! ==========================================================
! This subroutine is used to insert a molecule into a box
! Created on 1-31-2017 by Kaihang Shi
! 
! Modified on Aug 6, 2017:
!   For slit pore geometry, we insert molecules specifically into
!   pores, not to walls (invalid insertion)
! ==========================================================

      Subroutine insert(ibox)

      Use global

      IMPLICIT NONE

      ! Passed
      Integer :: ibox

      ! Local
      Logical :: ACCEPT
      Integer :: itype, itype_insert, imol, iitype, jjtype
      Double Precision :: ran, boltzf
      Double Precision, Dimension(n_sites_max) :: rx_s_old, ry_s_old, rz_s_old
      Double Precision :: eng_new, eng_change, ewldchange, engtailnew, engtailold, vir_new, virhyp_new
      Double Precision :: boxz_access
        

      ! Assume initially reject this move
      ACCEPT = .FALSE.

      ! Initialize local variables
      eng_change = 0.0d0

      ! Generate a random number
      ran = random(idum)
      ! Randomly choose a molecule type to insert
      Do itype = 1, n_mol_types
        If (ran .lt. transfer_prob(itype)) Then
            itype_insert = itype
            goto 100
        End If
      End Do

100   Continue

      ! Doubel check: Assume no insert for external structure
      if (initstyle(itype_insert,ibox) .eq. 'coords') Then
        Write(*,*) 'FATAL ERROR: TRIAL INSERTION MOVE FOR EXTERNAL STRUCTURE'
        STOP
      End If
        
      ! Set molecule's array index 
      imol = n_mol_tot(ibox) + 1

      ! Set imol's type
      mol_type(imol,ibox) = itype_insert

      ! Generate a random position 
      ! Optimized for slit geometry
      If (lfield) Then

        boxz_access = vol_pore(ibox)/(box(1,ibox)*box(2,ibox))

        If ((field_type .eq. STEELE) .or. (field_type .eq. STEELE_SLIT_PORE) .or. &
                & (field_type .eq. STEELE_SLIT_FINITEX)) Then

            ! Generate new position
            rx(imol,ibox) = random(idum)*box(1,ibox)
            ry(imol,ibox) = random(idum)*box(2,ibox)
            rz(imol,ibox) = steele_position(1) + random(idum)*boxz_access

        Else if((field_type .eq. HARD_WALL) .OR. (field_type .eq. HARD_SLIT_FINITEX)) Then

            ! Generate new position
            rx(imol,ibox) = random(idum)*box(1,ibox)
            ry(imol,ibox) = random(idum)*box(2,ibox)
            rz(imol,ibox) = wall_radius + random(idum)*boxz_access

        Else if((field_type .eq. CG_WALL) .or. (field_type .eq. CG_WALL_FFPW) &
                & .or. (field_type .eq. CG_WALL_COS) .or. (field_type .eq. CG_WALL_STRC)) Then

            ! Generate new position
            rx(imol,ibox) = random(idum)*box(1,ibox)
            ry(imol,ibox) = random(idum)*box(2,ibox)
            rz(imol,ibox) = cg_wall_position(1) + random(idum)*boxz_access

        End If

      Else

        ! Normal situation
        rx(imol,ibox) = random(idum)*box(1,ibox)
        ry(imol,ibox) = random(idum)*box(2,ibox)
        rz(imol,ibox) = random(idum)*box(3,ibox)

      End If


      ! Generate a random orientation
      CALL rand_quat(q1(imol,ibox),q2(imol,ibox),q3(imol,ibox),q4(imol,ibox))
              
      ! Calculate site coordinates
      CALL site_coords(imol,ibox) 

      If (lclist) Then
        ! Cell list method
        Call eng_mol_clist(ibox,imol,eng_new,vir_new,virhyp_new)
      Else
        ! Calculate energy of new configuration
        Call eng_mol(ibox,imol,eng_new,vir_new,virhyp_new)
      End If

      
      ! Calculate the energy change 
      eng_change = eng_new - 0.0d0

      ! Update the number of insert trials 
      in_stat(1,itype_insert) = in_stat(1,itype_insert) + 1

      ! Check for overlap
      if(OVERLAP) goto 200

      ! Check for Ewald summation
      IF(lewld) THEN

        ! Initialize dummy arrays (redundant but compiler may complain)
        rx_s_old = 0.0d0
        ry_s_old = 0.0d0
        rz_s_old = 0.0d0

        ! Calculate the Ewald contribution (including reciprocal, self and intra term change)
        CALL ewld_change(2,ibox,imol,rx_s_old,ry_s_old,rz_s_old,ewldchange)

        ! Update the energy change [K]
        eng_change = eng_change + ewldchange

      ENDIF

      ! Check for vdW tail correction
      ! Because n_mol changes, so energy_tail changes accordingly
      If (ltailc) Then

        ! Initialize new tail energy
        engtailnew = 0.0d0

        ! Store old vdW tail energy
        engtailold = energy_tail(ibox)

        ! Temporarily increase itype_insert's mol number
        n_mol(itype_insert,ibox) = n_mol(itype_insert,ibox) + 1

        ! Loop over all molecule types
        Do iitype = 1, n_mol_types
            Do jjtype = 1, n_mol_types

                ! Apply prefactor N*rho and convert to unit of [K]
                engtailnew = engtailnew + &
                & n_mol(iitype,ibox)*n_mol(jjtype,ibox)/vol(ibox)*vdw_tail(iitype,jjtype)
            End Do
        End Do

        ! Update the energy change [K]
        eng_change = eng_change + engtailnew - engtailold

        ! Restore itype_insert's mol number
        n_mol(itype_insert,ibox) = n_mol(itype_insert,ibox) - 1

      End If


      ! Check for large energy change
      if(eng_change/temp .gt. 50.0d0) goto 200


      ! Calculate determinant factor
      boltzf = mol_act(itype_insert)*dEXP(-eng_change/temp)/DBLE(n_mol(itype_insert,ibox)+1)

      ! Apply Metropolis scheme
      If (boltzf .ge. 1.0d0) Then
        ACCEPT = .TRUE.
      Else
        if(random(idum) .lt. boltzf) ACCEPT = .TRUE.
      End If

      
      ! Continue here if overlap and large energy change
200   Continue

      If (ACCEPT) Then
        
        ! Update accepted insert trial counter
        in_stat(2,itype_insert) = in_stat(2,itype_insert) + 1

        ! Update energy 
        energy(ibox) = energy(ibox) + eng_change
        energy_tail(ibox) = engtailnew

        ! Update virial 
        If(ldumpvir) Then
            vir(ibox) = vir(ibox) + vir_new
            vir_hyp(ibox) = vir_hyp(ibox) + virhyp_new
        End If

        ! Increase the number of molecules in the box
        n_mol_tot(ibox) = n_mol_tot(ibox) + 1

        ! Increase the number of molecules of the selected type in the box
        n_mol(itype_insert,ibox) = n_mol(itype_insert,ibox) + 1

        ! Update the Ewald components
        IF(lewld) CALL update_ewld(2,ibox,imol)

        ! Update cell list
        If (lclist) Call update_clist
            

        

      Else

        Continue

      End If
      
        Return
      
      End Subroutine 

















