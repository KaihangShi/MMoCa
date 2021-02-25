! ==========================================================
! This subroutine is used to make a translational move
! Created on 12-29-2016 by Kaihang Shi
! Last modified on Jan, 2020 for virial calculation
! ==========================================================

      Subroutine translate

      Use global

      IMPLICIT NONE

      ! Local
      Logical :: ACCEPT
      Double Precision :: ran, boltzf
      Double Precision :: rx_old, ry_old, rz_old
      Double Precision, Dimension(n_sites_max) :: rx_s_old, ry_s_old, rz_s_old
      Double Precision :: eng_old, eng_new, eng_change, ewldchange, vir_old, vir_new, virhyp_old, virhyp_new
      Integer :: ibox, itype, itype_move, imol



      ! Assume initially reject this move
      ACCEPT = .FALSE.

      ! Randomly select a simulation box
      ibox = INT(random(idum)*n_box) + 1

      ! Generate a random number
      ran = random(idum)
      ! Randomly choose a molecule type to move
      ! Follow detailed balance
      Do itype = 1, n_mol_types
        If (ran .lt. trans_prob(itype)) Then
            itype_move = itype
            goto 100
        End If
      End Do

100   Continue

      ! Doubel check: Assume no move for external structure
      if (initstyle(itype_move,ibox) .eq. 'coords') Then
        Write(*,*) 'FATAL ERROR: TRIAL TRANSLATIONAL MOVE FOR EXTERNAL STRUCTURE'
        STOP
      End If

      ! Make sure at least one molecule of itype_move in the box
      if(n_mol(itype_move,ibox) .le. 0) Return

      ! Randomly choose a molecule of itype_move
      Do

        ! Randomly select a molecule 
        imol = INT(random(idum)*n_mol_tot(ibox)) + 1

        ! Exit if it is the right type
        IF(mol_type(imol,ibox) .EQ. itype_move) EXIT

      END Do


      ! Store old COM position
      rx_old = rx(imol,ibox)
      ry_old = ry(imol,ibox)
      rz_old = rz(imol,ibox)

      ! Store old sites position
      rx_s_old(:) = rx_s(:,imol,ibox)
      ry_s_old(:) = ry_s(:,imol,ibox)
      rz_s_old(:) = rz_s(:,imol,ibox)

      ! cell list
      If (lclist) Then
        Call eng_mol_clist(ibox,imol,eng_old,vir_old,virhyp_old)
      Else
        ! Calculate energy of the old configuration
        Call eng_mol(ibox,imol,eng_old,vir_old,virhyp_old)
      End If


      ! Generate a new position
      rx(imol,ibox) = rx(imol,ibox) + (2.0d0*random(idum) - 1.0d0)*max_trans(itype_move,ibox)
      ry(imol,ibox) = ry(imol,ibox) + (2.0d0*random(idum) - 1.0d0)*max_trans(itype_move,ibox)
      rz(imol,ibox) = rz(imol,ibox) + (2.0d0*random(idum) - 1.0d0)*max_trans(itype_move,ibox)

      ! Check boundary condition
      If (lno_pbc) Then

        ! Hard wall boundary
        if((rx(imol,ibox) .LT. 0.0d0) .or. (rx(imol,ibox) .GT. box(1,ibox)) .or. &
                & (ry(imol,ibox) .LT. 0.0d0) .or. (ry(imol,ibox) .GT. box(2,ibox)) .or. &
                & (rz(imol,ibox) .LT. 0.0d0) .or. (rz(imol,ibox) .GT. box(3,ibox))) goto 200

      ! with PBC
      Else
      
        ! Put NEW coordinates into central box 
        rx(imol,ibox) = rx(imol,ibox) - FLOOR(rx(imol,ibox)/box(1,ibox))*box(1,ibox) 
        ry(imol,ibox) = ry(imol,ibox) - FLOOR(ry(imol,ibox)/box(2,ibox))*box(2,ibox) 
        rz(imol,ibox) = rz(imol,ibox) - FLOOR(rz(imol,ibox)/box(3,ibox))*box(3,ibox) 


      End If

      
      ! Compute the new sites coordinates
      Call site_coords(imol,ibox)


      ! Cell list
      If (lclist) Then
        Call eng_mol_clist(ibox,imol,eng_new,vir_new,virhyp_new)
      Else
        ! Calculate dispersion energy of the new configuration
        Call eng_mol(ibox,imol,eng_new,vir_new,virhyp_new)

      End If

     

      ! Update the number of trial translational move 
      trans_stat(1,itype_move,ibox) = trans_stat(1,itype_move,ibox) + 1

      ! Check if overlap
      If(OVERLAP) goto 200

      ! Calculate energy change [K] due to translational move
      eng_change = eng_new - eng_old

      ! Check for Ewald Sum 
      If (lewld) Then
        
        ! Calculate the reciprocal energy change, self/intra correction change & slab correction change
        ! in unit of [K]
        Call ewld_change(1,ibox,imol,rx_s_old,ry_s_old,rz_s_old,ewldchange)

        ! Update total energy change
        eng_change = eng_change + ewldchange

      End If


      ! Check for large energy change
      if(eng_change/temp .gt. 50.0d0) goto 200


      ! Apply the Metropolis criteria (NVT)
      If (eng_change .le. 0.0d0) Then
        ACCEPT = .TRUE.
      Else
        boltzf = dEXP(-eng_change/temp)
        if(random(idum) .lt. boltzf) ACCEPT = .TRUE.
      End If

      
      ! Continue here if overlap and large energy change
200   Continue

      If (ACCEPT) Then

        ! Update the number of accepted translation
        trans_stat(2,itype_move,ibox) = trans_stat(2,itype_move,ibox) + 1

        ! Update total energy  [K]
        energy(ibox) = energy(ibox) + eng_change

        ! Update total dispersion energy [K]
!       energy_disp(ibox) = energy_disp(ibox) + eng_change 

        ! Update virial [K]
        If (ldumpvir) Then
            vir(ibox) = vir(ibox) -vir_old + vir_new
            vir_hyp(ibox) = vir_hyp(ibox) -virhyp_old + virhyp_new
        End If
        
        ! Update Ewald component
        if(lewld) Call update_ewld(1,ibox,imol)

        ! Update cell list
        If (lclist) Call update_clist
            


      Else

        ! Restore old COM coordinates
        rx(imol,ibox) = rx_old
        ry(imol,ibox) = ry_old
        rz(imol,ibox) = rz_old

        ! Restore old sites coordinates
        rx_s(:,imol,ibox) = rx_s_old(:)
        ry_s(:,imol,ibox) = ry_s_old(:)
        rz_s(:,imol,ibox) = rz_s_old(:)

    
      End If
      
        Return
      
      End Subroutine 














