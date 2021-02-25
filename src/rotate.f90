! ==========================================================
! This subroutine is used to do a rotational move for a 
! molecule
! Created on 1/6/2017 by Kaihang Shi
! Last modified in Jan, 2020
! ==========================================================

      Subroutine rotate

      Use global

      IMPLICIT NONE

      ! Local
      Logical :: ACCEPT
      Double Precision :: ran, boltzf, norm
      Double Precision :: eng_old, eng_new, eng_change, ewldchange, vir_new, vir_old, virhyp_old, virhyp_new
      Integer :: ibox, itype, itype_rotat, imol
      Double Precision :: q1_old, q2_old, q3_old, q4_old
      Double Precision :: q1r, q2r, q3r, q4r
      Double Precision, Dimension(n_sites_max) :: rx_s_old, ry_s_old, rz_s_old


      ! Assume initially reject this move
      ACCEPT = .FALSE.

      ! Randomly select a simulation box
      ibox = INT(random(idum)*n_box) + 1

      ! Generate a random number
      ran = random(idum)
      ! Randomly choose a molecule type to rotate
      ! Follow detailed balance
      Do itype = 1, n_mol_types
        If (ran .lt. rotat_prob(itype)) Then
            itype_rotat = itype
            goto 100
        End If
      End Do

100   Continue

      ! Doubel check: Assume no move for external structure
      if (initstyle(itype_rotat,ibox) .eq. 'coords') Then
        Write(*,*) 'FATAL ERROR: TRIAL ROTATIONAL MOVE FOR EXTERNAL STRUCTURE'
        STOP
      End If

      ! Make sure itype_rotat is not a single site molecule
      if(n_sites(itype_rotat) .eq. 1) Return

      ! Make sure at least one molecule of itype_rotat in the box
      if(n_mol(itype_rotat,ibox) .le. 0) Return


      ! Randomly choose a molecule of itype_rotat
      Do

        ! Randomly select a molecule 
        imol = INT(random(idum)*n_mol_tot(ibox)) + 1

        ! Exit if it is the right type
        IF(mol_type(imol,ibox) .EQ. itype_rotat) EXIT

      END Do


      ! Store old quaternions
      q1_old = q1(imol,ibox)
      q2_old = q2(imol,ibox)
      q3_old = q3(imol,ibox)
      q4_old = q4(imol,ibox)

      ! Store old sites position
      rx_s_old(:) = rx_s(:,imol,ibox)
      ry_s_old(:) = ry_s(:,imol,ibox)
      rz_s_old(:) = rz_s(:,imol,ibox)


      If (lclist) Then
        Call eng_mol_clist(ibox,imol,eng_old,vir_old,virhyp_old)
      Else
        ! Calculate dispersion energy of the old configuration
        Call eng_mol(ibox,imol,eng_old,vir_old,virhyp_old)
      End If

      

      ! Generate a random quaternion
      Call rand_quat(q1r,q2r,q3r,q4r)

      ! Generate a new orientation
      q1(imol,ibox) = q1_old + max_rotat(itype_rotat,ibox) * q1r
      q2(imol,ibox) = q2_old + max_rotat(itype_rotat,ibox) * q2r
      q3(imol,ibox) = q3_old + max_rotat(itype_rotat,ibox) * q3r
      q4(imol,ibox) = q4_old + max_rotat(itype_rotat,ibox) * q4r

      ! Calculate the norm of the new quaternion
      norm = dSQRT(q1(imol,ibox)**2 + q2(imol,ibox)**2 + q3(imol,ibox)**2 + q4(imol,ibox)**2)  

      ! Normalize the new quaternion
      q1(imol,ibox) = q1(imol,ibox)/norm
      q2(imol,ibox) = q2(imol,ibox)/norm
      q3(imol,ibox) = q3(imol,ibox)/norm
      q4(imol,ibox) = q4(imol,ibox)/norm

      ! Compute the new sites coordinates
      Call site_coords(imol,ibox)

      If (lclist) Then
        Call eng_mol_clist(ibox,imol,eng_new,vir_new,virhyp_new)
      Else
        ! Calculate dispersion energy of the new configuration
        Call eng_mol(ibox,imol,eng_new,vir_new,virhyp_new)
      End If  

      ! Update the number of trial rotational move 
      rotat_stat(1,itype_rotat,ibox) = rotat_stat(1,itype_rotat,ibox) + 1

      ! Check if overlap
      If(OVERLAP) goto 200

      ! Calculate energy change [K] due to rotational move
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

        ! Update the number of accepted rotation
        rotat_stat(2,itype_rotat,ibox) = rotat_stat(2,itype_rotat,ibox) + 1

        ! Update total energy  [K]
        energy(ibox) = energy(ibox) + eng_change

        ! Update total dispersion energy [K]
!       energy_disp(ibox) = energy_disp(ibox) + eng_change

        ! Update ewld components
        if(lewld) Call update_ewld(1,ibox,imol)

        ! Update virial [K]
        If(ldumpvir) Then
            vir(ibox) = vir(ibox) - vir_old + vir_new
            vir_hyp(ibox) = vir_hyp(ibox) - virhyp_old + virhyp_new
        End If


      Else

        ! Restore old quaternion
        q1(imol,ibox) = q1_old
        q2(imol,ibox) = q2_old
        q3(imol,ibox) = q3_old
        q4(imol,ibox) = q4_old

        ! Restore old sites coordinates
        rx_s(:,imol,ibox) = rx_s_old(:)
        ry_s(:,imol,ibox) = ry_s_old(:)
        rz_s(:,imol,ibox) = rz_s_old(:)

    
      End If

      
      Return
      
      End Subroutine 












