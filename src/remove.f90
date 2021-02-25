! ==========================================================
! This subroutine is used to remove a molecule from a box
! Created on 1-31-2017 by Kaihang Shi
! ==========================================================


      Subroutine remove(ibox)
      
      Use global

      IMPLICIT NONE

      ! Passed 
      Integer :: ibox

      ! Local
      Logical :: ACCEPT
      Integer :: itype, itype_remove, imol, iitype, jjtype
      Double Precision :: ran, boltzf
      Double Precision, Dimension(n_sites_max) :: rx_s_old, ry_s_old, rz_s_old
      Double Precision :: eng_old, eng_change, ewldchange, engtailnew, engtailold, vir_old, virhyp_old
        

      ! Assume initially reject this move
      ACCEPT = .FALSE.

      ! Initialize local variables
      eng_change = 0.0d0

      ! Generate a random number
      ran = random(idum)
      ! Randomly choose a molecule type to remove
      Do itype = 1, n_mol_types
        If (ran .lt. transfer_prob(itype)) Then
            itype_remove = itype
            goto 100
        End If
      End Do

100   Continue

      ! Doubel check: Assume no remove for external structure
      if (initstyle(itype_remove,ibox) .eq. 'coords') Then
        Write(*,*) 'FATAL ERROR: TRIAL REMOVAL MOVE FOR EXTERNAL STRUCTURE'
        STOP
      End If
        
      ! Make sure at least one molecule of itype_remove in the box
      if(n_mol(itype_remove,ibox) .le. 0) Return

      ! Randomly choose a molecule of itype_remove
      Do

        ! Randomly select a molecule 
        imol = INT(random(idum)*n_mol_tot(ibox)) + 1

        ! Exit if it is the right type
        IF(mol_type(imol,ibox) .EQ. itype_remove) EXIT

      END Do


      If (lclist) Then
        Call eng_mol_clist(ibox,imol,eng_old,vir_old,virhyp_old)
      Else
        ! Calculate energy of old configuration
        Call eng_mol(ibox,imol,eng_old,vir_old,virhyp_old)
      End If


      ! Calculate the energy change 
      eng_change = 0.0d0 - eng_old

      ! Update the number of remove trials 
      rem_stat(1,itype_remove) = rem_stat(1,itype_remove) + 1


      ! Check for Ewald Sum 
      If (lewld) Then

        ! Get old site positions (redundant but compiler may complain)
        rx_s_old(:) = rx_s(:,imol,ibox)
        ry_s_old(:) = ry_s(:,imol,ibox)
        rz_s_old(:) = rz_s(:,imol,ibox)
        
        ! Calculate the reciprocal energy change, self/intra correction change & slab correction change
        ! in unit of [K]
        Call ewld_change(3,ibox,imol,rx_s_old,ry_s_old,rz_s_old,ewldchange)

        ! Update total energy change
        eng_change = eng_change + ewldchange

      End If

      ! Check for vdW tail correction
      ! Because n_mol changes, so energy_tail changes accordingly
      If (ltailc) Then

        ! Initialize new tail energy
        engtailnew = 0.0d0

        ! Store old vdW tail energy
        engtailold = energy_tail(ibox)

        ! Temporarily decrease itype_remove's mol number
        n_mol(itype_remove,ibox) = n_mol(itype_remove,ibox) - 1

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

        ! Restore itype_remove's mol number
        n_mol(itype_remove,ibox) = n_mol(itype_remove,ibox) + 1

      End If


      ! Check for large energy change
      if(eng_change/temp .gt. 50.0d0) goto 200


      ! Calculate determinant factor
      boltzf = DBLE(n_mol(itype_remove,ibox))*dEXP(-eng_change/temp)/mol_act(itype_remove)

      ! Apply Metropolis scheme
      If (boltzf .ge. 1.0d0) Then
        ACCEPT = .TRUE.
      Else
        if(random(idum) .lt. boltzf) ACCEPT = .TRUE.
      End If

      
      ! Continue here if overlap and large energy change
200   Continue

      If (ACCEPT) Then
        
        ! Update accepted remove trial counter
        rem_stat(2,itype_remove) = rem_stat(2,itype_remove) + 1

        ! Update energy 
        energy(ibox) = energy(ibox) + eng_change
        energy_tail(ibox) = engtailnew

        ! Update virial
        If(ldumpvir) Then
            vir(ibox) = vir(ibox) - vir_old
            vir_hyp(ibox) = vir_hyp(ibox) - virhyp_old
        End If

        ! Swap the molecule's type with the last molecule in the array
        mol_type(imol,ibox) = mol_type(n_mol_tot(ibox),ibox)

        ! Swap the molecule's position with the last molecule in the array
        rx(imol,ibox) = rx(n_mol_tot(ibox),ibox)
        ry(imol,ibox) = ry(n_mol_tot(ibox),ibox)
        rz(imol,ibox) = rz(n_mol_tot(ibox),ibox)
      
        ! Swap the molecule's orientation with the last molecule in the array
        q1(imol,ibox) = q1(n_mol_tot(ibox),ibox)
        q2(imol,ibox) = q2(n_mol_tot(ibox),ibox)
        q3(imol,ibox) = q3(n_mol_tot(ibox),ibox)
        q4(imol,ibox) = q4(n_mol_tot(ibox),ibox)
        
        ! Compute the site coordinates using the new orientation 
        CALL site_coords(imol,ibox)

        ! Update the Ewald components
        IF(lewld) CALL update_ewld(3,ibox,imol)

        ! Decrease the number of molecules in the box
        n_mol_tot(ibox) = n_mol_tot(ibox) - 1

        ! Decrease the number of molecules of the selected type in the box
        n_mol(itype_remove,ibox) = n_mol(itype_remove,ibox) - 1


        ! Update cell list
        If (lclist) Call update_clist 
            

      Else

        Continue

      End If


        Return
      
      End Subroutine 

























