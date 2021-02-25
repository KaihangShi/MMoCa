! ==========================================================
! This subroutine is used to do a volume change move in NPT
! MC simulation
! Created on 1-17-2017 by Kaihang Shi
!
! Reference:
! [1]. Eppinga and Frenkel. Monte Carlo study of the isotropic 
!      nad nematic phases of infinitely thin hard platelets. 
!      Mol. Phys., 52:1303-1334, 1984. 
! Comments: we didn't scale r_cut here for vdW interaction, 
!           so cutoff radius should maintain a value which makes 
!           longer-distance interactions negligible. 
!           But we should recalculate rcelect for real space coulomb interaction
!           in 'ewald_fix_kmax' scheme.
! ==========================================================


      Subroutine volchange_npt(ibox)

      Use global

      IMPLICIT NONE

      ! Passed
      Integer :: ibox

      ! Local 
      Logical :: ACCEPT
      COMPLEX*16,DIMENSION(:,:,:,:),ALLOCATABLE :: eikx_old, eiky_old, eikz_old
      COMPLEX*16,DIMENSION(:,:),ALLOCATABLE :: skewld_old
      Double Precision, Dimension(:), Allocatable :: ewld_surfc_old
      Double Precision :: eng_old, vol_old, box_old, eng_new, engtail_old, vir_old, virhyp_old
      Double Precision :: lnvol, scalef, arg, boltzf
      Integer :: idirec, imol
      Integer :: ierr


      ! Assume initially reject this move
      ACCEPT = .FALSE.

      ! Store old configuration properties
      eng_old = energy(ibox)
      vol_old = vol(ibox)
      box_old = box(1,ibox)
      engtail_old = energy_tail(ibox)
      If (ldumpvir) Then
        vir_old = vir(ibox)
        virhyp_old = vir_hyp(ibox)
      End If 

      ! Compute random volume displacement in ln(V) (Eppinga and Frenkel)
      lnvol = dLOG(vol_old) + (random(idum) - 0.5d0)*max_vol(ibox)

      ! Compute the new volumes
      vol(ibox) = dEXP(lnvol)

      ! Compute the new box dimension 
      ! Assume cubic box 
      Do idirec = 1,3
        box(idirec,ibox) = vol(ibox)**(1.0d0/3.0d0)
      End Do

      ! Check cutoff radius. cutoff raidus for vdW interaction remain unchanged during volume change move
      If (r_cut .gt. 0.5d0*box(1,ibox)) Then
        Write(*,*) 'VOLCHANGE_NPT: CUTOFF RADIUS IS LARGER THAN HALF OF BOX LENGTH. BAD BOX SETTING.'
        STOP
      End If

      ! Calculate scale factor 
      scalef = box(1,ibox)/box_old
      
      ! Loop over the molecules in the ibox
      ! Assume no external structures (coords type) in NPT simulation
      Do imol = 1, n_mol_tot(ibox)
        
        ! Scale the molecules center of mass (NEW)
        rx(imol,ibox) = rx(imol,ibox)*scalef
        ry(imol,ibox) = ry(imol,ibox)*scalef
        rz(imol,ibox) = rz(imol,ibox)*scalef

        ! Calculate the site position
        Call site_coords(imol,ibox)

      End Do

      ! Check for ewald sum and reset ewald parameters
      If (lewld) Then
        ! Initialize the error handler
        ierr = 0

        ! Allocate the temporary arrays
        ALLOCATE(eikx_old(n_box,n_mol_max,n_sites_max,0:MAXVAL(k_max)), STAT=ierr)
        ALLOCATE(eiky_old(n_box,n_mol_max,n_sites_max,-MAXVAL(k_max):MAXVAL(k_max)), STAT=ierr)
        ALLOCATE(eikz_old(n_box,n_mol_max,n_sites_max,-MAXVAL(k_max):MAXVAL(k_max)), STAT=ierr)
        ALLOCATE(skewld_old(n_box,MAXVAL(maxk)), STAT=ierr)
        ALLOCATE(ewld_surfc_old(n_box), STAT=ierr)

        ! Check for allocation error
        IF(ierr .NE. 0) THEN
          WRITE(*,*) 'FATAL ERROR: ALLOCATION OF ARRAYS FAILED IN SUBROUTINE VOLCHANGE_NPT'
          STOP
        ENDIF

        ! Store the old variables
        eikx_old = eikx
        eiky_old = eiky
        eikz_old = eikz
        skewld_old = skewld
        ewld_surfc_old(ibox) = ewld_surfc(ibox)

        ! Reset basic Ewald sum parameters(alpha, rcelect etc.)
        Call set_ewld(2,ibox)

        ! Setup the Ewald vectors using the new box dimensions   
        CALL set_ewld(1,ibox)
        
      End If

      ! Cell list
      If (lclist) Then
        ! Reset basic parameters of the cell list
        Call set_clist
        Call update_clist
      End If

      ! Calculate the total energy for ibox
      Call eng_total(ibox,eng_new)

      ! Update volume change trial array
      vol_stat(1,ibox) = vol_stat(1,ibox) + 1

      ! Check for overlap due to shrink of box
      If (OVERLAP) goto 100

      ! Compute the argument factor
      arg = ((eng_new-eng_old) + &
        & PVTOK*press*(vol(ibox)-vol_old) - &
        & temp*DBLE(n_mol_tot(ibox)+1)*dLOG(vol(ibox)/vol_old))

      ! Apply the Metropolis criteria
      If (arg .le. 0.0d0) Then
        ACCEPT = .TRUE.
      Else
        boltzf = dEXP(-arg/temp)
        if(random(idum) .lt. boltzf) ACCEPT = .TRUE.
      End If


100   If (ACCEPT) Then
        
        ! Update volume change acceptance counter
        vol_stat(2,ibox) = vol_stat(2,ibox) + 1

        ! Update total energy [K]
        energy(ibox) = eng_new

      Else

        ! Restore the old box volume
        vol(ibox) = vol_old

        ! Restore the old box dimensions
        DO idirec = 1,3
          box(idirec,ibox) = box_old 
        ENDDO

        ! Loop over the molecules in the box
        DO imol = 1,n_mol_tot(ibox)
        
          ! Rescale the molecules center of mass to original positions
          rx(imol,ibox) = rx(imol,ibox)/scalef
          ry(imol,ibox) = ry(imol,ibox)/scalef
          rz(imol,ibox) = rz(imol,ibox)/scalef

          ! Calculate the site position
          CALL site_coords(imol,ibox)

        ENDDO

        ! Restore tail correction for use in insert/widom/remove subroutines
        energy_tail(ibox) = engtail_old

        ! Restore old virial
        If(ldumpvir) Then
            vir(ibox) = vir_old
            vir_hyp(ibox) = virhyp_old
        End If

        ! Reset the Ewald vectors
        IF(lewld) THEN

          ! Restore the old variables
          eikx = eikx_old
          eiky = eiky_old
          eikz = eikz_old
          skewld = skewld_old
          ewld_surfc(ibox) = ewld_surfc_old(ibox)

          ! Reset basic Ewald Sum parameters (alpha, rcelect)
          Call set_ewld(2,ibox)

          ! Reset the Ewald vectors using the old box dimensions   
          CALL set_ewld(1,ibox)

        ENDIF

        ! Reset cell list
        If (lclist) Then
            Call set_clist
            Call update_clist
        End If

        
      ENDIF

      ! Deallocate the Ewald arrays
      IF(lewld) DEALLOCATE(eikx_old,eiky_old,eikz_old,skewld_old,ewld_surfc_old)
    
      
      Return
      
      End Subroutine 














