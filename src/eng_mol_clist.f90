! ==========================================================
! This subroutine is used to calculate energy and virial of imol using 
!   the cell list 
! Created on 7-15-2019 by Kaihang Shi
! For now, cell list is not compatible with Ewald (ewld_fix_kmax)
! ==========================================================


      Subroutine eng_mol_clist(ibox,imol,eng_update,vir_update,virhyp_update)

      Use global

      IMPLICIT NONE

      ! Passed
      Integer :: ibox, imol
      Double Precision :: eng_update, vir_update, virhyp_update

      ! Local
      Integer :: jmol, ibinx, ibiny, ibinz, icel, jcel, inei, itype, jtype, isite, jsite, isitetype, jsitetype
      Double Precision :: engdisp, engfield, lj_factor
      Double Precision :: rxijs, ryijs, rzijs, rijssq


      ! Initialize variables
      eng_update = 0.0d0
      

      ! Calculate the contribution from the external field
      If (lfield) Then
            
        ! Calculate the energy with the external field  
        Call eng_field(ibox,imol,engfield)

        ! Check for overlap
        if(OVERLAP) RETURN

        ! Update energy 
        eng_update = eng_update + engfield

      End If


      ! Using cell list to calculate imol's interactions with other molecules
      Call eng_ij_clist(ibox,imol,engdisp,vir_update,virhyp_update)

      ! Check with overlap
      If (OVERLAP) Return

      ! Dispersion contribution
      eng_update = eng_update + engdisp



      Return
      
      End Subroutine 









