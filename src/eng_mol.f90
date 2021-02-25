! ==========================================================
! This subroutine is used to calculate energy & virial 
! of imol in the system in the box (real space)
! Created on 12-30-2016 by Kaihang Shi
! Last modified on Jan 2020 for virial calculations
! ==========================================================


      Subroutine eng_mol(ibox,imol,eng_update,vir_update,virhyp_update)

      Use global

      IMPLICIT NONE

      ! Passed
      ! imol's initstyle is 'simple_cubic'
      Integer :: ibox, imol
      Double Precision :: eng_update, vir_update, virhyp_update

      ! Local
      Integer :: jmol
      Double Precision :: engdisp, engelec, engfield, virij, virhypij


      ! Initialize variables
      eng_update = 0.0d0
      vir_update = 0.0d0
      virhyp_update = 0.0d0

      ! Calculate the contribution from the external field
      If (lfield) Then
            
        ! Calculate the energy with the external field  
        Call eng_field(ibox,imol,engfield)

        ! Check for overlap
        if(OVERLAP) RETURN

        ! Update energy 
        eng_update = eng_update + engfield

      End If

      ! Loop over other molecule
      ! Assume external structure (jmol=1, if any) will not move
      Do jmol = 1, n_mol_tot(ibox)
        
        ! Avoid interation with itself
        if(imol == jmol) CYCLE

        ! Calculate dispersion energy & real space Coulombic energy
        ! Be careful about argument position of imol and jmol
        Call eng_ij(ibox,jmol,imol,engdisp,engelec,virij,virhypij)

        ! Check with overlap
        If (OVERLAP) Return

        ! Dispersion contribution and apply prefactor of LJ potential
        eng_update = eng_update + engdisp

        ! Update energy of electrostaic contribution and convert unit from [EE] to [K]
        if(lewld) eng_update = eng_update + EETOK*engelec

        ! Accumulate virial and hypervirial [K]
        If (ldumpvir) Then
            vir_update = vir_update + virij
            virhyp_update = virhyp_update + virhypij
        End If


      End Do
      
        
      
      Return
      
      End Subroutine 