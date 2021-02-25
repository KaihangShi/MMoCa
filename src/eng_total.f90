! ==========================================================
! This subroutine is used to calculate the total energy of the 
! system in each box
! Created on 12-26-2016 by Kaihang Shi 
!
! eng_tot (output): total energy of the system
! ==========================================================


      Subroutine eng_total(ibox,eng_tot)

      Use global

      IMPLICIT NONE

      ! Passed
      Integer :: ibox
      Double Precision :: eng_tot

      ! Local
      Integer :: imol, jmol, itype, jtype
      Double Precision :: engdisp, engelec, engewld, engfield, virtot, virhyp


      ! Initialize energies in ibox
      eng_tot = 0.0d0
      energy_disp(ibox) = 0.0d0
      ewld_real(ibox) = 0.0d0
      ewld_tot(ibox) = 0.0d0
      energy_tail(ibox) = 0.0d0
      energy_field(ibox) = 0.0d0

      ! Initialize virial (bc this module recalculate all energies & virials)
      If (ldumpvir) Then
        vir(ibox) = 0.0d0
        vir_hyp(ibox) = 0.0d0
      End If 


      ! Calculate the contribution from the external field
      If (lfield) Then
        
        ! Loop over all of the molecules
        Do imol = 1, n_mol_tot(ibox)
            
            ! Calculate the energy with the external field  
            Call eng_field(ibox,imol,engfield)

            ! Check for overlap
            if(OVERLAP) RETURN

            ! Update external field energy [K]
            energy_field(ibox) = energy_field(ibox) + engfield

        End Do

      End If

      ! Calculate short-range interactions using cell list, scaled with N
      ! Added on July 15, 2019
      If (lclist) Then

        Do imol = 1, n_mol_tot(ibox)

            Call eng_ij_clist(ibox,imol,engdisp,virtot,virhyp)

            ! Check if structure overlaps
            If (OVERLAP) Return

            ! Update dispersion energy of the system in ibox [K]
            ! Consider double counting
            energy_disp(ibox) = energy_disp(ibox) + 0.5d0*engdisp
            ! Update virial [K]
            If(ldumpvir) Then
                vir(ibox) = vir(ibox) + 0.5d0*virtot
                vir_hyp(ibox) = vir_hyp(ibox) + 0.5d0*virhyp
            End If

        End Do

      Else

        ! Calculate the vdw/electrostatic INTER-molecular energy contribution
        ! Loop over all pairs of molecules 
        Do imol = 1, n_mol_tot(ibox)-1
            Do jmol = imol+1, n_mol_tot(ibox)

                ! Calculate the dispersion energy & real space contribution to Ewald sum
                Call eng_ij(ibox,imol,jmol,engdisp,engelec,virtot,virhyp)

                ! Check if structure overlaps
                If (OVERLAP) Then
 !                  Write (*,*) 'FATAL ERROR: Initial/Final STRUCTURE OVERLAPS'
 !                  STOP
                    Return
                End If

                ! Update dispersion energy of the system in ibox [K]
                energy_disp(ibox) = energy_disp(ibox) + engdisp

                ! Update real space Ewald sum [EE]
                if(lewld) ewld_real(ibox) = ewld_real(ibox) + engelec

                ! Update virial [K]
                If(ldumpvir) Then
                    vir(ibox) = vir(ibox) + virtot
                    vir_hyp(ibox) = vir_hyp(ibox) + virhyp
                End If

            End Do
        End Do

      Endif

      ! Convert real space energy unit to [K]
      if(lewld) ewld_real(ibox) = ewld_real(ibox)*EETOK


      ! Calculate vdW tail correction (12-6 Lennard-Jones for now)
      ! Note: Tail correction is a constant as long as r_cut and 
      !       volume (number density or n_mol) remain unchanged
      If (ltailc) Then
        
        ! Loop over all molecule types
        Do itype = 1, n_mol_types
            Do jtype = 1, n_mol_types
                
                ! Calculate tail correction between itype and jtype [K*A^3]
                Call eng_tail(ibox,itype,jtype)

                ! Apply prefactor N*rho and convert to unit of [K]
                energy_tail(ibox) = energy_tail(ibox) + &
                    & n_mol(itype,ibox)*n_mol(jtype,ibox)/vol(ibox)*vdw_tail(itype,jtype)
            End Do
        End Do
      End If

      ! Calculate Fourier space and self/intra correction part of Ewald Sum [K]
      If (lewld) Then

        Call ewld_total(ibox,engewld)

        ! Update ewald energies [K]
        ewld_tot(ibox) = ewld_real(ibox) + engewld
      End If

      ! Update total energy [K]
      eng_tot = energy_disp(ibox) + ewld_tot(ibox) + energy_tail(ibox) + energy_field(ibox)
        
      
        Return
      
      End Subroutine 



