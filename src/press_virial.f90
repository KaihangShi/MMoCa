! ==========================================================
! This subroutine is used to calculate
! the pressure tensor for general planar surface geometry from virial route
! using Irving-Kirkwood and/or Harasima definition of contour
! Adapted from Deepti's subroutine
! Created on March 1 2019 by Kaihang Shi
! Last update: March, 2019
! ==========================================================



SUBROUTINE press_virial(selector,iblock,istep)

USE global
IMPLICIT NONE

! Passed
Integer :: selector, iblock, istep

!Local variables
INTEGER :: imol, jmol, i, j, ibin, ibox
INTEGER :: isite, jsite, itype, jtype, isitetype, jsitetype

DOUBLE PRECISION :: rxi, ryi, rzi, xij, yij, zij, xijsq, yijsq, zijsq, rij, rxj, ryj, rzj
DOUBLE PRECISION :: rxis, ryis, rzis, rxjs, ryjs, rzjs, sig_s, eps_s
DOUBLE PRECISION :: zijbottom, zijtop, z, x2y2tr, z2tr, a, b, c
DOUBLE PRECISION :: sr2, sr3, sr4, sr6, sr10, sr12, phi, phitz
DOUBLE PRECISION :: pnw_1, pnw_2, pnc_1, ptc_1, ptc_2, dudr
DOUBLE PRECISION :: rijs, rijstop, rijsbottom, xijs, yijs, zijs, zik1, zik2, zjk1, zjk2
DOUBLE PRECISION :: zijstop, zijsbottom, xijssq, yijssq, zijssq
Double precision :: posz, zii, zjj
Double Precision :: area


Select Case (selector)
  
  ! Perform one-time calculations of pressure tensor from virial route
  ! Input: iblock
  Case(1)
    
   ! Loop over boxes
   Do ibox = 1, n_box

     ! Check if the calculation should be attempted
     IF((virialpress_freq .EQ. 0) .OR. (MOD(istep,virialpress_freq) .NE. 0)) CYCLE

     ! Update counter (double precision)
     virialpress_stat(ibox) = virialpress_stat(ibox) + 1.0d0

     ! Calculate area in xy direction 
     ! For NPT ensemble, this value is changing along the simultion
     dz = box(3,1)/zden_bins
     area = box(1,ibox)*box(2,ibox)

     !Loop over molecules
     DO imol=1,n_mol_tot(ibox)

        !Get the molecule type
        itype = mol_type(imol,ibox)

        !Loop over molecule sites
        DO isite=1,n_sites(itype)

          ! Get isite type 
          isitetype = site_type(isite,itype)

          !Get x,y,z position of site
          rxis = rx_s(isite,imol,ibox)
          ryis = ry_s(isite,imol,ibox)
          rzis = rz_s(isite,imol,ibox)

          !Fluid-fluid interaction, avoiding double count of i,j and j,i so just
          !loop by using i<j     
          IF (imol .LT. n_mol_tot(ibox)) THEN
           DO jmol = imol+1, n_mol_tot(ibox)

           !Get molecule type
           jtype = mol_type(jmol,ibox)

           !Get x,y,z position of molecule
           rxj = rx(jmol,ibox)
           ryj = ry(jmol,ibox)
           rzj = rz(jmol,ibox)

           !Loop over molecule sites
           DO jsite=1,n_sites(jtype)

             ! Get jsite type
             jsitetype = site_type(jsite,jtype)

             !Get x,y,z position of molecule
             rxjs = rx_s(jsite,jmol,ibox)
             ryjs = ry_s(jsite,jmol,ibox)
             rzjs = rz_s(jsite,jmol,ibox)


             !Calculate vector between sites
             xijs = rxjs - rxis
             yijs = ryjs - ryis
             zijs = rzjs - rzis
             
             !Apply minimum image convention
             xijs = xijs - dNINT(xijs/box(1,ibox))*box(1,ibox)
             yijs = yijs - dNINT(yijs/box(2,ibox))*box(2,ibox)
             zijs = zijs - dNINT(zijs/box(3,ibox))*box(3,ibox)


             ! Recalculate the z-position of the particle j by keeping particle i at its original place
             ! This is used for particles interacting across the boundary and 
             !  fixing the pressure abnormality near the box boundary
             zjj = rzis + zijs

             
             !Square the values
             xijssq=xijs*xijs
             yijssq=yijs*yijs
             zijssq=zijs*zijs

             !Calculate distance between i and j sites
             rijs=dsqrt(xijssq+yijssq+zijssq)

             ! Check with r_cutoff
             If (rijs .gt. r_cut) CYCLE

             !Calculate local parameters
             x2y2tr=(xijssq+yijssq)/rijs
             z2tr = zijssq/rijs

             !Calculate force (based on sites)
             sr3=(sigma(isitetype,itype,jsitetype,jtype)/rijs)**3
             sr6=sr3*sr3
             sr12=sr6*sr6
             phitz = 24.0*epsilon(isitetype,itype,jsitetype,jtype)/rijs
             dudr = phitz*(sr6 - 2.0d0*sr12)

             ! Definition of integral contour
             ! Only Irving-Kirkwood definition or both
             If (virialpress_ctype .ne. 2) Then

              !Calculate configurational contribution
              ptc_1=x2y2tr*dudr/(dABS(zijs)*area)
              pnc_1=z2tr*dudr/(dABS(zijs)*area)

              !Loop through the bins
              DO ibin= 1, zden_bins

                
                !Get z-position of bin
                z = (DBLE(ibin)-0.5d0)*dz 


                ! Special treatment for particles interacting across the boundary
                ! For j particle lower than bottom surface
                If (zjj .LT. 0.0d0) Then
                  !Unit step function
                  IF (z .LT. rzis) THEN

                    virialpress_pn(ibin,iblock,ibox) = virialpress_pn(ibin,iblock,ibox) + pnc_1
                    virialpress_pt(1,ibin,iblock,ibox) = virialpress_pt(1,ibin,iblock,ibox) + ptc_1 

                  ! End Heavise unit step function
                  ENDIF

                  !Unit step function
                  IF ((z-box(3,1)) .GT. zjj) THEN

                    virialpress_pn(ibin,iblock,ibox) = virialpress_pn(ibin,iblock,ibox) + pnc_1
                    virialpress_pt(1,ibin,iblock,ibox) = virialpress_pt(1,ibin,iblock,ibox) + ptc_1 

                  ! End Heavise unit step function
                  ENDIF

                  CYCLE

                ! For particle j higher than top surface
                Else If (zjj .GE. box(3,1)) Then
                  
                  !Unit step function
                  IF (z .GT. rzis) THEN

                    virialpress_pn(ibin,iblock,ibox) = virialpress_pn(ibin,iblock,ibox) + pnc_1
                    virialpress_pt(1,ibin,iblock,ibox) = virialpress_pt(1,ibin,iblock,ibox) + ptc_1 

                  ! End Heavise unit step function
                  ENDIF

                  !Unit step function
                  IF ((z+box(3,1)) .LT. zjj) THEN

                    virialpress_pn(ibin,iblock,ibox) = virialpress_pn(ibin,iblock,ibox) + pnc_1
                    virialpress_pt(1,ibin,iblock,ibox) = virialpress_pt(1,ibin,iblock,ibox) + ptc_1 

                  ! End Heavise unit step function
                  ENDIF

                  CYCLE
                
                End If

                ! Deal with particles interacting within the box
                !Unit step function
                IF ((z - rzis)/zijs .GT. 0.0d0) THEN
                  IF ((rzjs - z)/zijs .GT.  0.0d0) THEN

                    virialpress_pn(ibin,iblock,ibox) = virialpress_pn(ibin,iblock,ibox) + pnc_1
                    virialpress_pt(1,ibin,iblock,ibox) = virialpress_pt(1,ibin,iblock,ibox) + ptc_1 
                     
                  ENDIF 
                ! End Heavise unit step function
                ENDIF
              !End bin cycle
              ENDDO
             ! End IK 
             ENDIF

             ! Harasima definition
             ! Normal pressure is turned off for Harasima definition
             If (virialpress_ctype .ne. 1) Then

              ! Set up variable for i site bin
              zik1 = rzis - 0.5d0*dz
              zik2 = rzis + 0.5d0*dz
              zjk1 = rzjs - 0.5d0*dz
              zjk2 = rzjs + 0.5d0*dz

              !Calculate configurational contribution to tangential component
              ptc_2=0.5d0*x2y2tr*dudr/(dABS(dz)*area)


              !Loop through the bins
              DO ibin= 1, zden_bins

                
                !Get z-position of bin
                z = (DBLE(ibin)-0.5d0)*dz 

               
                ! Unit step function (Tangential component)
                ! 1/2 contributes to i site
                If ((z - zik1) .GT. 0.0d0) Then
                  If ((zik2-z) .GT. 0.0d0) Then

                    
                    ! Update tangential pressure from the contribution of fluid-fluid interactions           
                    virialpress_pt(2,ibin,iblock,ibox) = virialpress_pt(2,ibin,iblock,ibox) + ptc_2
 
                    
                  End If
                End If

                ! 1/2 contributes to j site
                If ((z - zjk1) .GT. 0.0d0) Then
                  If ((zjk2-z) .GT. 0.0d0) Then
                    
                    ! Update tangential pressure from the contribution of fluid-fluid interactions             
                    virialpress_pt(2,ibin,iblock,ibox) = virialpress_pt(2,ibin,iblock,ibox) + ptc_2
                   
                  End If
                End If
       
              !End bin cycle
              ENDDO

             ! End definition of contour
             End If
             
            
           ENDDO    !End loop over jmol sites

           ENDDO    !End loop over jmol
          ENDIF    !End IF imol < n_mol_tot
        !End loop over imol sites
        ENDDO  
     !End loop over imol
     ENDDO  



   ! End loop over boxes
   End do


  ! Average the pressure tensor for iblock
  ! Input: iblock
  Case(2)

    ! loop over boxes
    Do ibox = 1, n_box

      ! Loop over bins 
      Do ibin = 1, zden_bins

        ! Average normal pressure (fluid-fluid)
        virialpress_pn(ibin,iblock,ibox) = virialpress_pn(ibin,iblock,ibox)/virialpress_stat(ibox)

        ! Average tangential pressure 
        ! IK definition
        If (virialpress_ctype .ne. 2) Then
          virialpress_pt(1,ibin,iblock,ibox) = virialpress_pt(1,ibin,iblock,ibox)/virialpress_stat(ibox)
        ENDIF

        ! Harasima definition 
        If (virialpress_ctype .ne. 1) Then
          virialpress_pt(2,ibin,iblock,ibox) = virialpress_pt(2,ibin,iblock,ibox)/virialpress_stat(ibox) 
        End If
        
     
      ! End loop over bins
      End Do
    End Do
    
  ! See sample.f90 subroutine for the output of pressure tensor

  ! Error message
  Case default
    Write(*,*) 'PRESS_VIRIAL_SLIT: INVALID CASE FLAG IN PRESS_VIRIAL_SLIT.F90'
    STOP
          
      
End Select

   RETURN

END SUBROUTINE press_virial


















