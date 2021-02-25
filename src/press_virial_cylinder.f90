! ==========================================================
! This subroutine is used to calculate
! the pressure tensor for cylindrical geometry from virial route
! using Irving-Kirkwood and/or Harasima definition of contour
! ATTENTION: No periodic boundary conditions are applied. We 
!           assume there is no PBC in x- and y-directions.
! Created on Feb 27 2019 by Kaihang Shi
! Last update: March 2019
! ==========================================================



SUBROUTINE press_virial_cylinder(selector,iblock,istep)

USE global
IMPLICIT NONE

! Passed
Integer :: selector, iblock, istep

!Local variables
INTEGER :: imol, jmol, ibin, ibox
INTEGER :: isite, jsite, itype, jtype, isitetype, jsitetype

DOUBLE PRECISION :: rxi, ryi, rzi, xij, yij, zij, xijsq, yijsq, zijsq, rij, rxj, ryj, rzj
DOUBLE PRECISION :: rxis, ryis, rzis, rxjs, ryjs, rzjs
DOUBLE PRECISION :: sr3, sr6, sr12, phitz
DOUBLE PRECISION :: pnrc, ptzc, dudr, pttci, pttcj
DOUBLE PRECISION :: rijs, xijs, yijs, zijs
DOUBLE PRECISION :: xijssq, yijssq, zijssq
Double Precision :: clri, clrj, clrij, posr, alpk, cos_thetaij




Select Case (selector)
  
  ! Perform one-time calculations of pressure tensor from virial route
  ! Input: iblock
  Case(1)
    
   ! Loop over boxes
   Do ibox = 1, n_box

     ! Check if the calculation should be attempted
     IF((virialpress_cylin_freq .EQ. 0) .OR. (MOD(istep,virialpress_cylin_freq) .NE. 0)) CYCLE

     ! Update counter (double precision)
     virialpress_cylin_stat(ibox) = virialpress_cylin_stat(ibox) + 1.0d0

     !Loop over molecules
     DO imol=1,n_mol_tot(ibox)

        !Get the molecule type
        itype = mol_type(imol,ibox)

        !Loop over molecule sites
        DO isite=1,n_sites(itype)

          ! Get isite type 
          isitetype = site_type(isite,itype)

          !Get x,y,z position of site and move x,y-position to the image box with
          ! (0,0,0) as the center of its bottom surface
          rxis = rx_s(isite,imol,ibox) - 0.5d0*box(1,ibox)
          ryis = ry_s(isite,imol,ibox) - 0.5d0*box(2,ibox)
          rzis = rz_s(isite,imol,ibox)

          ! Calculate R-distance of site i in cylindrical coordiantes
          clri = dSQRT(rxis**2 + ryis**2)

          ! Determine if carrying on (rden_lim = rden_cut + r_cut)
          If (clri .GT. rden_lim) CYCLE


          !Fluid-fluid interaction, avoiding double count of i,j and j,i so just
          !loop by using i<j     
          IF (imol .LT. n_mol_tot(ibox)) THEN
            DO jmol = imol+1, n_mol_tot(ibox)

              !Get molecule type
              jtype = mol_type(jmol,ibox)

              !Loop over molecule sites
              DO jsite=1,n_sites(jtype)

                ! Get jsite type
                jsitetype = site_type(jsite,jtype)

                !Get x,y,z position of molecule and move x,y-position to the image box with
                ! (0,0,0) as the center of its bottom surface
                rxjs = rx_s(jsite,jmol,ibox) - 0.5d0*box(1,ibox)
                ryjs = ry_s(jsite,jmol,ibox) - 0.5d0*box(2,ibox)
                rzjs = rz_s(jsite,jmol,ibox)

                ! Calculate R-distance of site j 
                clrj = dSQRT(rxjs**2+ryjs**2)

                ! Harasima
                If (virialpress_ctype .eq. 2) Then
                  ! Determine if carrying on
                  If ((clri .GT. rden_cut) .AND. (clrj .GT. rden_cut)) CYCLE
                ENDIF

                !Calculate vector between sites
                xijs = rxjs - rxis
                yijs = ryjs - ryis
                zijs = rzjs - rzis
                
                ! Apply minimum image convention
                ! Only apply to z-direction
                zijs = zijs - dNINT(zijs/box(3,ibox))*box(3,ibox)

                !Square the values
                xijssq=xijs*xijs
                yijssq=yijs*yijs
                zijssq=zijs*zijs

                !Calculate distance between i and j sites
                rijs=dSQRT(xijssq+yijssq+zijssq)

                ! Check with force cutoff
                If (rijs .GT. r_cut) CYCLE

                ! Difference
                clrij = clrj - clri

                ! calculate COS(theta)
                cos_thetaij = (rxis*rxjs+ryis*ryjs)/(clri*clrj)


                !Calculate 12-6LJ force (based on sites)
                sr3=(sigma(isitetype,itype,jsitetype,jtype)/rijs)**3
                sr6=sr3*sr3
                sr12=sr6*sr6
                phitz = 24.0*epsilon(isitetype,itype,jsitetype,jtype)/rijs
                dudr = phitz*(sr6 - 2.0d0*sr12)

                ! Definition of integral contour
                ! Harasima contour
                If (virialpress_ctype .eq. 2) Then

                  ! Calculate normal pressure in radial direction
                  pnrc = dudr*dABS(clrij)*(1.0d0+cos_thetaij)/(rijs*box(3,ibox))

                  ! Calculate tangential pressure in theta-direction
                  pttci = 0.5d0*clri*((rxjs/clrj - rxis/clri)*xijs + (ryjs/clrj - ryis/clri)*yijs)*dudr/(rijs*box(3,ibox))
                  pttcj = pttci*clrj/clri
                  
                  ! Calculate tangential pressure in z-direction
                  ptzc = 0.5d0*dudr*zijs**2/(rijs*box(3,ibox))


                  !Loop through the bins
                  DO ibin= 1, rden_bins

                    ! Get r-distance of ibin in cylindrical system
                    posr = (DBLE(ibin)-0.5d0)*rden_dr

                    ! Calculate alpha_k
                    alpk = (posr - clri)/clrij
                    
                    ! P_RR
                    !Unit step function
                    IF ((alpk .GT. 0.0d0) .and. (alpk .LE.  1.0d0)) THEN

                        ! Update normal pressure from the contribution of fluid-fluid interactions (2), harasima(2)
                        virialpress_cylin_pnr(2,2,ibin,iblock,ibox) = virialpress_cylin_pnr(2,2,ibin,iblock,ibox) + pnrc

                    ! End Heavisde step function                         
                    ENDIF

                    ! P_ThetaTheta
                    ! unit step function
                    ! Half contribute to particle i
                    If ((posr-clri+delrr) .gt. 0.0d0) Then
                      If ((clri+delrr-posr) .gt. 0.0d0) Then
                        
                        virialpress_cylin_ptt(2,2,ibin,iblock,ibox) = virialpress_cylin_ptt(2,2,ibin,iblock,ibox) + pttci

                      End If
                    End If

                    ! Half contribute to particle j
                    If ((posr-clrj+delrr) .gt. 0.0d0) Then
                      If ((clrj+delrr-posr) .gt. 0.0d0) Then
                        
                        virialpress_cylin_ptt(2,2,ibin,iblock,ibox) = virialpress_cylin_ptt(2,2,ibin,iblock,ibox) + pttcj

                      End If
                    End If


                    ! P_zz
                    ! unit step function
                    ! Half contribute to particle i
                    If ((posr-clri+delrr) .gt. 0.0d0) Then
                      If ((clri+delrr-posr) .gt. 0.0d0) Then
                        
                        virialpress_cylin_ptz(2,2,ibin,iblock,ibox) = virialpress_cylin_ptz(2,2,ibin,iblock,ibox) + ptzc

                      End If
                    End If

                    ! Half contribute to particle j
                    If ((posr-clrj+delrr) .gt. 0.0d0) Then
                      If ((clrj+delrr-posr) .gt. 0.0d0) Then
                        
                        virialpress_cylin_ptz(2,2,ibin,iblock,ibox) = virialpress_cylin_ptz(2,2,ibin,iblock,ibox) + ptzc

                      End If
                    End If

                    


                  !End bin cycle
                  ENDDO
               ! End Haraisma 
               ENDIF
              !End loop over jmol sites 
              ENDDO    
            !End loop over jmol 
            ENDDO    
          !End IF imol < n_mol_tot  
          ENDIF    
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
      Do ibin = 1, rden_bins

        ! Average pressure component
        ! Harasima definition
        If (virialpress_ctype .eq. 2) Then
          virialpress_cylin_pnr(2,2,ibin,iblock,ibox) = virialpress_cylin_pnr(2,2,ibin,iblock,ibox)/virialpress_cylin_stat(ibox)
          virialpress_cylin_ptt(2,2,ibin,iblock,ibox) = virialpress_cylin_ptt(2,2,ibin,iblock,ibox)/virialpress_cylin_stat(ibox)
          virialpress_cylin_ptz(2,2,ibin,iblock,ibox) = virialpress_cylin_ptz(2,2,ibin,iblock,ibox)/virialpress_cylin_stat(ibox)
         

        ENDIF
        
      ! End loop over bins
      End Do
      
      
    ! End loop over box
    End Do
    
  ! See sample.f90 subroutine for the output of pressure tensor

  ! Error message
  Case default
    Write(*,*) 'PRESS_VIRIAL_SLIT: INVALID CASE FLAG IN PRESS_VIRIAL_SLIT.F90'
    STOP
          
      
End Select

   RETURN

END SUBROUTINE press_virial_cylinder







