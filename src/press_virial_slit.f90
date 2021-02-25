! ==========================================================
! This subroutine is used to calculate
! the pressure tensor for slit geometry from virial route
! using Irving-Kirkwood and/or Harasima definition of contour
! Adapted from Deepti's subroutine
! Created on 8-8-2017 by Kaihang Shi
! Last modified on 6-10-2020: Now for 'HARD_SLIT_FINITEX' field
!   type, the averaging region can be set freely and no extension
!   region is used (partial contribution is now fully accounted for 
!   within the cutoff raidus)
! ==========================================================



SUBROUTINE press_virial_slit(selector,iblock,istep)

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
DOUBLE PRECISION :: pnw_1, pnw_2, pnc_1, ptc_1, ptc_2, dudr, ptc_3, ptc_4
DOUBLE PRECISION :: rijs, rijstop, rijsbottom, xijs, yijs, zijs, zik1, zik2, zjk1, zjk2, zmk1, zmk2
DOUBLE PRECISION :: zijstop, zijsbottom, xijssq, yijssq, zijssq
Double precision :: posz, zii, zjj, frac_x, alp, x_alp, alp_ikvr1, alp_ikvr2
Logical :: INAVG


! Variables for the IK-VR1 contour
Double Precision :: poszA, poszB, rzAs, rzBs


Select Case (selector)
  
  ! Perform one-time calculations of pressure tensor from virial route
  ! Input: iblock
  Case(1)
    
   ! Loop over boxes
   Do ibox = 1, n_box

     ! Check if the calculation should be attempted
     IF((virialpress_slit_freq .EQ. 0) .OR. (MOD(istep,virialpress_slit_freq) .NE. 0)) CYCLE

     ! Update counter (double precision)
     virialpress_slit_stat(ibox) = virialpress_slit_stat(ibox) + 1.0d0

     !Loop over molecules
     DO imol=1,n_mol_tot(ibox)

        !Get the molecule type
        itype = mol_type(imol,ibox)

        ! Check which type of contributions should be accounted to speed up simulation
        ! Fluid-wall
        If (calc_type .EQ. 1) Then
          If (initstyle(itype,ibox) .NE. 'coords') CYCLE
        ! Fluid-fluid
        Else If (calc_type .EQ. 2) Then
          If (initstyle(itype,ibox) .EQ. 'coords') CYCLE
        ! Both
        Else If (calc_type .EQ. 3) Then
          Continue
        Endif

        !Loop over molecule sites
        DO isite=1,n_sites(itype)

          ! Get isite type 
          isitetype = site_type(isite,itype)

          !Get x,y,z position of site
          rxis = rx_s(isite,imol,ibox)
          ryis = ry_s(isite,imol,ibox)
          rzis = rz_s(isite,imol,ibox)

          ! Check if within averaging region
          If (lfield .and. (field_type .eq. STEELE_SLIT_FINITEX)) Then

            ! Include the extended region like Yun did (4sigma_ff)
            if((rxis .lt. (steele_avgx(1)-13.62d0)) .or. ((rxis .gt. steele_avgx(2)+13.62d0))) CYCLE

            if((rxis .lt. steele_avgx(1)) .or. (rxis .gt. steele_avgx(2))) goto 100

          ! Modified on 6-10-2020
          ! Now the averaing region can be set freely and no extension region is used
          Else if (lfield .and. (field_type .eq. HARD_SLIT_FINITEX)) Then

            !if((rxis .lt. (steele_avgx(1)-13.62d0)) .or. ((rxis .gt. steele_avgx(2)+13.62d0))) CYCLE
            Continue

          End If

          
          !If Steele slit pore option is on
          IF (lfield .and. ((field_type .eq. STEELE_SLIT_PORE) .or. &
              & (field_type .eq. STEELE_SLIT_FINITEX))) THEN

            !Get x,y,z position of molecule
            rxi = rx(imol,ibox)
            ryi = ry(imol,ibox)
            rzi = rz(imol,ibox)

            ! Get local, sigma and epsilon values
            sig_s = steele_sigmasf(isitetype,itype)
            eps_s = steele_epsilonsf(isitetype,itype)

            !f-w configurational contribution (For Steele potential, this will only
            !contribute to the normal pressure tensor)

            !Calculate distance of molecule COM from top wall
            zijtop = DBLE(steele_position(2) - rzi)

            !Calculate distance of molecule COM from bottom wall
            zijbottom = DBLE(rzi - steele_position(1))

            !Calculate distance of site from top wall
            zijstop = DBLE(steele_position(2) - rzis)

            !Calculate distance of site from bottom wall
            zijsbottom = DBLE(rzis - steele_position(1))

            !Calculate parameters (for top wall)
            zij = zijstop              !site z distance to the top surface
            rijs = dABS(zij)           !site distance
            z2tr = zij*zij/rijs        !Based on site-site distance

            !Calculate force (based on sites)
            sr2 = (sig_s/rijs)**2     
            sr4 = sr2*sr2
            sr10 = sr4*sr4*sr2
            phi = two_Pi*eps_s*steele_rhos*sig_s**2*steele_delta
            a = -4.0d0/rijs*sr10
            b = 4.0/rijs*sr4
            c = sig_s**4.0d0/(steele_delta*(rijs+0.61d0*steele_delta)**4.0d0)
            dudr = phi*(a+b+c)
            
            !Calculate configurational contribution to pressure
            pnw_1 = z2tr*dudr/dABS(zij)
           
            !Loop over the bins
            DO ibin=1, virialpress_slit_bins

               !Get z-position of bin
               z = (DBLE(ibin)-0.5d0)*virialpress_slit_dz + steele_position(1)

               !Test if the joining line crosses the whole z-plane
               ! two Heaviside step function
               IF ((z-rzi) .GT. 0.0d0) THEN
                IF ((steele_position(2)-z) .GT. 0.0d0) THEN
               
                !If yes, then add to pressure tensor 
                ! Contribution due to fluid-wall interactions
                virialpress_slit_pn(1,ibin,iblock,ibox) = virialpress_slit_pn(1,ibin,iblock,ibox) + pnw_1
       

                ENDIF
               ! End Heviside step function
               ENDIF   
            ! End loop over bins
            ENDDO  

            !Fluid-wall interaction with the bottom wall
            !Calculate local parameters (for bottom wall)
            zij = zijsbottom
            rijs = dABS(zij)
            z2tr = zij*zij/rijs

            !Calculate force
            sr2 = (sig_s/rijs)**2
            sr4 = sr2*sr2
            sr10 = sr4*sr4*sr2
            phi = two_Pi*eps_s*steele_rhos*sig_s**2*steele_delta
            a = -4.0d0/rijs*sr10
            b = 4.0/rijs*sr4
            c = sig_s**4.0d0/(steele_delta*(rijs+0.61*steele_delta)**4.0d0)
            dudr = phi*(a+b+c)
            
            !Calculate configurational contribution to pressure
            pnw_2 = z2tr*dudr/dABS(zij)

            !Loop over the bins
            DO ibin=1, virialpress_slit_bins

               !Get z-position of bin
               z = (DBLE(ibin)-0.5d0)*virialpress_slit_dz +steele_position(1)

               !Test if the joining line crosses the whole z-plane
               ! two Heaviside step function
               IF ((rzi-z) .GT. 0.0d0) THEN
                IF ((z-steele_position(1)) .GT. 0.0d0) THEN
               
                !If yes, then add to pressure tensor 
                ! Contribution due to fluid-wall interactions (1)
                virialpress_slit_pn(1,ibin,iblock,ibox) = virialpress_slit_pn(1,ibin,iblock,ibox) + pnw_2
       
                ENDIF
               ! End Heviside step function
               ENDIF   
            ! End loop over bins
            ENDDO  

          ! End Steele potential
          ENDIF   

100       Continue 

          !Fluid-fluid interaction, avoiding double count of i,j and j,i so just
          !loop by using i<j     
          IF (imol .LT. n_mol_tot(ibox)) THEN

            DO jmol = imol+1, n_mol_tot(ibox)

            ! Assume both molecules are not within the averaging region initially
            !INAVG = .false. 

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

              ! Check if ij pair contributes 
              ! Modified on 6-11-2020 for 'HARD_SLIT_FINITEX'. Now the averaging region can be set freely
              If (lfield .and. ((field_type .eq. STEELE_SLIT_FINITEX) .OR. (field_type .eq. HARD_SLIT_FINITEX))) Then
                ! Include the extended region like Yun did (3sigma_ff)
                !if((rxjs .lt. (steele_avgx(1)-13.62d0)) .or. ((rxjs .gt. steele_avgx(2)+13.62d0))) CYCLE

                if((rxis .lt. steele_avgx(1)) .and. (rxjs .gt. steele_avgx(2))) CYCLE
                if((rxjs .lt. steele_avgx(1)) .and. (rxis .gt. steele_avgx(2))) CYCLE

                if((rxis .lt. steele_avgx(1)) .and. (rxjs .lt. steele_avgx(1))) CYCLE
                if((rxis .gt. steele_avgx(2)) .and. (rxjs .gt. steele_avgx(2))) CYCLE
                if(rxis .lt. (steele_avgx(1) - r_cut)) CYCLE
                if(rxjs .lt. (steele_avgx(1) - r_cut)) CYCLE
                if(rxis .gt. (steele_avgx(2) + r_cut)) CYCLE
                if(rxjs .gt. (steele_avgx(2) + r_cut)) CYCLE
                !if((rxjs .lt. steele_avgx(1)) .and. (rxis .gt. steele_avgx(2))) CYCLE

            
              End If


              !Calculate vector between sites
              xijs = rxjs - rxis
              yijs = ryjs - ryis
              zijs = rzjs - rzis

              !Apply periodic boundaries
              If (lfield .and. ((field_type .eq. STEELE_SLIT_FINITEX) .OR. (field_type .eq. HARD_SLIT_FINITEX))) Then
                ! Only apply PBC in y-direction
                yijs = yijs - dNINT(yijs/box(2,ibox))*box(2,ibox)
              Else
                xijs = xijs - dNINT(xijs/box(1,ibox))*box(1,ibox)
                yijs = yijs - dNINT(yijs/box(2,ibox))*box(2,ibox)
              End If

              ! Set up basic parameters for each contour definition
              ! Calculate the z-position for point A and A' (B) in IK-VR contour
              If ((virialpress_ctype .EQ. 5) .or. (virialpress_ctype .EQ. 6)) Then
                rzAs = (rzis + 2* rzjs)/3.0
                rzBs = (rzjs + 2* rzis)/3.0
              End if

              
              If (lfield .and. ((field_type .eq. STEELE_SLIT_FINITEX) .OR. (field_type .eq. HARD_SLIT_FINITEX))) Then

                ! Set up prefactor for Harasima and H-VR
                If((virialpress_ctype .NE. 1) .and. (virialpress_ctype .NE. 5)) Then

                  ! The averaging region is between ij pair
                  If ( ((rxis .le. steele_avgx(1)) .and. (rxjs .ge. steele_avgx(2))) .OR. &
                      & ((rxjs .lt. steele_avgx(1)) .and. (rxis .ge. steele_avgx(2))) ) Then

                    ! calculate the prefactor 
                    frac_x = dABS(l_avgx)/dABS(xijs)

                  ! if both in the averaing region
                  Else If (  (rxis .gt. steele_avgx(1)) .and. (rxis .lt. steele_avgx(2)) .and. &
                      & (rxjs .gt. steele_avgx(1)) .and. (rxjs .lt. steele_avgx(2)) ) Then

                    frac_x = 1.0d0

                  ! Only one particle in the averaging region
                  Else if ( (rxis .lt. steele_avgx(1)) .and. (rxjs .gt. steele_avgx(1)) ) Then

                    frac_x = dABS(rxjs - steele_avgx(1))/dABS(xijs)

                  Else if ( (rxjs .lt. steele_avgx(1)) .and. (rxis .gt. steele_avgx(1)) ) Then

                    frac_x = dABS(rxis - steele_avgx(1))/dABS(xijs)

                  Else if ( (rxis .lt. steele_avgx(2)) .and. (rxjs .gt. steele_avgx(2)) ) Then

                    frac_x = dABS(steele_avgx(2) - rxis)/dABS(xijs)

                  Else if ( (rxjs .lt. steele_avgx(2)) .and. (rxis .gt. steele_avgx(2)) ) Then

                    frac_x = dABS(steele_avgx(2) - rxjs)/dABS(xijs)

                  Else 

                    Write(*,*) "PRESS_VIRIAL_SLIT: THERE IS AN ERROR IN CODE: line 318"
                    STOP

                  End If 
                ! END for Harasima and H-VR
                End If 

              !            
              End If
              
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


              !!!! Definition of integral contour !!!!
              ! ------------------------------------------------------------------
              ! (1) Only Irving-Kirkwood definition or (3) IK&H or (6) ALL
              If ((virialpress_ctype .EQ. 1) .OR. (virialpress_ctype .EQ. 3) .OR. (virialpress_ctype .EQ. 6)) Then

                !Calculate configurational contribution
                ptc_1=x2y2tr*dudr/dABS(zijs)
                pnc_1=z2tr*dudr/dABS(zijs)

                !Loop through the bins
                DO ibin= 1, virialpress_slit_bins

                  If (lfield .and. ((field_type .eq. STEELE_SLIT_FINITEX) .OR. (field_type .eq. STEELE_SLIT_PORE))) Then

                    !Get z-position of bin
                    z = (DBLE(ibin)-0.5d0)*virialpress_slit_dz + steele_position(1)

                    !Unit step function
                    IF ((z - rzis)/zijs .GT. 0.0d0) THEN
                      IF ((rzjs - z)/zijs .GT.  0.0d0) THEN

                        ! For atomistic model, fluid-wall
                        If (initstyle(itype,ibox) .eq. 'coords') Then
                          ! Update tangential pressure from the contribution of fluid-wall interaction (1)
                          virialpress_slit_pt(1,1,ibin,iblock,ibox) = virialpress_slit_pt(1,1,ibin,iblock,ibox) + ptc_1  
                          ! Update normal pressure from the contribution of fluid-wall interactions (1)
                          virialpress_slit_pn(1,ibin,iblock,ibox) = virialpress_slit_pn(1,ibin,iblock,ibox) + pnc_1

                        Else
                          ! Update normal pressure from the contribution of fluid-fluid interactions (2)
                          virialpress_slit_pn(2,ibin,iblock,ibox) = virialpress_slit_pn(2,ibin,iblock,ibox) + pnc_1
                          ! Update tangential pressure from the contribution of fluid-fluid interactions (2)
                          ! In slit geometry, only fluid-fluid interactions contribute to the tangential pressure
                          virialpress_slit_pt(2,1,ibin,iblock,ibox) = virialpress_slit_pt(2,1,ibin,iblock,ibox) + ptc_1 

                        End If

                      ENDIF 
                    ! End Heavise unit step function
                    ENDIF


                  ! Modified on 6-11-2020
                  Else if (lfield .and. (field_type .eq. HARD_SLIT_FINITEX) ) Then
                    !Get z-position of bin
                    z = (DBLE(ibin)-0.5d0)*virialpress_slit_dz + wall_radius

                    ! calculate alpha value
                    alp = (z - rzis)/zijs

                    ! Exclude contribution out of the averaging region
                    If ((alp .lt. 0.0d0) .or. (alp .gt. 1.0d0)) CYCLE

                    x_alp = rxis + alp*xijs

                    If ((x_alp .ge. steele_avgx(1)) .and. (x_alp .le. steele_avgx(2))) Then

                      ! For atomistic model, fluid-wall
                      If (initstyle(itype,ibox) .eq. 'coords') Then
                        ! Update tangential pressure from the contribution of fluid-wall interaction (1)
                        virialpress_slit_pt(1,1,ibin,iblock,ibox) = virialpress_slit_pt(1,1,ibin,iblock,ibox) + ptc_1  
                        ! Update normal pressure from the contribution of fluid-wall interactions (1)
                        virialpress_slit_pn(1,ibin,iblock,ibox) = virialpress_slit_pn(1,ibin,iblock,ibox) + pnc_1

                      Else
                        ! Update normal pressure from the contribution of fluid-fluid interactions (2)
                        virialpress_slit_pn(2,ibin,iblock,ibox) = virialpress_slit_pn(2,ibin,iblock,ibox) + pnc_1
                        ! Update tangential pressure from the contribution of fluid-fluid interactions (2)
                        virialpress_slit_pt(2,1,ibin,iblock,ibox) = virialpress_slit_pt(2,1,ibin,iblock,ibox) + ptc_1 

                      End If

                    ENDIF 
                    
                  Else 
                    Write (*,*) 'PRESS_VIRIAL_SLIT: ONLY WORKS FOR STEELE AND HARD_SLIT_FINITEX'
                    STOP
                  End If
                  
                !End bin cycle
                ENDDO
              ! End IK 
              ENDIF

              ! ---------------------------------------------------
              ! Harasima definition (modifed the algorithm for averaing in the region on 5/2/2020)
              ! Normal pressure is turned off for Harasima definition
              If ((virialpress_ctype .EQ. 2) .OR. (virialpress_ctype .EQ. 3) .OR. (virialpress_ctype .EQ. 6)) Then

                ! Set up variable for i site bin
                zik1 = rzis - 0.5d0*virialpress_slit_dz
                zik2 = rzis + 0.5d0*virialpress_slit_dz
                zjk1 = rzjs - 0.5d0*virialpress_slit_dz
                zjk2 = rzjs + 0.5d0*virialpress_slit_dz

                !Calculate configurational contribution to tangential component
                ptc_2=0.5d0*x2y2tr*dudr/dABS(virialpress_slit_dz)
                ! Normal component is the same as that for IK definition
                !pnc_1=z2tr*dudr/dABS(zijs)


                !Loop through the bins
                DO ibin= 1, virialpress_slit_bins

                  If (lfield .and. ((field_type .eq. STEELE_SLIT_FINITEX) .OR. (field_type .eq. STEELE_SLIT_PORE))) Then
                    !Get z-position of bin
                    z = (DBLE(ibin)-0.5d0)*virialpress_slit_dz + steele_position(1)

                  Else if (lfield .and. (field_type .eq. HARD_SLIT_FINITEX) ) Then
                    !Get z-position of bin
                    z = (DBLE(ibin)-0.5d0)*virialpress_slit_dz + wall_radius

                  Else 
                    Write (*,*) 'PRESS_VIRIAL_SLIT: ONLY WORKS FOR STEELE AND HARD_SLIT_FINITEX'
                    STOP
                  End If
               
                  If (lfield .and. ((field_type .eq. STEELE_SLIT_FINITEX) .OR. (field_type .eq. HARD_SLIT_FINITEX)))  Then
                  
                    ! Multiply a factor (modified on 5/2/2020)
                    ! Unit step function (Tangential component)
                    ! 1/2 contributes to i site
                    If ((z - zik1) .GT. 0.0d0) Then
                      If ((zik2-z) .GT. 0.0d0) Then

                        If (initstyle(itype,ibox) .eq. 'coords') Then
                          ! Update tangential pressure from the contribution of fluid-wall interactions (1)
                          virialpress_slit_pt(1,2,ibin,iblock,ibox) = virialpress_slit_pt(1,2,ibin,iblock,ibox) + frac_x * ptc_2
                        Else 
                          ! Update tangential pressure from the contribution of fluid-fluid interactions (2)              
                          virialpress_slit_pt(2,2,ibin,iblock,ibox) = virialpress_slit_pt(2,2,ibin,iblock,ibox) + frac_x * ptc_2
                        End If
                        
                      End If
                    End If

                    ! 1/2 contributes to j site
                    If ((z - zjk1) .GT. 0.0d0) Then
                      If ((zjk2-z) .GT. 0.0d0) Then

                        If (initstyle(itype,ibox) .eq. 'coords') Then
                          ! Update tangential pressure from the contribution of fluid-wall interactions (1)
                          virialpress_slit_pt(1,2,ibin,iblock,ibox) = virialpress_slit_pt(1,2,ibin,iblock,ibox) + frac_x * ptc_2
                        Else 
                          ! Update tangential pressure from the contribution of fluid-fluid interactions (2)              
                          virialpress_slit_pt(2,2,ibin,iblock,ibox) = virialpress_slit_pt(2,2,ibin,iblock,ibox) + frac_x * ptc_2
                        End If
                        
                      End If
                    End If

                  ! Not steele_slit_finitex
                  Else

                    ! Unit step function (Tangential component)
                    ! 1/2 contributes to i site
                    If ((z - zik1) .gt. 0.0d0) Then
                      If ((zik2-z) .gt. 0.0d0) Then
                        ! Update tangential pressure from the contribution of fluid-fluid interactions (2)
                        ! In slit geometry, only fluid-fluid interactions contribute to the tangential pressure
                        virialpress_slit_pt(2,2,ibin,iblock,ibox) = virialpress_slit_pt(2,2,ibin,iblock,ibox) + ptc_2
                      End If
                    End If
                    ! 1/2 contributes to j site
                    If ((z - zjk1) .gt. 0.0d0) Then
                      If ((zjk2-z) .gt. 0.0d0) Then
                        ! Update tangential pressure from the contribution of fluid-fluid interactions (2)
                        ! In slit geometry, only fluid-fluid interactions contribute to the tangential pressure
                        virialpress_slit_pt(2,2,ibin,iblock,ibox) = virialpress_slit_pt(2,2,ibin,iblock,ibox) + ptc_2
                      End If
                    End If

                  End If

                !End bin cycle
                ENDDO

              ! End definition of Harasima contour
              End If

              ! ---------------------------------------------------------
              ! Added on May 29, 2019, last modified on June 12, 2020
              ! Variation of the Harasima type (H-VR) 
              ! Normal pressure is turned off 
              If ((virialpress_ctype .EQ. 4) .OR. (virialpress_ctype .EQ. 6))  Then

                ! Set up variable for bin
                zmk1 = 0.5d0*(rzis+rzjs) - 0.5d0*virialpress_slit_dz
                zmk2 = 0.5d0*(rzis+rzjs) + 0.5d0*virialpress_slit_dz

                !Calculate configurational contribution to tangential component
                ptc_3=x2y2tr*dudr/dABS(virialpress_slit_dz)

                !Loop through the bins
                DO ibin= 1, virialpress_slit_bins

                  ! Calculate z-position of each bin
                  If (lfield .and. ((field_type .eq. STEELE_SLIT_FINITEX) .OR. (field_type .eq. STEELE_SLIT_PORE))) Then
                    !Get z-position of bin
                    z = (DBLE(ibin)-0.5d0)*virialpress_slit_dz + steele_position(1)
                  Else if (lfield .and. (field_type .eq. HARD_SLIT_FINITEX) ) Then
                    !Get z-position of bin
                    z = (DBLE(ibin)-0.5d0)*virialpress_slit_dz + wall_radius
                  Else 
                    Write (*,*) 'PRESS_VIRIAL_SLIT: ONLY WORKS FOR STEELE AND HARD_SLIT_FINITEX'
                    STOP
                  End If
               
                  If (lfield .and. ((field_type .eq. STEELE_SLIT_FINITEX) .OR. (field_type .eq. HARD_SLIT_FINITEX)))  Then
      
                    ! unit function
                    If ((z - zmk1) .GT. 0.0d0) Then
                      If ((zmk2-z) .GT. 0.0d0) Then

                        If (initstyle(itype,ibox) .eq. 'coords') Then
                          ! Update tangential pressure from the contribution of fluid-wall interactions (1)
                          virialpress_slit_pt(1,3,ibin,iblock,ibox) = virialpress_slit_pt(1,3,ibin,iblock,ibox) + frac_x * ptc_3

                        Else 
                          ! Update tangential pressure from the contribution of fluid-fluid interactions (2)              
                          virialpress_slit_pt(2,3,ibin,iblock,ibox) = virialpress_slit_pt(2,3,ibin,iblock,ibox) + frac_x * ptc_3
                        End If
                        
                      End If
                    End If

                  Else 
                    Write(*,*) 'H-VR only works for STEELE_SLIT_FINITEX and HARD_SLIT_FINITEX now.'
                    STOP
                    
                  ! End lfield
                  End If

                !End bin cycle
                ENDDO

              ! End definition of H-VR contour
              End If

              ! ------------------------------------------------------------------------
              ! Added on May 30, 2019, last modified on June 12 ,2020
              ! New contour, a variation of the IK type (IK-VR)
              If ((virialpress_ctype .EQ. 5) .OR. (virialpress_ctype .EQ. 6)) Then

                !Calculate configurational contribution
                ptc_4 = 1.5d0*x2y2tr*dudr/dABS(zijs)


                !Loop through the bins
                DO ibin= 1, virialpress_slit_bins

                  If (lfield .and. ((field_type .eq. STEELE_SLIT_FINITEX) .OR. (field_type .eq. STEELE_SLIT_PORE))) Then
                    !Get z-position of bin
                    z = (DBLE(ibin)-0.5d0)*virialpress_slit_dz + steele_position(1)

                    !Unit step function
                    IF ((z - rzAs)/zijs .GT. 0.0d0) THEN
                      IF ((rzjs - z)/zijs .GT.  0.0d0) THEN

                        ! For atomistic model, fluid-wall
                        If (initstyle(itype,ibox) .eq. 'coords') Then
                          ! Update tangential pressure from the contribution of fluid-wall interaction (1)
                          virialpress_slit_pt(1,4,ibin,iblock,ibox) = virialpress_slit_pt(1,4,ibin,iblock,ibox) + ptc_4
                        Else
                          ! Update tangential pressure from the contribution of fluid-fluid interactions (2)
                          virialpress_slit_pt(2,4,ibin,iblock,ibox) = virialpress_slit_pt(2,4,ibin,iblock,ibox) + ptc_4
                        End If
                      ENDIF 
                    ! End Heavise unit step function
                    ENDIF

                    !Unit step function
                    IF ((rzBs - z)/zijs .GT. 0.0d0) THEN
                      IF ((z - rzis)/zijs .GT.  0.0d0) THEN

                        ! For atomistic model, fluid-wall
                        If (initstyle(itype,ibox) .eq. 'coords') Then
                          ! Update tangential pressure from the contribution of fluid-wall interaction (1)
                          virialpress_slit_pt(1,4,ibin,iblock,ibox) = virialpress_slit_pt(1,4,ibin,iblock,ibox) + ptc_4
                        Else
                          ! Update tangential pressure from the contribution of fluid-fluid interactions (2)
                          virialpress_slit_pt(2,4,ibin,iblock,ibox) = virialpress_slit_pt(2,4,ibin,iblock,ibox) + ptc_4
                        End If
                      ENDIF 
                    ! End Heavise unit step function
                    ENDIF

                  ! Modified 6-12-2020
                  Else if (lfield .and. (field_type .eq. HARD_SLIT_FINITEX) ) Then
                    !Get z-position of bin
                    z = (DBLE(ibin)-0.5d0)*virialpress_slit_dz + wall_radius

                    ! calculate alpha value for i->A>j contour
                    alp_ikvr1 = (z - rzAs)/(rzjs - rzAs)
                    x_alp = rxis + alp_ikvr1*xijs

                    If ( (alp_ikvr1 .ge. 0.0d0) .and. (alp_ikvr1 .le. 1.0d0) .and. (x_alp .ge. steele_avgx(1)) &
                        & .and. (x_alp .le. steele_avgx(2)) ) Then

                      ! For atomistic model, fluid-wall
                      If (initstyle(itype,ibox) .eq. 'coords') Then
                        ! Update tangential pressure from the contribution of fluid-wall interaction (1)
                        virialpress_slit_pt(1,4,ibin,iblock,ibox) = virialpress_slit_pt(1,4,ibin,iblock,ibox) + ptc_4
                      Else
                        ! Update tangential pressure from the contribution of fluid-fluid interactions (2)
                        virialpress_slit_pt(2,4,ibin,iblock,ibox) = virialpress_slit_pt(2,4,ibin,iblock,ibox) + ptc_4
                      End If

                    ! End Heavise unit step function
                    ENDIF

                    ! calculate alpha value for j->B>i contour
                    alp_ikvr2 = (z - rzBs)/(rzis - rzBs)
                    x_alp = rxjs - alp_ikvr2*xijs

                    If ( (alp_ikvr2 .ge. 0.0d0) .and. (alp_ikvr2 .le. 1.0d0) .and. (x_alp .ge. steele_avgx(1)) &
                        & .and. (x_alp .le. steele_avgx(2)) ) Then

                        ! For atomistic model, fluid-wall
                        If (initstyle(itype,ibox) .eq. 'coords') Then
                          ! Update tangential pressure from the contribution of fluid-wall interaction (1)
                          virialpress_slit_pt(1,4,ibin,iblock,ibox) = virialpress_slit_pt(1,4,ibin,iblock,ibox) + ptc_4
                        Else
                          ! Update tangential pressure from the contribution of fluid-fluid interactions (2)
                          virialpress_slit_pt(2,4,ibin,iblock,ibox) = virialpress_slit_pt(2,4,ibin,iblock,ibox) + ptc_4
                        End If
                    End If 

                  Else 
                    Write (*,*) 'PRESS_VIRIAL_SLIT: ONLY WORKS FOR STEELE AND HARD_SLIT_FINITEX'
                    STOP
                  End If
              
                !End bin cycle
                ENDDO
              ! End IK 
              ENDIF
             
            
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
      Do ibin = 1, virialpress_slit_bins

        ! Average normal pressure (fluid-wall)
        virialpress_slit_pn(1,ibin,iblock,ibox) = virialpress_slit_pn(1,ibin,iblock,ibox)/virialpress_slit_stat(ibox)

        ! Average normal pressure (fluid-fluid)
        virialpress_slit_pn(2,ibin,iblock,ibox) = virialpress_slit_pn(2,ibin,iblock,ibox)/virialpress_slit_stat(ibox)

        ! Average tangential pressure 
        ! IK definition
        If ((virialpress_ctype .EQ. 1) .OR. (virialpress_ctype .EQ. 3) .OR. (virialpress_ctype .EQ. 6)) Then
          virialpress_slit_pt(1,1,ibin,iblock,ibox) = virialpress_slit_pt(1,1,ibin,iblock,ibox)/virialpress_slit_stat(ibox)
          virialpress_slit_pt(2,1,ibin,iblock,ibox) = virialpress_slit_pt(2,1,ibin,iblock,ibox)/virialpress_slit_stat(ibox)
        ENDIF

        ! Harasima definition 
        If ((virialpress_ctype .EQ. 2) .OR. (virialpress_ctype .EQ. 3) .OR. (virialpress_ctype .EQ. 6)) Then
          virialpress_slit_pt(1,2,ibin,iblock,ibox) = virialpress_slit_pt(1,2,ibin,iblock,ibox)/virialpress_slit_stat(ibox)
          virialpress_slit_pt(2,2,ibin,iblock,ibox) = virialpress_slit_pt(2,2,ibin,iblock,ibox)/virialpress_slit_stat(ibox) 
        End If

        ! Added on May 29, 2019
        ! H-VR definition 
        If ((virialpress_ctype .EQ. 4) .OR. (virialpress_ctype .EQ. 6)) Then
          virialpress_slit_pt(1,3,ibin,iblock,ibox) = virialpress_slit_pt(1,3,ibin,iblock,ibox)/virialpress_slit_stat(ibox)
          virialpress_slit_pt(2,3,ibin,iblock,ibox) = virialpress_slit_pt(2,3,ibin,iblock,ibox)/virialpress_slit_stat(ibox) 
        End If

        ! Added on May 30, 2019
        ! IK-VR definition 
        If ((virialpress_ctype .EQ. 5) .OR. (virialpress_ctype .EQ. 6)) Then
          virialpress_slit_pt(1,4,ibin,iblock,ibox) = virialpress_slit_pt(1,4,ibin,iblock,ibox)/virialpress_slit_stat(ibox)
          virialpress_slit_pt(2,4,ibin,iblock,ibox) = virialpress_slit_pt(2,4,ibin,iblock,ibox)/virialpress_slit_stat(ibox) 
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

END SUBROUTINE press_virial_slit


















