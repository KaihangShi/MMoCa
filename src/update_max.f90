! ==========================================================
! This subroutine is used to update max values (max_trans,
! max_rotat etc)
! Created on 1-1-2017 by Kaihang Shi
! ==========================================================


      Subroutine update_max

      Use global

      IMPLICIT NONE

      ! Local
      Integer :: ibox, itype
      Double Precision :: accratio


      ! Loop over boxes
      Do ibox = 1, n_box

        ! Loop over molecule types
        Do itype = 1, n_mol_types

            ! Update max displacement
            If (trans_stat(1,itype,ibox) .gt. 0) Then
                ! Calculate acceptance ratio
                accratio = DBLE(trans_stat(2,itype,ibox))/DBLE(trans_stat(1,itype,ibox))

                ! Adjust maximum displacement (Allen & Tildesley)
                If (accratio .ge. 0.5d0) Then
                    ! Increase max_trans
                    IF(1.05d0*max_trans(itype,ibox) .LT. MINVAL(box(:,ibox))) &
                    & max_trans(itype,ibox) = 1.05d0*max_trans(itype,ibox)

                Else
                    ! Decrease max_trans
                    max_trans(itype,ibox) = 0.95d0*max_trans(itype,ibox)

                End If

            End If

            ! Update the maximum rotation
            IF(rotat_stat(1,itype,ibox) .GT. 0) THEN

                ! Calculate the acceptance ratio
                accratio = DBLE(rotat_stat(2,itype,ibox))/DBLE(rotat_stat(1,itype,ibox))

                ! Adjust the maximum rotation
                IF(accratio .GE. 0.5d0) THEN
                  IF(1.05d0*max_rotat(itype,ibox) .LT. 10.0) &
                  & max_rotat(itype,ibox) = 1.05d0*max_rotat(itype,ibox)
                ELSE
                  IF(0.95d0*max_rotat(itype,ibox) .GT. 1.0d-5) &
                  & max_rotat(itype,ibox) = 0.95d0*max_rotat(itype,ibox)
                ENDIF
         
            ENDIF

            
        End Do

        ! Update the maximum volume change
        If (vol_stat(1,ibox) .gt. 0) Then

            ! Calculate the acceptance ratio
            accratio = DBLE(vol_stat(2,ibox))/DBLE(vol_stat(1,ibox))

            ! Adjust the maximum volume change
            If (accratio .ge. 0.5d0) Then
                If (max_vol(ibox) .le. 0.5d0) max_vol(ibox) = 1.05d0*max_vol(ibox)
            ELSE
                max_vol(ibox) = 0.95d0*max_vol(ibox)
            End If
            
        End If


100     Continue

      ! End loop over boxes
      End Do

      
        
      
      Return
      
      End Subroutine 