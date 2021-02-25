! ==================================================================
! This subroutine is used to set basic parameters for the cell list
!   and assign neighbors to each cell
! Created on 7-12-2019 by Kaihang Shi 
! ==================================================================


 	    Subroutine set_clist

 	    Use global

 	    IMPLICIT NONE


   	  ! Local
   	  Integer :: ibox, ibinx, ibiny, ibinz, icel, jcel, inei
   	  Double Precision, Dimension(:), Allocatable :: rxc, ryc, rzc
   	  Double Precision :: cutsq, xijc, yijc, zijc, rijcsq
   	  Integer :: ierr


   	  ! Assume only one box
   	  ibox = 1

   	  ! Initialize local variables
   	  icel = 0


   	  ! Number of cells in each dimension
   	  clist_nx = INT(box(1,ibox)/r_cut)
   	  clist_ny = INT(box(2,ibox)/r_cut)
      clist_nz = INT(box(3,ibox)/r_cut)

   	  ! Cell list must have at least 3 cells in each dimension
   	  If ((clist_nx .lt. 3) .or. (clist_ny .lt. 3) .or. (clist_nz .lt. 3)) Then
   	  	Write(*,*) 'System size is too small for cell list. Minimum of 3 cells in each dimension!'
   	  	STOP
   	  End If

   	  ! Total number of cells
   	  clist_ncel = clist_nx*clist_ny*clist_nz

      ! In current version, maximum of 125000 cells are allowed
      If (clist_ncel .GT. 125000) Then
        Write(*,*) 'Number of cell in cell list exceed maximum 125000'
        STOP
      End If

   	  ! Allocate position array for cells
   	  Allocate(rxc(clist_ncel),STAT=ierr)
   	  Allocate(ryc(clist_ncel),STAT=ierr)
      Allocate(rzc(clist_ncel),STAT=ierr)

   	  ! Lateral size for each cell
   	  clist_dx = box(1,ibox)/DBLE(clist_nx)
   	  clist_dy = box(2,ibox)/DBLE(clist_ny)
      clist_dz = box(3,ibox)/DBLE(clist_nz)

   	  ! Calculate cutoff radius for assigning neighboring cells
   	  cutsq = clist_dx**2 + clist_dy**2 + clist_dz**2

   	  ! Get position of each cell
  	  Do ibinx = 1, clist_nx
    		Do ibiny = 1, clist_ny
          Do ibinz = 1, clist_nz

            icel = icel + 1
            ! Recover the x-position of the cell
            rxc(icel) = (DBLE(ibinx)-0.5d0)*clist_dx
            ! Recover the y-position
            ryc(icel) = (DBLE(ibiny)-0.5d0)*clist_dy
            ! z-position
            rzc(icel) = (DBLE(ibinz)-0.5d0)*clist_dz

            ! Setup cell locator
            clist_loca(ibinx,ibiny,ibinz) = icel

          EndDo
        EndDo
      EndDo

      ! Loop over cells and assign its neighbors
      Do icel = 1, clist_ncel

      	! Initalize variable
      	inei = 0

      	Do jcel = 1, clist_ncel

      		! Calculate distance between cells
      		xijc = rxc(icel) - rxc(jcel)
      		yijc = ryc(icel) - ryc(jcel)
          zijc = rzc(icel) - rzc(jcel)

      		! Apply minimum image conversion 
          xijc = xijc - dNINT(xijc/box(1,ibox))*box(1,ibox)
          yijc = yijc - dNINT(yijc/box(2,ibox))*box(2,ibox)
          zijc = zijc - dNINT(zijc/box(3,ibox))*box(3,ibox)

          ! square distance
          rijcsq = xijc**2 + yijc**2 + zijc**2

          ! Assign neighbors to icell
          If (rijcsq .le. (cutsq+0.1d0)) Then

            inei = inei + 1
            clist_neigh(icel,inei) = jcel

          End If
               
      	End Do  

      	! Check neighbours
      	If (inei .ne. 27) Then
           	Write(*,*) 'SET_CLIST: neighbors are not 27 for 3D!'
           	Write(*,*) inei, icel
            STOP
        End If 	
      ! Finish assigning neighbors
      End Do

  	  Return
   	  
   	  End Subroutine 







