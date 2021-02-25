! ==========================================================
! This subroutine is used to build a brand new cell list
! Created on 7-15-2019 by Kaihang Shi 
! Considering the external structures
! ==========================================================


 	  Subroutine update_clist

 	  Use global

 	  IMPLICIT NONE


 	  ! Local
 	  Integer :: imol, ibox, icel, ibinx, ibiny, ibinz, itype, isite



 	  ! Assume only one simulation box
 	  ibox = 1

 	  ! Initialize head of chain for each cell
 	  Do icel = 1, clist_ncel

 	  	! head of chain for normal molecules
 	  	clist_hoc(icel) = 0
 	  	! head of chain for substrate sites
 	  	clist_hoc_sub(icel) = 0

 	  End Do


 	  ! Loop over all molecules
 	  Do imol = 1, n_mol_tot(ibox)

 	  	! Get the molecule types 
	  	itype = mol_type(imol,ibox)

	  	! Check initial style
	  	If (initstyle(itype,ibox) .eq. 'coords') Then

	  		! Directly loop over imol's sites
	  		Do isite = 1, n_sites(itype)

	  			! Determine the location of substrate sites
 	  			ibinx = FLOOR(rx_s(isite,imol,ibox)/clist_dx) + 1
      			ibiny = FLOOR(ry_s(isite,imol,ibox)/clist_dy) + 1
      			ibinz = FLOOR(rz_s(isite,imol,ibox)/clist_dz) + 1

      			icel = clist_loca(ibinx,ibiny,ibinz)

      			clist_llist_sub(isite) = clist_hoc_sub(icel)
      			clist_hoc_sub(icel) = isite
      		EndDo

      	Else if ((initstyle(itype,ibox) .eq. 'simple_cubic') .OR. (initstyle(itype,ibox) .eq. 'random')) Then

 	  		! Determine the location of imol
 	  		ibinx = FLOOR(rx(imol,ibox)/clist_dx) + 1
      		ibiny = FLOOR(ry(imol,ibox)/clist_dy) + 1
      		ibinz = FLOOR(rz(imol,ibox)/clist_dz) + 1

      		icel = clist_loca(ibinx,ibiny,ibinz)

      		clist_llist(imol) = clist_hoc(icel)
      		clist_hoc(icel) = imol
      	Else

      		Write(*,'(A,I5,A,I5,A,I5)') 'UPDATE_CLIST: invalid initial style', initstyle(itype,ibox),' for mol', imol, ' type', itype
      		STOP

      	Endif

 	  ! End looping over all molecules	
 	  End Do
 	  
 	  
 	  Return
 	  
 	  End Subroutine 













 	  