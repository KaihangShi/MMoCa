! ==========================================================
! This subroutine is used to update Ewald components if MC
! move is accepted 
! Created on 1-11-2017 by Kaihang Shi
! Last modified on 1-11-2017 by Kaihang Shi
! ==========================================================

	  Subroutine update_ewld(selector,ibox,imol)

	  Use global

	  IMPLICIT NONE

	  ! Passed 
	  Integer :: selector, ibox, imol

	  ! Local
	  INTEGER :: lmol

	  ! Update the molecule's vectors

	  ! If not a deletion
	  IF(selector .NE. 3)  THEN

	    ! Update the vectors
	    eikx(ibox,imol,:,:) = eikx_mol(:,:)
	    eiky(ibox,imol,:,:) = eiky_mol(:,:)
	    eikz(ibox,imol,:,:) = eikz_mol(:,:)

	  ELSE

	    ! Get the last molecule in the array
	    lmol = n_mol_tot(ibox) 

	    ! Swap the vectors 
	    eikx(ibox,imol,:,:) = eikx(ibox,lmol,:,:)
	    eiky(ibox,imol,:,:) = eiky(ibox,lmol,:,:)
	    eikz(ibox,imol,:,:) = eikz(ibox,lmol,:,:)

	  ENDIF
	   
	  ! Update the surface correction components
	  ewld_surfc(ibox) = ewld_surfc_new(ibox)

	  ! Update the structure factor
	  skewld = skewld_new
	  skewld_ex = skewld_ex_new
	  
	  
	  	Return
	  
	  End Subroutine 











	  