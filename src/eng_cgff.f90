! ==========================================================
! This subroutine handles the effective fluid-fluid interactions
! Created on 11/15/2017 by Kaihang Shi 
! Last modified on 11/20/2017
! ==========================================================

	  Subroutine eng_cgff(ibox,engcgff,rijssq,isite,imol,jsite,jmol)

	  Use global

	  IMPLICIT NONE

	  ! Passed
	  Integer :: ibox,isite,imol,jsite,jmol
	  Double Precision :: engcgff, rijssq

	  ! Local
	  Double Precision :: lj_factor, sig12sq, iloop, sig1, sig2, sig12
	  Integer :: itype, jtype, isitetype, jsitetype, ibin



	  ! Initialize energy
	  engcgff = 0.0d0

	  ! Get the molecule types 
	  itype = mol_type(imol,ibox)
	  jtype = mol_type(jmol,ibox)
	  ! Get isite type 
	  isitetype = site_type(isite,itype)
	  ! Get jsite type 
	  jsitetype = site_type(jsite,jtype)



	  If ((potential .eq. LENNARD_JONES) .and. (epsilon(isitetype,itype,jsitetype,jtype) .ne. 0.0d0)) Then

	  	! Determine sigma value for imol
	  	If (rz_s(isite,imol,ibox) .le. box(3,ibox)/2.0d0) Then

	  		! Loop to see where is the imol 
		  	Do ibin = 1, cg_ff_maxnum

		  		If (rz_s(isite,imol,ibox) .le. eff_sig_z(ibin)) Then

		  			If (ibin .eq. 1) Then
		  				sig1 = eff_sig(1)
		  				EXIT
		  			Else

		  				sig1 = (rz_s(isite,imol,ibox)-eff_sig_z(ibin-1))*&
			 					&((eff_sig(ibin)-eff_sig(ibin-1))/(eff_sig_z(ibin)-eff_sig_z(ibin-1))) + &
			 					& eff_sig(ibin-1)
			 			EXIT

		  			End If

		  		Else If (rz_s(isite,imol,ibox) .gt. eff_sig_z(cg_ff_maxnum)) Then

		  			sig1 = sigma(isitetype,itype,isitetype,itype)
		  			EXIT
		  			
		  		End If

		  			
		  	End Do

		Else 

			! Loop to see where is the imol 
		  	Do ibin = 1, cg_ff_maxnum

		  		If (box(3,ibox)-rz_s(isite,imol,ibox) .le. eff_sig_z(ibin)) Then

		  			If (ibin .eq. 1) Then
		  				sig1 = eff_sig(1)
		  				EXIT
		  			Else

		  				sig1 = (box(3,ibox)-rz_s(isite,imol,ibox)-eff_sig_z(ibin-1))*&
			 					&((eff_sig(ibin)-eff_sig(ibin-1))/(eff_sig_z(ibin)-eff_sig_z(ibin-1))) + &
			 					& eff_sig(ibin-1)
			 			EXIT

		  			End If

		  		Else If (box(3,ibox)-rz_s(isite,imol,ibox) .gt. eff_sig_z(cg_ff_maxnum)) Then

		  			sig1 = sigma(isitetype,itype,isitetype,itype)
		  			EXIT
		  			
		  		End If

		  			
		  	End Do


	  	End If


		! Determine sigma value for jmol
	  	If (rz_s(jsite,jmol,ibox) .le. box(3,ibox)/2.0d0) Then

	  		! Loop to see where is the imol 
		  	Do ibin = 1, cg_ff_maxnum

		  		If (rz_s(jsite,jmol,ibox) .le. eff_sig_z(ibin)) Then

		  			If (ibin .eq. 1) Then
		  				sig2 = eff_sig(1)
		  				EXIT
		  			Else

		  				sig2 = (rz_s(jsite,jmol,ibox)-eff_sig_z(ibin-1))*&
			 					&((eff_sig(ibin)-eff_sig(ibin-1))/(eff_sig_z(ibin)-eff_sig_z(ibin-1))) + &
			 					& eff_sig(ibin-1)
			 			EXIT

		  			End If

		  		Else If (rz_s(jsite,jmol,ibox) .gt. eff_sig_z(cg_ff_maxnum)) Then

		  			sig2 = sigma(jsitetype,jtype,jsitetype,jtype)
		  			EXIT
		  			
		  		End If

		  			
		  	End Do

		Else 

			! Loop to see where is the imol 
		  	Do ibin = 1, cg_ff_maxnum

		  		If (box(3,ibox)-rz_s(jsite,jmol,ibox) .le. eff_sig_z(ibin)) Then

		  			If (ibin .eq. 1) Then
		  				sig2 = eff_sig(1)
		  				EXIT
		  			Else

		  				sig2 = (box(3,ibox)-rz_s(jsite,jmol,ibox)-eff_sig_z(ibin-1))*&
			 					&((eff_sig(ibin)-eff_sig(ibin-1))/(eff_sig_z(ibin)-eff_sig_z(ibin-1))) + &
			 					& eff_sig(ibin-1)
			 			EXIT

		  			End If

		  		Else If (box(3,ibox)-rz_s(jsite,jmol,ibox) .gt. eff_sig_z(cg_ff_maxnum)) Then

		  			sig2 = sigma(jsitetype,jtype,jsitetype,jtype)
		  			EXIT
		  			
		  		End If

		  			
		  	End Do


	  	End If

	  	sig12 = (sig1+sig2)/2.0d0
	  	If ((sig12 .lt. 8.0d0) .and. (sig12 .gt. 4.0d0)) Then

	  		sig12sq = sigmasq(isitetype,itype,jsitetype,jtype)

	  	Else
	  		! Lorentz combining rule 
	  		sig12sq = (sig1+sig2)**2/4.0d0
	  		
	  	End If

	  	


		! Calculate lJ factor
		lj_factor = (sig12sq/rijssq)**3
		
		! Calculate dispersion energy using 12-6 Lennard-Jones potential
		engcgff = engcgff + &
		& 4.0d0*epsilon(isitetype,itype,jsitetype,jtype)*lj_factor*(lj_factor-1.0d0)

      End If


	  ! If (rzs1 .lt. box(3,1)/2.0d0) Then

   !      ! Check if in the same slab
   !      ibin1 = FLOOR((rzs1-cg_ff_lob)/cg_ff_res)
   !      ibin2 = FLOOR((rzs2-cg_ff_lob)/cg_ff_res)

   !      If (ibin1 .ne. ibin2) Then
   !      	Write(*,*) 'Code is wrong in eng_cgff'
   !      	STOP     
   !      End If

   !      If (ibin1+1 .gt. cg_ff_maxnum) Then
   !      	Write(*,*) 'Code is wrong in eng_cgff'
   !      	STOP 
   !      End If

   !      engcgff = cg_ff_egrid(ibin1+1)

   !    Else 

   !      ! Check if in the same slab
   !      ibin1 = FLOOR((box(3,1)-rzs1-cg_ff_lob)/cg_ff_res)
   !      ibin2 = FLOOR((box(3,1)-rzs2-cg_ff_lob)/cg_ff_res)

   !      If (ibin1 .ne. ibin2) Then
   !      	Write(*,*) 'Code is wrong in eng_cgff'
   !      	STOP     
   !      End If

   !      If (ibin1+1 .gt. cg_ff_maxnum) Then
   !      	Write(*,*) 'Code is wrong in eng_cgff'
   !      	STOP 
   !      End If

   !      engcgff = cg_ff_egrid(ibin1+1)


            
      ! End If


	  ! Check if isite is inside the wall
		! If ((rzs(isite,imol,ibox) .le. cg_wall_position(1)) .or. (rzs(isite,imol,ibox) .ge. cg_wall_position(2))) Then
		! 	OVERLAP = .TRUE.
		! 	RETURN
		! End If



		! If (rzs .lt. box(3,ibox)/2.0d0) Then
		! 	! Double check
		! 	If(rzs .lt. cg_wall_initz) Then
		! 		engcgff = engcgff + 1.0d4
		! 		RETURN
		! 	End if

		! 	! Get distance to the surface
	 !  		zz = rzs - cg_wall_position(1)


		! 	! Get bin number for energy grid
		! 	ibin = CEILING((rzs-cg_wall_initz)/cg_wall_res)
		! 	! case 1: ibin=0
		! 	If (ibin .eq. 0) Then
		! 		engcgff = engcgff + cg_ff_egrid(1)
		! 		RETURN
		! 	End If
		! 	! check if beyond the cutoff to avoid segmentation fault
		! 	If ((ibin+1) .gt. cg_ff_maxnum) Then
		! 		Write(*,*) "ibin beyond cg_Ff_maxnum. Double check!!"
		! 		STOP
		! 	End if

		! 	! Do linear interpolation to calculate the energy [K]
		! 	engcgff = engcgff + (zz-cg_ff_dgrid(ibin))*&
		! 		&((cg_ff_egrid(ibin+1)-cg_ff_egrid(ibin))/(cg_ff_dgrid(ibin+1)-cg_ff_dgrid(ibin))) + &
		! 		& cg_ff_egrid(ibin)

		! Else

		! 	! Double check
	 !  		If((box(3,ibox)-rzs) .lt. cg_wall_initz) Then
	 !  			engcgff = engcgff + 1.0d4
	 !  			RETURN
	 !  		End if 

	 ! 		! Get distance to the surface
	 !  		zz =  cg_wall_position(2) - rzs

	  		 
	 !  		! Get bin number for energy grid
	 !  		ibin = CEILING((box(3,ibox)-rzs-cg_wall_initz)/cg_wall_res)
	 !  		! case 1: ibin=0
	 !  		If (ibin .eq. 0) Then
	 !  			engcgff = engcgff + cg_ff_egrid(1)
	 !  			Return
	 !  		End If
	 !  		! check if beyond the cutoff to avoid segmentation fault
		! 	If ((ibin+1) .gt. cg_ff_maxnum) Then
		! 		Write(*,*) "ibin beyond cg_Ff_maxnum. Double check!!"
		! 		STOP
		! 	End if
	  		
	 !  		! Do linear interpolation to calculate the energy [K]
	 !  		engcgff = engcgff + (zz-cg_ff_dgrid(ibin))*&
	 !  			&((cg_ff_egrid(ibin+1)-cg_ff_egrid(ibin))/(cg_ff_dgrid(ibin+1)-cg_ff_dgrid(ibin))) + &
	 !  			& cg_ff_egrid(ibin)
			
		! End If
		
		


	  
	  	
	  
	  	Return
	  
	  End Subroutine 