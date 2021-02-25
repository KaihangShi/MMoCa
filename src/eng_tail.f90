! ==========================================================
! This subroutine is used to calculate vdW tail correction
! between any two molecule types
! Only apply to 12-6 Lennard-Jones now
! Created on 1-24-2017 by Kaihang Shi
! 
! U_tail = 2 * Pi * rho_x * N_y * Integral (rc,infinity,r^2*u(r))
! Note: In this subroutine we didn't include the N*rho prefactor 
! 		where N is number of molecules, rho is number density of 
! 		one molecule type, but we will include all the cross term 
! 		prefactor when we calculate total tail correction energy.
! 
! Reference:
! [1]. Scott Shell's lecture notes
! [2]. Towhee subroutine 'tail.F'
! ==========================================================

      Subroutine eng_tail(ibox,itype,jtype)

      Use global

      IMPLICIT NONE

      ! Passed 
      Integer :: ibox
      Integer :: itype,jtype

      ! Local
      Integer :: isite, jsite, isitetype, jsitetype
      Double Precision :: sigma3, rci3


      ! Initialize global variable
      vdw_tail(itype,jtype) = 0.0d0

      ! Loop over all itype's sites 
      Do isite = 1, n_sites(itype)
      	! Get isite type
      	isitetype = site_type(isite,itype)

      	! Loop over all jtype's sites
      	Do jsite = 1, n_sites(jtype)
      	     ! Get jsite type
      	     jsitetype = site_type(jsite,jtype)

      	     ! Calculate sigma^3
      	     sigma3 = sigma(isitetype,itype,jsitetype,jtype)**3
      	     ! Calculate (sigma/rc)^3
      	     rci3 = (sigma(isitetype,itype,jsitetype,jtype)/r_cut)**3

      	     ! Calculate 2pi*Integral value in the unit of [K*A^3]
      	     vdw_tail(itype,jtype) = vdw_tail(itype,jtype) + &
      		& two_Pi*4.0d0*epsilon(isitetype,itype,jsitetype,jtype)*sigma3*(rci3**3/9.0d0 - rci3/3.0d0)

      	! End loop over jsite
      	End Do
      ! End loop over isite
      End Do

      
      	
      
      	Return
      
      End Subroutine 











