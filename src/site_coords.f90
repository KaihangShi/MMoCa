! ==========================================================
!  This subroutine is used to calculate and set the coordinates
!  of sites on the molecule
!  Created on 1/6/2017 by Kaihang Shi
!  Reference:
!  [1] Vesely, Journal of Computational Physics 47.2 (1982): 291-296.
!  [2] Charles F. F. Karney, "Quaternions in molecular modeling"
!
!  --comments on different rotation matrix---
!  Through test, different matrix gives the same final results.
!  Concept of 'orientation' of the same quaternion is different 
!  for different set of internal coordinates. But it works equivalently
!  ----------
! ==========================================================


	  Subroutine site_coords(imol,ibox)

	  Use global 

	  IMPLICIT NONE

	  ! Passed 
	  Integer :: imol, ibox

	  ! Local
	  INTEGER :: itype,isite, jsite
	  DOUBLE PRECISION :: a11, a12, a13
	  DOUBLE PRECISION :: a21, a22, a23
	  DOUBLE PRECISION :: a31, a32, a33


	  ! Compute rotation matrix elements (follow the matrix in Vesely paper)
	  a11 = q1(imol,ibox)**2 - q2(imol,ibox)**2 - q3(imol,ibox)**2 + q4(imol,ibox)**2
	  a12 = 2.0d0 * ( q1(imol,ibox) * q2(imol,ibox) - q3(imol,ibox) * q4(imol,ibox) )
	  a13 = 2.0d0 * ( q3(imol,ibox) * q1(imol,ibox) + q2(imol,ibox) * q4(imol,ibox) )

	  a21 = 2.0d0 * ( q1(imol,ibox) * q2(imol,ibox) + q3(imol,ibox) * q4(imol,ibox) )
	  a22 = q2(imol,ibox)**2 - q3(imol,ibox)**2 - q1(imol,ibox)**2 + q4(imol,ibox)**2
	  a23 = 2.0d0 * ( q2(imol,ibox) * q3(imol,ibox) - q1(imol,ibox) * q4(imol,ibox) )

	  a31 = 2.0d0 * ( q3(imol,ibox) * q1(imol,ibox) - q2(imol,ibox) * q4(imol,ibox) )
	  a32 = 2.0d0 * ( q2(imol,ibox) * q3(imol,ibox) + q1(imol,ibox) * q4(imol,ibox) )
	  a33 = q3(imol,ibox)**2 - q1(imol,ibox)**2 - q2(imol,ibox)**2 + q4(imol,ibox)**2

	  ! Compute rotation matrix elements (follow matrix in Charles F. F. Karney's paper)
	  ! Through test, it's equivalent to the matrix in Vesely paper
!	  a11 = q1(imol,ibox)**2 + q2(imol,ibox)**2 - q3(imol,ibox)**2 - q4(imol,ibox)**2
!	  a12 = 2.0d0 * ( q2(imol,ibox) * q3(imol,ibox) - q1(imol,ibox) * q4(imol,ibox) )
!	  a13 = 2.0d0 * ( q2(imol,ibox) * q4(imol,ibox) + q1(imol,ibox) * q3(imol,ibox) )

!	  a21 = 2.0d0 * ( q2(imol,ibox) * q3(imol,ibox) + q1(imol,ibox) * q4(imol,ibox) )
!	  a22 = q1(imol,ibox)**2 - q2(imol,ibox)**2 + q3(imol,ibox)**2 - q4(imol,ibox)**2
!	  a23 = 2.0d0 * ( q3(imol,ibox) * q4(imol,ibox) - q1(imol,ibox) * q2(imol,ibox) )

!	  a31 = 2.0d0 * ( q4(imol,ibox) * q1(imol,ibox) - q1(imol,ibox) * q4(imol,ibox) )
!	  a32 = 2.0d0 * ( q4(imol,ibox) * q3(imol,ibox) + q1(imol,ibox) * q2(imol,ibox) )
!	  a33 = q1(imol,ibox)**2 - q2(imol,ibox)**2 - q3(imol,ibox)**2 + q4(imol,ibox)**2

	  ! Get imol's type
	  itype = mol_type(imol,ibox)

	  ! Calculate the Cartesian coordinates of the sites
	  DO isite = 1,n_sites(itype)

	    rx_s(isite,imol,ibox) = rx(imol,ibox) + a11 * rx_i(isite,itype) + a12 * ry_i(isite,itype) + a13 * rz_i(isite,itype)
	    ry_s(isite,imol,ibox) = ry(imol,ibox) + a21 * rx_i(isite,itype) + a22 * ry_i(isite,itype) + a23 * rz_i(isite,itype)
	    rz_s(isite,imol,ibox) = rz(imol,ibox) + a31 * rx_i(isite,itype) + a32 * ry_i(isite,itype) + a33 * rz_i(isite,itype)

	    ! Initialize site_type
		site_type(isite,itype) = -1

		! Set site type
		Do jsite = 1, n_site_types(itype)
			If (Trim(sites_name(isite,itype)) .eq. site_type_name(jsite,itype)) site_type(isite,itype)=jsite
		End Do

		! Throw an error if no site type is found
		If (site_type(isite,itype) .eq. -1) Then
			Write(*,*) 'FATAL ERROR: NO SITE TYPE IS FOUND FOR SITES LISTED IN INTERNAL COORDINATES.'
			Write(*,*) 'FAILED TO SET SITE TYPE IN site_coords SUBROUTINE'
			STOP
		End If
	        
	  END DO
	  
	  	
	  	Return
	  
	  End Subroutine 