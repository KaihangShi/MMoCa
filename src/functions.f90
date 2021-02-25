! ==========================================================
! This file contains some useful function utilities
! Created on Dec. 4th, 2016 by Kaihang Shi
! ==========================================================



! ==================================================
! Funtion to check labels in input files
! Created on Dec. 4th, 2016 by Kaihang Shi 
! Last modified on Dec. 4th, 2016 by Kaihang Shi
! ==================================================
      Logical Function label_check(ifile,label)


      IMPLICIT NONE

      ! Passed
      Integer :: ifile
      Character(Len=*) :: label

      ! Local
      Character(Len=128) :: string

      ! Read string from ifile
      Read(ifile,'(A)',Err=1) string

      If (label .eq. string) Then
      	label_check = .true.
      Else
      	label_check = .false.
      Endif

      Return	

1     Write(*,*) 'ERROR IN READING FILE ', ifile, ". Please check 'preproc.h' for file number."
	  label_check = .false.

	  Return

	  End Function




! ==================================================
! Subroutine to skip lines while reading data
! Created on Dec. 5th, 2016 by Kaihang Shi
! Last modified on Dec. 5th, 2016
! ==================================================
	  Subroutine skip_lines(ifile,nlines)

	  Use global

	  IMPLICIT NONE

	  ! Passed
	  Integer :: ifile
	  Integer :: nlines

	  ! Local
	  Integer :: iline

	  Do iline = 1,nlines
	  	Read(ifile,*)
	  End Do

	  Return

	  End Subroutine skip_lines



! ==================================================
! Funtion to check if the particle is within the specific range in z direction
! For CG_WALL_FFPW field type specifically
! Created on 11/18/2017 by Kaihang Shi 
! ==================================================
      Subroutine range_check(WITHIN,rzs1,rzs2)

      Use global

      IMPLICIT NONE

      ! Passed
      Double Precision :: rzs1, rzs2
      Logical :: WITHIN

      ! Local
      Double Precision :: zz
      Integer :: ibin1, ibin2

      ! Initialize variable
      ! Assume not in the effective range initially
      WITHIN = .false. 

      ! zz = ABS(rzs1-rzs2)

      ! ! Return if distance larger than a slab size 
      ! If (zz .ge. cg_ff_res) Return

      ! Check if is outside the range
      If (rzs1 .lt. box(3,1)/2.0d0) Then

      	If ((rzs1 .lt. cg_ff_lob) .or. (rzs1 .ge. cg_ff_upb)) Return

      Else
      	zz = box(3,1)-rzs1
      	If ((zz .lt. cg_ff_lob) .or. (zz .ge. cg_ff_upb)) Return

      Endif

      If (rzs2 .lt. box(3,1)/2.0d0) Then

            If ((rzs2 .lt. cg_ff_lob) .or. (rzs2 .ge. cg_ff_upb)) Return

      Else
            zz = box(3,1)-rzs2
            If ((zz .lt. cg_ff_lob) .or. (zz .ge. cg_ff_upb)) Return
      
      Endif

      WITHIN = .TRUE.
  

      ! If (rzs1 .lt. box(3,1)/2.0d0) Then

      !       ! Check if in the same slab
      !       ibin1 = FLOOR((rzs1-cg_ff_lob)/cg_ff_res)
      !       ibin2 = FLOOR((rzs2-cg_ff_lob)/cg_ff_res)
      !       If (ibin1 .eq. ibin2) Then
      !             WITHIN = .true.
      !             Return
      !       Else 
      !             Return      
      !       End If

      ! Else 

      !       ! Check if in the same slab
      !       ibin1 = FLOOR((box(3,1)-rzs1-cg_ff_lob)/cg_ff_res)
      !       ibin2 = FLOOR((box(3,1)-rzs2-cg_ff_lob)/cg_ff_res)
      !       If (ibin1 .eq. ibin2) Then
      !             WITHIN = .true.
      !             Return
      !       Else 
      !             Return      
      !       End If


            
      ! End If
      



      Return
      
	End Subroutine




! ==================================================

! ==================================================




