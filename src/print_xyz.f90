! ==========================================================
! This subroutine is used to print *.xyz configuration data
! Borrowed from Jeremy Palmer's write_config.f90 subroutine
! 
! ==========================================================

      Subroutine print_xyz(ibox,nstep,selector)

      Use global

      IMPLICIT NONE


      ! Passed
      INTEGER :: ibox,nstep
      Character (LEN=*) :: selector

      ! Local
      INTEGER :: imol,itype,isite,n_sites_tot, isitetype, istep
      CHARACTER(LEN=30) :: filename


      Select Case (selector)
            
            
        Case('sites')

            ! Write site positions
      
              ! Initialize the total number of sites
              n_sites_tot = 0

              ! Calculate the total number of atomic sites  
              DO itype = 1,n_mol_types
                n_sites_tot = n_sites_tot + n_sites(itype)*n_mol(itype,ibox)
              ENDDO

              ! Create an appropriate file name
              WRITE(filename,'(A,I0,A,I0,A)') 'config_sites_b',ibox,'_',nstep,'.xyz'

              OPEN(UNIT=FILE_XYZ_SITES,FILE=filename,STATUS='UNKNOWN',ACCESS='SEQUENTIAL',ACTION='WRITE')

                ! Write the total number of molecules
                WRITE(FILE_XYZ_SITES,'(I10)') n_sites_tot

                ! Write the total number of molecules and box dimensions
                WRITE(FILE_XYZ_SITES,'(I0,3F17.9)') n_mol_tot(ibox),box(1,ibox),box(2,ibox),box(3,ibox)

                ! Write the coordinates to file (sorting by type)
                DO itype = 1,n_mol_types

                  ! Loop over the molecules
                  DO imol = 1,n_mol_tot(ibox)

                    ! Check that it is the right type
                    IF(mol_type(imol,ibox) .EQ. itype) THEN

                      ! Loop over imol's sites 
                      DO isite = 1,n_sites(itype)

                        ! Get isite type
                        isitetype = site_type(isite,itype)

                        ! Write the type and coordinates to file
                        WRITE(FILE_XYZ_SITES,'(A4,3F17.9)') &
                        &site_type_name(isitetype,itype),&
                        &rx_s(isite,imol,ibox),ry_s(isite,imol,ibox),rz_s(isite,imol,ibox)
            
                      ENDDO

                    ENDIF

                  ENDDO

                ENDDO

              CLOSE(UNIT=FILE_XYZ_SITES)

        Case('all')

            ! Write site positions
      
              ! Initialize the total number of sites
              n_sites_tot = 0

              ! Calculate the total number of atomic sites  
              DO itype = 1,n_mol_types
                n_sites_tot = n_sites_tot + n_sites(itype)*n_mol(itype,ibox)
              ENDDO

              ! Open file
              OPEN(UNIT=FILE_XYZ_ALL,FILE='dump_xyz.txt',STATUS='UNKNOWN',ACCESS='SEQUENTIAL',ACTION='WRITE',POSITION='APPEND')

                ! Write the total number of molecules
                WRITE(FILE_XYZ_ALL,'(I10,I10)') n_sites_tot, nstep

                ! Write the total number of molecules and box dimensions
                WRITE(FILE_XYZ_ALL,'(I0,3F17.9)') n_mol_tot(ibox),box(1,ibox),box(2,ibox),box(3,ibox)

                ! Write the coordinates to file (sorting by type)
                DO itype = 1,n_mol_types

                  ! Loop over the molecules
                  DO imol = 1,n_mol_tot(ibox)

                    ! Check that it is the right type
                    IF(mol_type(imol,ibox) .EQ. itype) THEN

                      ! Loop over imol's sites 
                      DO isite = 1,n_sites(itype)

                        ! Get isite type
                        isitetype = site_type(isite,itype)

                        ! Write the type and coordinates to file
                        WRITE(FILE_XYZ_ALL,'(A4,3F17.9)') &
                        &site_type_name(isitetype,itype),&
                        &rx_s(isite,imol,ibox),ry_s(isite,imol,ibox),rz_s(isite,imol,ibox)
            
                      ENDDO

                    ENDIF

                  ENDDO

                ENDDO

              CLOSE(UNIT=FILE_XYZ_ALL)

        Case('mol')

            ! Dump specified molecule type to files
            OPEN(UNIT=FILE_XYZ_MOL,FILE='dump_xyz_mol.txt',STATUS='UNKNOWN',ACCESS='SEQUENTIAL',ACTION='WRITE',POSITION='APPEND')
                
                ! Write the number of molecules for specified type
                DO itype = 1,n_mol_types
                    If (dump_moltype .eq. itype) Then
                        WRITE(FILE_XYZ_MOL,'(I10)',ADVANCE='NO') n_mol(itype,ibox)
                    End If
                ENDDO  

                ! Write number of steps
                WRITE(FILE_XYZ_MOL,'(I10)') nstep

                ! Write the box dimensions
                WRITE(FILE_XYZ_MOL,'(3F17.9)') box(1,ibox),box(2,ibox),box(3,ibox) 

                  ! Loop over the molecules
                  DO imol = 1,n_mol_tot(ibox)

                    ! Check that it is the right type
                    IF(mol_type(imol,ibox) .EQ. dump_moltype) THEN

                        ! Only dump center-of-mass of the selected molecule
                        If(dumpxyz_mode .eq. 'center') Then
                    
                            ! Write the type and coordinates to file
                            WRITE(FILE_XYZ_MOL,'(A,3F17.9)') &
                            &mol_type_name(dump_moltype),rx(imol,ibox),ry(imol,ibox),rz(imol,ibox)

                        ! Dump molecular sites 
                        Else if(dumpxyz_mode .eq. 'sites') Then
                            ! Loop over imol's sites 
                            DO isite = 1,n_sites(dump_moltype)

                                ! Get isite type
                                isitetype = site_type(isite,dump_moltype)

                                ! Write the type and coordinates to file
                                WRITE(FILE_XYZ_MOL,'(A4,3F17.9)') &
                                &site_type_name(isitetype,dump_moltype),&
                                &rx_s(isite,imol,ibox),ry_s(isite,imol,ibox),rz_s(isite,imol,ibox)
            
                            ENDDO

                        ENDIF
                    ENDIF

                  ENDDO
                  

              CLOSE(UNIT=FILE_XYZ_MOL)

        ! Added on Oct. 13, 2019
        Case('restart')

            ! Write site positions
      
              ! Initialize the total number of sites
              n_sites_tot = 0

              ! Calculate the total number of atomic sites  
              DO itype = 1,n_mol_types
                n_sites_tot = n_sites_tot + n_sites(itype)*n_mol(itype,ibox)
              ENDDO

              OPEN(UNIT=FILE_RST,FILE='rst1.xyz',STATUS='UNKNOWN',ACCESS='SEQUENTIAL',ACTION='WRITE')

                ! Write the total number of molecules
                WRITE(FILE_RST,'(I10)') n_sites_tot

                ! Write the total number of molecules and box dimensions
                WRITE(FILE_RST,'(I0,3F17.9)') n_mol_tot(ibox),box(1,ibox),box(2,ibox),box(3,ibox)

                ! Write the coordinates to file (sorting by type)
                DO itype = 1,n_mol_types

                  ! Loop over the molecules
                  DO imol = 1,n_mol_tot(ibox)

                    ! Check that it is the right type
                    IF(mol_type(imol,ibox) .EQ. itype) THEN

                      ! Loop over imol's sites 
                      DO isite = 1,n_sites(itype)

                        ! Get isite type
                        isitetype = site_type(isite,itype)

                        ! Write the type and coordinates to file
                        WRITE(FILE_RST,'(A4,3F17.9)') &
                        &site_type_name(isitetype,itype),&
                        &rx_s(isite,imol,ibox),ry_s(isite,imol,ibox),rz_s(isite,imol,ibox)
            
                      ENDDO

                    ENDIF

                  ENDDO

                ENDDO

              CLOSE(UNIT=FILE_RST)

              ! Copy file
              !Call EXECUTE_COMMAND_LINE("cp rst1.xyz rst2.xyz")


        Case default

            Write(*,*) 'FATAL ERROR: INVALID SELECTOR IN print_xyz SUBROUTINE'
            STOP


      
      End Select
        
      
        Return
      
      End Subroutine 








