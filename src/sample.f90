! ==========================================================
! This subroutine is used to perform optional sampling tasks
! Created on 2-4-2017 by Kaihang Shi
! ==========================================================

    
      Subroutine sample(selector,iblock,istep)
      
      Use global

      IMPLICIT NONE

      ! Passed 
      Integer :: selector, iblock, istep

      ! Local
      Integer :: itype, imol, isite, jtype, jmol, jsite, ibin, istp, iblk, ibox, inum
      Double Precision, Dimension(zden_bins,n_mol_types) :: zdenavg
      Double Precision, Dimension(rden_bins,n_mol_types) :: rdenavg
      Double Precision, Dimension(n_mol_types) :: surfexavg, surfexstd
      Double Precision, Dimension(zden_bins) :: pkin
      Double Precision, Dimension(3,zden_bins,n_box) :: virialpnavg
      Double Precision, Dimension(3,5,zden_bins,n_box) :: virialptavg
      Double Precision :: rxis, ryis, rzis, rxijs, ryijs, rzijs, rijssq, clrsq, clr, clrx, clry
      Double Precision :: fluc_N, fluc_U
      Double Precision, Dimension(4) :: qstavg
      Double Precision :: rden_dvol, posr
      Double Precision, Dimension(rden_bins) :: prkin
      ! first rank 1-fluid-wall 2-fluid-fluid 3-total;
      ! second rank 1-IK 2-Harasima
      Double Precision, Dimension(3,2,rden_bins,n_box) :: virialpnravg, virialpttavg, virialptzavg
      ! Virial
      Double Precision :: instvir





      Select Case (selector)
        
        ! Reset sampling statistics
        ! Input: none
        Case(1)
            
            ! z_density
            If (lzdensity) Then
                blk_zden = 0.0d0
                zden_stat = 0.0d0
            End If

            ! r_density
            If (lrdensity) Then
                blk_rden = 0.0d0
                rden_stat = 0.0d0
            End If

            ! surface_excess
            If (lsurfex) Then
                blk_surfex = 0.0d0
            End If

            ! Lattice constant
            If (llattconst) Then
                blk_lattconst = 0.0d0
                lattconst_stat = 0.0d0
            End If

            ! Isosteric heat
            If (lqst) Then
                blk_qst = 0.0d0
            End If

            ! Reset pressure tensor counter
            if(lthermopress_slit) thermopress_slit_stat = 0.0d0
            if(lvirialpress_slit) virialpress_slit_stat = 0.0d0
            if(lvirialpress_cylin) virialpress_cylin_stat = 0.0d0
            if(lvirialpress) virialpress_stat = 0.0d0



        ! Perform optional sampling/accumulate statistics
        ! Input: istep, iblock
        Case(2)
                
            ! Sampling Z-density
            ! ibox = 1, only sample for box 1 right now
            If (lzdensity) Then

                ! Check if the calculation should be attempted
                IF((zden_freq .EQ. 0) .OR. (MOD(istep,zden_freq) .NE. 0)) goto 100

                zden_stat = zden_stat + 1.0d0


                ! Calculate related parameters
                ! For NPT ensemble, the bulk volume is changing
                ! For lfield, parameters have already been initialized in 'initialize.f90'
                If (.NOT. lfield) Then

                    dz = box(3,1)/zden_bins
                    dvol = box(1,1)*box(2,1)*dz
                    
                End If
                
                
                ! Loop over all molecules
                Do imol = 1, n_mol_tot(1)

                    ! Get imol's type
                    itype = mol_type(imol,1)
                    
                    ! Cycle if it is external structure
                    If (initstyle(itype,1) .eq. 'coords') CYCLE

                    IF (lfield .and. (field_type .eq. STEELE_SLIT_PORE)) THEN
                        ! Calculate ibin number
                        ibin = FLOOR((rz(imol,1)-steele_position(1))/dz) + 1

                        ! Accumulate number and will convert to density in case(3)
                        blk_zden(ibin,itype,istep) = blk_zden(ibin,itype,istep) + 1.0d0

                    Else if (lfield .and. (field_type .eq. STEELE_SLIT_FINITEX)) THEN

                        If ((rx(imol,1) .lt. steele_avgx(1)) .or. &
                            & (rx(imol,1) .gt. steele_avgx(2))) CYCLE
                        

                        ! Calculate ibin number
                        ibin = FLOOR((rz(imol,1)-steele_position(1))/dz) + 1

                        ! Accumulate number and will convert to density in case(3)
                        blk_zden(ibin,itype,istep) = blk_zden(ibin,itype,istep) + 1.0d0

                    Else if (lfield .and. (field_type .eq. HARD_SLIT_FINITEX)) Then

                        If ((rx(imol,1) .lt. steele_avgx(1)) .or. &
                            & (rx(imol,1) .gt. steele_avgx(2))) CYCLE
                        

                        ! Calculate ibin number
                        ibin = FLOOR((rz(imol,1)-wall_radius)/dz) + 1

                        ! Accumulate number and will convert to density in case(3)
                        blk_zden(ibin,itype,istep) = blk_zden(ibin,itype,istep) + 1.0d0

                    Else 

                        ! Calculate ibin number
                        ibin = FLOOR(rz(imol,1)/dz) + 1

                        ! Accumulate number and will convert to density in case(3)
                        blk_zden(ibin,itype,istep) = blk_zden(ibin,itype,istep) + 1.0d0
                    Endif

                        
                    
                End Do

                ! Loop over each bin
                Do ibin = 1, zden_bins
                    ! Loop over each molecule type
                    Do itype = 1, n_mol_types

                        ! Exclude external structure
                        If (initstyle(itype,1) .eq. 'coords') CYCLE

                        ! convert to number density (1/A^3)
                        blk_zden(ibin,itype,istep) = blk_zden(ibin,itype,istep)/dvol

                        ! Accumulate statistics
                        avg_zden(ibin,itype,iblock) = avg_zden(ibin,itype,iblock) + &
                                                    & blk_zden(ibin,itype,istep)
                    
                    End Do
                End Do
            ! End sampling z-density
            End If

100         Continue

            ! Sampling surface excess
            If (lsurfex) Then
                
                ! Perform calculation per 10 steps 
                If (MOD(istep,10) .ne. 0) goto 200

                ! Loop over molecule types
                Do itype = 1, n_mol_types
                    
                    ! Cycle if it is external structure
                    If (initstyle(itype,1) .eq. 'coords') CYCLE

                    ! Calculate surface excess [1/A^2]
                    blk_surfex(itype,istep) = (DBLE(n_mol(itype,1)) - surfex_n_mol(itype))/surfex_area

                End Do

            End If

200         Continue
            
            ! Sampling lattice constant
            If (llattconst) Then
                
                ! Check if the calculation should be attempted
                IF((lattconst_freq .EQ. 0) .OR. (MOD(istep,lattconst_freq) .NE. 0)) goto 300

                ! Loop over all molecules
                Do imol = 1, n_mol_tot(1)-1

                    ! Get imol's type
                    itype = mol_type(imol,1)

                    ! Cycle if it is external structure
                    If (initstyle(itype,1) .eq. 'coords') CYCLE   

                    Do isite = 1, n_sites(itype)

                        rxis = rx_s(isite,imol,1)
                        ryis = ry_s(isite,imol,1)
                        rzis = rz_s(isite,imol,1)                       

                        If ((rxis .lt. steele_avgx(1)) .or. &
                            & (rxis .gt. steele_avgx(2))) CYCLE

                        ! Loop over jmol
                        Do jmol = imol+1, n_mol_tot(1)

                            ! Get jmol's type
                            jtype = mol_type(jmol,1)

                            ! loop over jmol's sites (adsorbates)
                            Do jsite = 1, n_sites(jtype)

                                If ((rx_s(jsite,jmol,1) .lt. steele_avgx(1)) .or. &
                                    & (rx_s(jsite,jmol,1) .gt. steele_avgx(2))) CYCLE

                                ! Loop over each confined layer
                                Do inum = 1, lattconst_n

                                    ! If both are in the same layer
                                    If (((rzis .ge. lattconst_cut(1,inum)) .and. (rzis .le. lattconst_cut(2,inum))) &
                                        & .AND. ((rz_s(jsite,jmol,1) .ge. lattconst_cut(1,inum)) &
                                        & .and. (rz_s(jsite,jmol,1) .le. lattconst_cut(2,inum)))) Then                                          

                                        !Calculate the site separation vector 
                                        rxijs = rx_s(jsite,jmol,1) - rxis
                                        ryijs = ry_s(jsite,jmol,1) - ryis
                                        !rzijs = rz_s(jsite,jmol,1) - rzis                                      

                                        ! Apply minimum image convension
                                        ! For steele_slit_finitex, only apply to y-direction                        
                                        ryijs = ryijs - dNINT(ryijs/box(2,1))*box(2,1)

                                        ! Calculate the square of the site scalar separation distance
                                        rijssq = rxijs**2 + ryijs**2 

                                        ! 
                                        If (rijssq .le. 20.25d0) Then
                                            blk_lattconst(inum,istep) = blk_lattconst(inum,istep) + dSQRT(rijssq)
                                            lattconst_stat(inum) = lattconst_stat(inum) + 1.0d0
                                        Else

                                            EXIT
                                        End If

                                    ! If both are not in the same layer
                                    Else
                                        CYCLE                                                                                       
                                    End If
                                ! End loop over confined layers
                                End Do

                            ! End loop over jmol's site
                            End Do
                        ! End loop over jmols
                        End Do
                    ! End loop isite
                    End Do
                ! End loop imol
                End Do

            End If

300         Continue

            ! Sampling r-density
            ! ibox = 1, only sample for box 1 right now
            If (lrdensity) Then

                ! Check if the calculation should be attempted
                IF((rden_freq .EQ. 0) .OR. (MOD(istep,rden_freq) .NE. 0)) goto 400

                rden_stat = rden_stat + 1.0d0
                
                ! Loop over all molecules
                Do imol = 1, n_mol_tot(1)

                    ! Get imol's type
                    itype = mol_type(imol,1)
                    
                    ! Cycle if it is external structure
                    If (initstyle(itype,1) .eq. 'coords') CYCLE

                    ! Move molecules into a image box with (0,0,0) at the center of its bottom surface 
                    clrx = rx(imol,1) - 0.5d0*box(1,1) 
                    clry = ry(imol,1) - 0.5d0*box(2,1)                          

                    ! Calculate r-distance of imol in cylindrical coordiantes
                    clrsq = clrx**2 + clry**2
                    clr = dSQRT(clrsq)

                    If (clr .lt. rden_cut) Then
                        ! Calculate ibin number
                        ibin = FLOOR(clr/rden_dr) + 1

                        ! Accumulate number and will convert to density in case(3)
                        blk_rden(ibin,itype,istep) = blk_rden(ibin,itype,istep) + 1.0d0
                    End If

                ! End loop over all molecules   
                End Do


                ! Loop over bins
                Do ibin = 1, rden_bins
                    ! Loop over molecule types
                    Do itype = 1, n_mol_types

                        ! Exclude external structure
                        If (initstyle(itype,1) .eq. 'coords') CYCLE

                        ! Calculate volume for each bin (changing in NPT ensemble)
                        rden_dvol = Pi*DBLE(2*ibin-1)*box(3,1)*rden_drsq

                        ! Convert to number density (1/A^3)
                        blk_rden(ibin,itype,istep) = blk_rden(ibin,itype,istep)/rden_dvol
                        
                        ! Accumulate statistics
                        avg_rden(ibin,itype,iblock) = avg_rden(ibin,itype,iblock) + &
                                                    & blk_rden(ibin,itype,istep)
                    
                    End Do              
                End Do              
            

            End If

400         Continue
            
            ! Isosteric heat of adsorption
            ! Input: istep
            If (lqst) Then

                ! Assume one box now
                ibox = 1

                ! Loop over molecule types
                Do itype = 1, n_mol_types
                    
                    ! Cycle if it is external structure
                    If (initstyle(itype,1) .eq. 'coords') CYCLE

                    ! Accumulate statistics
                    ! Number of molecules <N>
                    !blk_qst(1,istep,ibox) = n_mol(itype,ibox)
                    ! Configurational energy of the system <U>
                    blk_qst(2,istep,ibox) = energy(ibox)
                    ! <N^2>
                    blk_qst(3,istep,ibox) = DBLE(n_mol(itype,ibox))*DBLE(n_mol(itype,ibox))
                    ! <U*N> 
                    blk_qst(4,istep,ibox) = energy(ibox)*DBLE(n_mol(itype,ibox))


                End Do

                
            End If

            ! Dump density to file
            If (ldumpdensity) Then
                ! Check if dump should be attempted
                IF((dumpdensity_freq .EQ. 0) .OR. (MOD(istep,dumpdensity_freq) .NE. 0)) goto 500

                ! Open file
                OPEN(UNIT=FILE_DENSITY,FILE='density.txt', &
                        & STATUS='UNKNOWN',ACCESS='SEQUENTIAL',ACTION='WRITE',POSITION='APPEND')
     
!               Write(FILE_DENSITY,*) ' Steps            density [#molecule (in box 1)/volume[A^3]]'


                Write(FILE_DENSITY,'(I10,7X,E15.7)') &
                    & (iblock-1)*block_size+istep, DBLE(n_mol_tot(1))/vol(1)

                
                ! Close file
                Close(UNIT=FILE_DENSITY)

            End If



500         Continue

            
            ! Dump total energy to file
            If (ldumpenergy) Then
                ! Check if dump should be attempted
                IF((dumpenergy_freq .EQ. 0) .OR. (MOD(istep,dumpenergy_freq) .NE. 0)) goto 600

                ! Open file
                OPEN(UNIT=FILE_ENERGY,FILE='energy.txt', &
                        & STATUS='UNKNOWN',ACCESS='SEQUENTIAL',ACTION='WRITE',POSITION='APPEND')
     
!               Write(FILE_ENERGY,*) ' Steps            Energy [J]'


                Write(FILE_ENERGY,'(I10,7X,E18.11)') &
                    & (iblock-1)*block_size+istep, energy(1) * Kb

                
                ! Close file
                Close(UNIT=FILE_ENERGY)

            End If

600         Continue
            

            ! Dump virial to file (Added on Jan, 2020)
            If (ldumpvir) Then
                ! Check if dump should be attempted
                IF((dumpvir_freq .EQ. 0) .OR. (MOD(istep,dumpvir_freq) .NE. 0)) goto 700


                ! Calculate instant virial [J] (A&T 2.61)
                instvir = -1.0/3.0*vir(1)*Kb

                ! Open file
                OPEN(UNIT=FILE_VIR,FILE='virial.txt', &
                        & STATUS='UNKNOWN',ACCESS='SEQUENTIAL',ACTION='WRITE',POSITION='APPEND')
     
!               Write(FILE_VIR,*) ' Instant hydrostatic pressure [bar]       Virial (A&T 2.61) [J]     Hypervirial (A&T 2.77) [J]'
                
                ! Assume only one simulation box
                Write(FILE_VIR,'(E18.11,7X,E18.11,7X,E18.11)') &
                    & n_mol_tot(1)/vol(1)*temp*PCOEFF + instvir/vol(1)*1.0d25, instvir, 1.0/9.0*vir_hyp(1)*Kb

                
                ! Close file
                Close(UNIT=FILE_VIR)

            End If

700         Continue



        ! Average statistics for iblock
        ! Input: iblock
        Case(3)

            ! z-density
            If (lzdensity) Then
                
                ! Loop over bins
                Do ibin = 1, zden_bins
                    ! Loop over molecule types
                    Do itype = 1, n_mol_types

                        ! Exclude external structure
                        If (initstyle(itype,1) .eq. 'coords') CYCLE
                        
                        ! average density 
                        avg_zden(ibin,itype,iblock) = avg_zden(ibin,itype,iblock)/zden_stat
                                                    
                    
                    End Do              
                End Do              

            End If


            ! Surface_excess
            If (lsurfex) Then
                
                ! Loop over block_size
                Do istp = 1, block_size

                    ! Cycle if not multiple of 10
                    If (MOD(istp,10) .ne. 0) CYCLE

                    ! Loop over molecule types
                    Do itype = 1, n_mol_types

                        ! Cycle if it is external structure
                        If (initstyle(itype,1) .eq. 'coords') CYCLE

                        ! Calculate average
                        avg_surfex(itype,iblock) = avg_surfex(itype,iblock) + blk_surfex(itype,istp)

                    End Do
                End Do

                ! Loop over moelcule types
                Do itype = 1, n_mol_types

                    ! Cycle if it is external structure
                    If (initstyle(itype,1) .eq. 'coords') CYCLE

                    ! Calculate average
                    avg_surfex(itype,iblock) = avg_surfex(itype,iblock)/(DBLE(block_size)/10.0d0)

                End Do

            End If

            ! Lattice constant
            If (llattconst) Then

                ! Loop over block_size
                Do istp = 1, block_size

                    ! Check if the calculation should be attempted
                    IF((lattconst_freq .EQ. 0) .OR. (MOD(istp,lattconst_freq) .NE. 0)) CYCLE

                    ! loop confined layers
                    Do inum = 1, lattconst_n
                            
                        ! Accumulate statistics
                        avg_lattconst(inum,iblock) = avg_lattconst(inum,iblock) + &
                                                        & blk_lattconst(inum,istp)
                                        
                    End Do              
                End Do

                ! Averaging 
                Do inum = 1, lattconst_n
                    avg_lattconst(inum,iblock) = avg_lattconst(inum,iblock)/lattconst_stat(inum)
                End Do

                
            End If


            ! r-density
            If (lrdensity) Then

                ! Loop over each bin
                Do ibin = 1, rden_bins
                    ! Loop over each molecule type
                    Do itype = 1, n_mol_types

                        ! Exclude external structure
                        If (initstyle(itype,1) .eq. 'coords') CYCLE

                        ! Average and convert to number density (1/A^3)
                        avg_rden(ibin,itype,iblock) = avg_rden(ibin,itype,iblock)/rden_stat
                        
                    End Do
                End Do

            End If


            ! Isosteric heat of adsorption
            ! Input: iblock
            If (lqst) Then
                
                ! Loop over every step in iblock
                Do istp = 1, block_size

                    ! Loop over each box
                    Do ibox = 1, n_box

                        ! Loop over each molecule type
                        Do itype = 1, n_mol_types

                            ! Cycle if it is external structure
                            If (initstyle(itype,1) .eq. 'coords') CYCLE

                            ! <U> in units of K
                            avg_qst(2,iblock,ibox) = avg_qst(2,iblock,ibox) + &
                                                        & blk_qst(2,istp,ibox)
                            ! <N^2>
                            avg_qst(3,iblock,ibox) = avg_qst(3,iblock,ibox) + &
                                                        & blk_qst(3,istp,ibox)
                            ! <U*N>
                            avg_qst(4,iblock,ibox) = avg_qst(4,iblock,ibox) + &
                                                        & blk_qst(4,istp,ibox)
                            
                        End Do

                    End Do
                    
                End Do

                ! Compute the average
                ! Loop over boxes
                Do ibox = 1, n_box

                    ! Loop over molecule types
                    Do itype = 1, n_mol_types

                        ! Cycle if it is external structure
                        If (initstyle(itype,1) .eq. 'coords') CYCLE

                        ! Number of molecules <N>
                        avg_qst(1,iblock,ibox) = avg_nmol(iblock,itype,ibox)
                        ! Configurational energy of the system <U> in units of [K]
                        avg_qst(2,iblock,ibox) = avg_qst(2,iblock,ibox)/DBLE(block_size)
                        ! <N^2>
                        avg_qst(3,iblock,ibox) = avg_qst(3,iblock,ibox)/DBLE(block_size)
                        ! <U*N>
                        avg_qst(4,iblock,ibox) = avg_qst(4,iblock,ibox)/DBLE(block_size)


                    End Do

        
                End Do
            End If




        ! Write statistics to file
        ! Input: none
        Case(4)
        
            ! z-density
            If (lzdensity) Then
                
                ! Open file
                OPEN(UNIT=FILE_ZDENSITY,FILE='z-density.txt',STATUS='UNKNOWN',ACCESS='SEQUENTIAL',ACTION='WRITE')

                ! Initialize 
                zdenavg = 0.0d0

                ! Loop over blocks 
                Do iblk = n_blocks_equil+1, n_blocks_tot
                    
                    ! Loop over bins
                    Do ibin = 1, zden_bins

                        ! Loop over molecule types
                        Do itype = 1, n_mol_types
                            
                            ! Exclude external structure
                            If (initstyle(itype,1) .eq. 'coords') CYCLE

                            zdenavg(ibin,itype) = zdenavg(ibin,itype) + avg_zden(ibin,itype,iblk)

                        End Do

                    End Do

                End Do

                ! Loop over bins
                Do ibin = 1, zden_bins

                    ! Loop over molecule types
                    Do itype = 1, n_mol_types
                        
                        ! Exclude external structure
                        If (initstyle(itype,1) .eq. 'coords') CYCLE

                        ! Calculate averages 
                        zdenavg(ibin,itype) = zdenavg(ibin,itype)/DBLE(n_blocks_prod)

                        ! Convert unit from [1/A^3] to [g/ml]
                        zdenavg(ibin,itype) = (mol_mass(itype)/(Na*1.0d-24))*zdenavg(ibin,itype)

                    End Do
                    
                End Do

                ! Write data to file
                ! loop over molecule types
                Do itype = 1, n_mol_types

                    ! Exclude external structure
                    If (initstyle(itype,1) .eq. 'coords') CYCLE

                    ! Write molecule info
                    Write(FILE_ZDENSITY,'(2A)') 'Molecule name: ', mol_type_name(itype)
                    Write(FILE_ZDENSITY,*) ' z             z-rho [g/ml]           z-rho [1/A^3]'

                    ! Loop over bins
                    Do ibin = 1, zden_bins

                        ! Specialize for graphene slit pore 
                        IF (lfield .and. ((field_type .eq. STEELE_SLIT_PORE) .OR. &
                                & (field_type .eq. STEELE_SLIT_FINITEX))) THEN

                            ! Write to file
                            Write(FILE_ZDENSITY,'(F8.4,7X,E15.7,8X,E15.7)')  (DBLE(ibin)-0.5d0)*dz+steele_position(1),&
                                                        & zdenavg(ibin,itype), (zdenavg(ibin,itype)/mol_mass(itype))*Na*1.0d-24

                        Else if (lfield .and. (field_type .eq. HARD_SLIT_FINITEX)) Then

                            ! Write to file
                            Write(FILE_ZDENSITY,'(F8.4,7X,E15.7,8X,E15.7)')  (DBLE(ibin)-0.5d0)*dz+wall_radius,&
                                                        & zdenavg(ibin,itype), (zdenavg(ibin,itype)/mol_mass(itype))*Na*1.0d-24

                        Else 
                            ! Write to file
                            ! now dz is calculated based on the final configuration
                            Write(FILE_ZDENSITY,'(F8.4,7X,E15.7,8X,E15.7)')  (DBLE(ibin)-0.5d0)*dz, zdenavg(ibin,itype), &
                                                                & (zdenavg(ibin,itype)/mol_mass(itype))*Na*1.0d-24
                        End if

                    End Do
                End Do

                ! Close file
                CLOSE(UNIT=FILE_ZDENSITY)
            
            ! End of z_density section  
            End If


            ! r-density
            If (lrdensity) Then
                
                ! Open file
                OPEN(UNIT=FILE_RDENSITY,FILE='r-density.txt',STATUS='UNKNOWN',ACCESS='SEQUENTIAL',ACTION='WRITE')

                ! Initialize 
                rdenavg = 0.0d0

                ! Loop over blocks 
                Do iblk = n_blocks_equil+1, n_blocks_tot
                    
                    ! Loop over bins
                    Do ibin = 1, rden_bins

                        ! Loop over molecule types
                        Do itype = 1, n_mol_types
                            
                            ! Exclude external structure
                            If (initstyle(itype,1) .eq. 'coords') CYCLE

                            rdenavg(ibin,itype) = rdenavg(ibin,itype) + avg_rden(ibin,itype,iblk)

                        End Do

                    End Do

                End Do

                ! Loop over bins
                Do ibin = 1, rden_bins

                    ! Loop over molecule types
                    Do itype = 1, n_mol_types
                        
                        ! Exclude external structure
                        If (initstyle(itype,1) .eq. 'coords') CYCLE

                        ! Calculate averages 
                        rdenavg(ibin,itype) = rdenavg(ibin,itype)/DBLE(n_blocks_prod)

                        ! Convert unit from [1/A^3] to [g/ml]
                        rdenavg(ibin,itype) = (mol_mass(itype)/(Na*1.0d-24))*rdenavg(ibin,itype)

                    End Do
                    
                End Do

                ! Write data to file
                ! loop over molecule types
                Do itype = 1, n_mol_types

                    ! Exclude external structure
                    If (initstyle(itype,1) .eq. 'coords') CYCLE

                    ! Write molecule info
                    Write(FILE_RDENSITY,'(2A)') 'Molecule name: ', mol_type_name(itype)
                    Write(FILE_RDENSITY,*) ' R             R-rho [g/ml]           R-rho [1/A^3]'

                    ! Loop over bins
                    Do ibin = 1, rden_bins

                        ! Write to file
                        Write(FILE_RDENSITY,'(F8.4,7X,E15.7,8X,E15.7)')  (DBLE(ibin)-0.5d0)*rden_dr, rdenavg(ibin,itype), &
                                                            & (rdenavg(ibin,itype)/mol_mass(itype))*Na*1.0d-24
                        

                    End Do
                End Do

                ! Close file
                CLOSE(UNIT=FILE_RDENSITY)
            
            ! End of r_density section  
            End If




            ! Cylindrical Pressure tensor from virial route
            If (lvirialpress_cylin) Then

                ! Initialize 
                virialpnravg = 0.0d0
                virialpttavg = 0.0d0
                virialptzavg = 0.0d0
                prkin = 0.0d0

                ! Only sample for box 1 right now
                ibox = 1

                ! Loop over iblocks
                Do iblk = n_blocks_equil+1, n_blocks_tot
                    ! Loop over bins
                    Do ibin = 1, rden_bins

                        ! Harasima route
                        If (virialpress_ctype .eq. 2) Then
                            ! Fluid-fluid contribution
                            virialpnravg(2,2,ibin,ibox) = virialpnravg(2,2,ibin,ibox) + virialpress_cylin_pnr(2,2,ibin,iblk,ibox)
                            virialpttavg(2,2,ibin,ibox) = virialpttavg(2,2,ibin,ibox) + virialpress_cylin_ptt(2,2,ibin,iblk,ibox)
                            virialptzavg(2,2,ibin,ibox) = virialptzavg(2,2,ibin,ibox) + virialpress_cylin_ptz(2,2,ibin,iblk,ibox)
                            
                        End If

                    End Do

                End Do

                ! Loop over bins
                Do ibin = 1, rden_bins

                    ! Get r-distance of ibin in cylindrical system
                    posr = (DBLE(ibin)-0.5d0)*rden_dr

                    ! Calculate averages 
                    virialpnravg(2,2,ibin,ibox) = virialpnravg(2,2,ibin,ibox)/DBLE(n_blocks_prod)
                    virialpttavg(2,2,ibin,ibox) = virialpttavg(2,2,ibin,ibox)/DBLE(n_blocks_prod)
                    virialptzavg(2,2,ibin,ibox) = virialptzavg(2,2,ibin,ibox)/DBLE(n_blocks_prod)

                    ! Calculate final pressure (configurational part) in unit of [K/A^3]
                    virialpnravg(2,2,ibin,ibox) = -1.0d0/(2.0d0*two_Pi*posr)*virialpnravg(2,2,ibin,ibox)
                    virialpttavg(2,2,ibin,ibox) = -1.0d0/(two_Pi*posr*rden_dr)*virialpttavg(2,2,ibin,ibox)
                    virialptzavg(2,2,ibin,ibox) = -1.0d0/(two_Pi*posr*rden_dr)*virialptzavg(2,2,ibin,ibox)


                    ! Calculate totoal pressure tensor in unit of [K/A^3]
                    virialpnravg(3,2,ibin,ibox) = virialpnravg(2,2,ibin,ibox)
                    virialpttavg(3,2,ibin,ibox) = virialpttavg(2,2,ibin,ibox)
                    virialptzavg(3,2,ibin,ibox) = virialptzavg(2,2,ibin,ibox)
                    
                End Do

                ! Caculate kinetic part contribution to the pressure
                ! Loop over bins
                Do ibin = 1, rden_bins
                    ! Loop over molecule types
                    Do itype = 1, n_mol_types
                        ! Pressure (Kinetic part) in unit of [K/A^3]
                        prkin(ibin) = prkin(ibin) + (rdenavg(ibin,itype)/mol_mass(itype))*Na*1.0d-24*temp
                    End Do
                End Do

                ! Harasima route
                If (virialpress_ctype .eq. 2) Then

                    ! Open file
                    OPEN(UNIT=FILE_VIRIALPRESS_CYLINH,FILE='press_virial_cylinH.txt', &
                         & STATUS='UNKNOWN',ACCESS='SEQUENTIAL',ACTION='WRITE')

                    ! Write file head
                    Write(FILE_VIRIALPRESS_CYLINH,*) 'Cylindrical pressure tensor from virial route using Harasima definition'
                    Write(FILE_VIRIALPRESS_CYLINH,*) 'Unit: pressure in [bar] and length in [Angstrom]'
                    Write(FILE_VIRIALPRESS_CYLINH,'(A)') &
                                    & '  R               Pkin          Pnr(ff)             Pnr(tot)          Ptz(ff)&
                                    &         Ptz(tot)         Ptt(ff)          Ptt(tot)'

                    ! Write data to file
                    ! Loop over bins
                    Do ibin = 1, rden_bins
                        
                        ! Write to file
                        Write(FILE_VIRIALPRESS_CYLINH, &
                            & '(F8.4,F16.4,F16.4,3X,F16.4,3X,F16.4,F16.4,F16.4,F16.4)')  &
                            & (DBLE(ibin)-0.5d0)*rden_dr, prkin(ibin)*PCOEFF, &
                            & virialpnravg(2,2,ibin,ibox)*PCOEFF, &
                            & (virialpnravg(3,2,ibin,ibox)+prkin(ibin))*PCOEFF, &
                            & virialptzavg(2,2,ibin,ibox)*PCOEFF, &
                            & (virialptzavg(3,2,ibin,ibox)+prkin(ibin))*PCOEFF, &
                            & virialpttavg(2,2,ibin,ibox)*PCOEFF, &
                            & (virialpttavg(3,2,ibin,ibox)+prkin(ibin))*PCOEFF

                    End Do
                
                ! End integral contour definition
                End If
                

            ! End virial Cylindrical pressure tensor    
            End If



            ! Planar surface pressure tensor from virial route
            If (lvirialpress) Then

                ! Initialize 
                virialpnavg = 0.0d0
                virialptavg = 0.0d0
                pkin = 0.0d0
                ! Only sample for box 1 right now
                ibox = 1

                ! Loop over iblocks
                Do iblk = n_blocks_equil+1, n_blocks_tot
                    ! Loop over bins
                    Do ibin = 1, zden_bins
                        ! Accumulate    
                        virialpnavg(1,ibin,ibox) = virialpnavg(1,ibin,ibox) + virialpress_pn(ibin,iblk,ibox)

                        ! IK route
                        If (virialpress_ctype .ne. 2) Then
                            virialptavg(1,1,ibin,ibox) = virialptavg(1,1,ibin,ibox) + virialpress_pt(1,ibin,iblk,ibox)
                            
                        Endif
                        ! Harasima route
                        If (virialpress_ctype .ne. 1) Then
                            virialptavg(1,2,ibin,ibox) = virialptavg(1,2,ibin,ibox) + virialpress_pt(2,ibin,iblk,ibox)
                            
                        End If

                    End Do

                End Do

                ! Loop over bins
                Do ibin = 1, zden_bins

                    ! Calculate averages 
                    virialpnavg(1,ibin,ibox) = virialpnavg(1,ibin,ibox)/DBLE(n_blocks_prod)
                    virialptavg(1,1,ibin,ibox) = virialptavg(1,1,ibin,ibox)/DBLE(n_blocks_prod)
                    virialptavg(1,2,ibin,ibox) = virialptavg(1,2,ibin,ibox)/DBLE(n_blocks_prod)
                    

                    ! Calculate final pressure (configurational part) in unit of [K/A^3]
                    virialpnavg(1,ibin,ibox) = -1.0d0*virialpnavg(1,ibin,ibox)
                    virialptavg(1,1,ibin,ibox) = -0.5d0*virialptavg(1,1,ibin,ibox)                      
                    virialptavg(1,2,ibin,ibox) = -0.5d0*virialptavg(1,2,ibin,ibox)                              

                    
                End Do

                ! Caculate kinetic part contribution to the pressure
                ! Loop over bins
                Do ibin = 1, zden_bins
                    ! Loop over molecule types
                    Do itype = 1, n_mol_types
                        ! Pressure (Kinetic part) in unit of [K/A^3]
                        pkin(ibin) = pkin(ibin) + (zdenavg(ibin,itype)/mol_mass(itype))*Na*1.0d-24*temp
                    End Do
                End Do

                ! IK route
                If (virialpress_ctype .ne. 2) Then
                    ! Open file
                    OPEN(UNIT=FILE_VIRIALPRESS_1,FILE='press_1.txt',STATUS='UNKNOWN',ACCESS='SEQUENTIAL',ACTION='WRITE')

                    ! Write file head
                    Write(FILE_VIRIALPRESS_1,*) 'Pressure tensor from virial route using Irving-Kirkwood definition'
                    Write(FILE_VIRIALPRESS_1,*) 'Unit: pressure in [bar] and length in [Angstrom]'
                    Write(FILE_VIRIALPRESS_1,'(A)') &
                                    & '  z              Pkin               Pt(ff)             Pt(tot)        &
                                    &    Pn(ff)             Pn(tot)'

                    ! Write data to file
                    ! Loop over bins
                    Do ibin = 1, zden_bins

                        ! Write to file
                        Write(FILE_VIRIALPRESS_1,'(F8.4,1X,F14.4,3X,F16.4,3X,F16.4,3X,F16.4,3X,F16.4)')&
                            & (DBLE(ibin)-0.5d0)*dz, pkin(ibin)*PCOEFF, &
                            & virialptavg(1,1,ibin,ibox)*PCOEFF, &
                            & (virialptavg(1,1,ibin,ibox)+pkin(ibin))*PCOEFF, virialpnavg(1,ibin,ibox)*PCOEFF, &
                            & (virialpnavg(1,ibin,ibox)+pkin(ibin))*PCOEFF

                    End Do
                Endif

                ! Harasima route
                If (virialpress_ctype .ne. 1) Then

                    ! Open file
                    OPEN(UNIT=FILE_VIRIALPRESS_2,FILE='press_2.txt', &
                         & STATUS='UNKNOWN',ACCESS='SEQUENTIAL',ACTION='WRITE')

                    ! Write file head
                    Write(FILE_VIRIALPRESS_2,*) 'Pressure tensor from virial route using Harasima definition'
                    Write(FILE_VIRIALPRESS_2,*) 'Unit: pressure in [bar] and length in [Angstrom]'
                    Write(FILE_VIRIALPRESS_2,'(A)') &
                                    & '  z              Pkin               Pt(ff)             Pt(tot)        &
                                    &    Pn(ff)             Pn(tot)'

                    ! Write data to file
                    ! Loop over bins
                    Do ibin = 1, zden_bins

                        ! Write to file
                        Write(FILE_VIRIALPRESS_2,'(F8.4,1X,F14.4,3X,F16.4,3X,F16.4,3X,F16.4,3X,F16.4)')&
                            & (DBLE(ibin)-0.5d0)*dz, pkin(ibin)*PCOEFF, &
                            & virialptavg(1,2,ibin,ibox)*PCOEFF, &
                            & (virialptavg(1,2,ibin,ibox)+pkin(ibin))*PCOEFF, virialpnavg(1,ibin,ibox)*PCOEFF, &
                            & (virialpnavg(1,ibin,ibox)+pkin(ibin))*PCOEFF

                    End Do
                
                ! End integral contour definition
                End If
                

            ! End virial pressure tensor of planar surface
            End If



            !-------------------------------------------------------
            ! slit pore pressure tensor from virial route
            If (lvirialpress_slit) Then

                ! If ((.not. lfield) .and. (field_type .ne. STEELE_SLIT_PORE)) Then
                !     Write(*,*) 'FATAL ERROR: Current virial pressure method can only be applied to steele slit pore model'
                !     Write(*,*) '             and box 1. Please revise the code if want to apply to other systems!'
                !     STOP
                ! End if

                ! Initialize 
                virialpnavg = 0.0d0
                virialptavg = 0.0d0
                pkin = 0.0d0
                ! Only sample for box 1 right now
                ibox = 1

                ! Loop over iblocks
                Do iblk = n_blocks_equil+1, n_blocks_tot
                    ! Loop over bins
                    Do ibin = 1, virialpress_slit_bins
                        ! Accumulate    
                        virialpnavg(1,ibin,ibox) = virialpnavg(1,ibin,ibox) + virialpress_slit_pn(1,ibin,iblk,ibox)
                        virialpnavg(2,ibin,ibox) = virialpnavg(2,ibin,ibox) + virialpress_slit_pn(2,ibin,iblk,ibox)

                        ! IK route
                        If ((virialpress_ctype .EQ. 1) .OR. (virialpress_ctype .EQ. 3) .OR. (virialpress_ctype .EQ. 6)) Then
                            virialptavg(1,1,ibin,ibox) = virialptavg(1,1,ibin,ibox) + virialpress_slit_pt(1,1,ibin,iblk,ibox)
                            virialptavg(2,1,ibin,ibox) = virialptavg(2,1,ibin,ibox) + virialpress_slit_pt(2,1,ibin,iblk,ibox)
                        Endif
                        ! Harasima route
                        If ((virialpress_ctype .EQ. 2) .OR. (virialpress_ctype .EQ. 3) .OR. (virialpress_ctype .EQ. 6)) Then
                            virialptavg(1,2,ibin,ibox) = virialptavg(1,2,ibin,ibox) + virialpress_slit_pt(1,2,ibin,iblk,ibox)
                            virialptavg(2,2,ibin,ibox) = virialptavg(2,2,ibin,ibox) + virialpress_slit_pt(2,2,ibin,iblk,ibox)
                        End If
                        ! H-VR1 route
                        If ((virialpress_ctype .EQ. 4) .OR. (virialpress_ctype .EQ. 6)) Then
                            virialptavg(1,3,ibin,ibox) = virialptavg(1,3,ibin,ibox) + virialpress_slit_pt(1,3,ibin,iblk,ibox)
                            virialptavg(2,3,ibin,ibox) = virialptavg(2,3,ibin,ibox) + virialpress_slit_pt(2,3,ibin,iblk,ibox)
                        End If
                        ! IK-VR1 route
                        If ((virialpress_ctype .EQ. 5) .OR. (virialpress_ctype .EQ. 6)) Then
                            virialptavg(1,4,ibin,ibox) = virialptavg(1,4,ibin,ibox) + virialpress_slit_pt(1,4,ibin,iblk,ibox)
                            virialptavg(2,4,ibin,ibox) = virialptavg(2,4,ibin,ibox) + virialpress_slit_pt(2,4,ibin,iblk,ibox)
                        End If

                    End Do

                End Do

                ! Loop over bins
                Do ibin = 1, virialpress_slit_bins

                    ! Calculate averages 
                    virialpnavg(1,ibin,ibox) = virialpnavg(1,ibin,ibox)/DBLE(n_blocks_prod)
                    virialpnavg(2,ibin,ibox) = virialpnavg(2,ibin,ibox)/DBLE(n_blocks_prod)
                    virialptavg(1,1,ibin,ibox) = virialptavg(1,1,ibin,ibox)/DBLE(n_blocks_prod)
                    virialptavg(2,1,ibin,ibox) = virialptavg(2,1,ibin,ibox)/DBLE(n_blocks_prod)
                    virialptavg(1,2,ibin,ibox) = virialptavg(1,2,ibin,ibox)/DBLE(n_blocks_prod)
                    virialptavg(2,2,ibin,ibox) = virialptavg(2,2,ibin,ibox)/DBLE(n_blocks_prod)
                    virialptavg(1,3,ibin,ibox) = virialptavg(1,3,ibin,ibox)/DBLE(n_blocks_prod)
                    virialptavg(2,3,ibin,ibox) = virialptavg(2,3,ibin,ibox)/DBLE(n_blocks_prod)
                    virialptavg(1,4,ibin,ibox) = virialptavg(1,4,ibin,ibox)/DBLE(n_blocks_prod)
                    virialptavg(2,4,ibin,ibox) = virialptavg(2,4,ibin,ibox)/DBLE(n_blocks_prod)

                    If (lfield .and. ((field_type .eq. STEELE_SLIT_FINITEX) .OR. (field_type .eq. HARD_SLIT_FINITEX))) Then

                        ! Calculate final pressure (configurational part) in unit of [K/A^3]
                        virialpnavg(1,ibin,ibox) = -1.0d0/(l_avgx*box(2,ibox))*virialpnavg(1,ibin,ibox)
                        virialpnavg(2,ibin,ibox) = -1.0d0/(l_avgx*box(2,ibox))*virialpnavg(2,ibin,ibox)

                        virialptavg(1,1,ibin,ibox) = -0.5d0/(l_avgx*box(2,ibox))*virialptavg(1,1,ibin,ibox)
                        virialptavg(2,1,ibin,ibox) = -0.5d0/(l_avgx*box(2,ibox))*virialptavg(2,1,ibin,ibox)

                        virialptavg(1,2,ibin,ibox) = -0.5d0/(l_avgx*box(2,ibox))*virialptavg(1,2,ibin,ibox)
                        virialptavg(2,2,ibin,ibox) = -0.5d0/(l_avgx*box(2,ibox))*virialptavg(2,2,ibin,ibox)

                        virialptavg(1,3,ibin,ibox) = -0.5d0/(l_avgx*box(2,ibox))*virialptavg(1,3,ibin,ibox)
                        virialptavg(2,3,ibin,ibox) = -0.5d0/(l_avgx*box(2,ibox))*virialptavg(2,3,ibin,ibox)

                        virialptavg(1,4,ibin,ibox) = -0.5d0/(l_avgx*box(2,ibox))*virialptavg(1,4,ibin,ibox)
                        virialptavg(2,4,ibin,ibox) = -0.5d0/(l_avgx*box(2,ibox))*virialptavg(2,4,ibin,ibox)

                    Else 
                        ! Calculate final pressure (configurational part) in unit of [K/A^3]
                        virialpnavg(1,ibin,ibox) = -1.0d0/(box(1,ibox)*box(2,ibox))*virialpnavg(1,ibin,ibox)
                        virialpnavg(2,ibin,ibox) = -1.0d0/(box(1,ibox)*box(2,ibox))*virialpnavg(2,ibin,ibox)
                        virialptavg(1,1,ibin,ibox) = -0.5d0/(box(1,ibox)*box(2,ibox))*virialptavg(1,1,ibin,ibox)                        
                        virialptavg(2,1,ibin,ibox) = -0.5d0/(box(1,ibox)*box(2,ibox))*virialptavg(2,1,ibin,ibox)
                        virialptavg(1,2,ibin,ibox) = -0.5d0/(box(1,ibox)*box(2,ibox))*virialptavg(1,2,ibin,ibox)                        
                        virialptavg(2,2,ibin,ibox) = -0.5d0/(box(1,ibox)*box(2,ibox))*virialptavg(2,2,ibin,ibox)                            
                    End If

                    

                    ! Calculate totoal pressure tensor in unit of [K/A^3]
                    virialpnavg(3,ibin,ibox) = virialpnavg(1,ibin,ibox) + virialpnavg(2,ibin,ibox)
                    virialptavg(3,1,ibin,ibox) = virialptavg(1,1,ibin,ibox) + virialptavg(2,1,ibin,ibox)
                    virialptavg(3,2,ibin,ibox) = virialptavg(1,2,ibin,ibox) + virialptavg(2,2,ibin,ibox)
                    virialptavg(3,3,ibin,ibox) = virialptavg(1,3,ibin,ibox) + virialptavg(2,3,ibin,ibox)
                    virialptavg(3,4,ibin,ibox) = virialptavg(1,4,ibin,ibox) + virialptavg(2,4,ibin,ibox)
                    
                End Do

                ! Caculate kinetic part contribution to the pressure
                ! Loop over bins
                Do ibin = 1, virialpress_slit_bins
                    ! Loop over molecule types
                    Do itype = 1, n_mol_types
                        ! Pressure (Kinetic part) in unit of [K/A^3]
                        pkin(ibin) = pkin(ibin) + (zdenavg(ibin,itype)/mol_mass(itype))*Na*1.0d-24*temp
                    End Do
                End Do

                ! IK route
                If ((virialpress_ctype .EQ. 1) .OR. (virialpress_ctype .EQ. 3) .OR. (virialpress_ctype .EQ. 6)) Then
                    ! Open file
                    OPEN(UNIT=FILE_VIRIALPRESS_SLIT_1,FILE='press_slit_1.txt',STATUS='UNKNOWN',ACCESS='SEQUENTIAL',ACTION='WRITE')

                    ! Write file head
                    Write(FILE_VIRIALPRESS_SLIT_1,*) 'Pressure tensor from virial route using Irving-Kirkwood definition'
                    Write(FILE_VIRIALPRESS_SLIT_1,*) 'Unit: pressure in [bar] and length in [Angstrom]'
                    Write(FILE_VIRIALPRESS_SLIT_1,'(A)') &
                                    & '  z              Pkin               Pt(fw)             Pt(ff)             Pt(tot)        &
                                    &    Pn(fw)             Pn(ff)             Pn(tot)'

                    ! Write data to file
                    ! Loop over bins
                    Do ibin = 1, virialpress_slit_bins

                        If (field_type .eq. HARD_SLIT_FINITEX) Then

                            ! Write to file
                            Write(FILE_VIRIALPRESS_SLIT_1,'(F8.4,1X,F14.4,3X,F16.4,3X,F16.4,3X,F16.4,3X,F16.4,3X,F16.4,3X,F16.4)')&
                                & (DBLE(ibin)-0.5d0)*virialpress_slit_dz+wall_radius, pkin(ibin)*PCOEFF, &
                                & virialptavg(1,1,ibin,ibox)*PCOEFF,virialptavg(2,1,ibin,ibox)*PCOEFF, &
                                & (virialptavg(3,1,ibin,ibox)+pkin(ibin))*PCOEFF, virialpnavg(1,ibin,ibox)*PCOEFF, &
                                & virialpnavg(2,ibin,ibox)*PCOEFF, &
                                & (virialpnavg(3,ibin,ibox)+pkin(ibin))*PCOEFF

                        Else 
                            ! Write to file
                            Write(FILE_VIRIALPRESS_SLIT_1,'(F8.4,1X,F14.4,3X,F16.4,3X,F16.4,3X,F16.4,3X,F16.4,3X,F16.4,3X,F16.4)')&
                                & (DBLE(ibin)-0.5d0)*virialpress_slit_dz+steele_position(1), pkin(ibin)*PCOEFF, &
                                & virialptavg(1,1,ibin,ibox)*PCOEFF, virialptavg(2,1,ibin,ibox)*PCOEFF, &
                                & (virialptavg(3,1,ibin,ibox)+pkin(ibin))*PCOEFF, virialpnavg(1,ibin,ibox)*PCOEFF, &
                                & virialpnavg(2,ibin,ibox)*PCOEFF, &
                                & (virialpnavg(3,ibin,ibox)+pkin(ibin))*PCOEFF

                        End if

                        

                    End Do
                Endif

                ! Harasima route
                If ((virialpress_ctype .EQ. 2) .OR. (virialpress_ctype .EQ. 3) .OR. (virialpress_ctype .EQ. 6)) Then

                    ! Open file
                    OPEN(UNIT=FILE_VIRIALPRESS_SLIT_2,FILE='press_slit_2.txt', &
                         & STATUS='UNKNOWN',ACCESS='SEQUENTIAL',ACTION='WRITE')

                    ! Write file head
                    Write(FILE_VIRIALPRESS_SLIT_2,*) 'Pressure tensor from virial route using Harasima definition'
                    Write(FILE_VIRIALPRESS_SLIT_2,*) 'Unit: pressure in [bar] and length in [Angstrom]'
                    Write(FILE_VIRIALPRESS_SLIT_2,'(A)') &
                                    & '  z              Pkin               Pt(fw)             Pt(ff)             Pt(tot)        &
                                    &    Pn(fw)             Pn(ff)             Pn(tot)'

                    ! Write data to file
                    ! Loop over bins
                    Do ibin = 1, virialpress_slit_bins

                        If (field_type .eq. HARD_SLIT_FINITEX) Then

                            ! Write to file
                            Write(FILE_VIRIALPRESS_SLIT_2, &
                                & '(F8.4,1X,F14.4,3X,F16.4,3X,F16.4,3X,F16.4,3X,F16.4,3X,F16.4,3X,F16.4)') &
                                & (DBLE(ibin)-0.5d0)*virialpress_slit_dz+wall_radius, pkin(ibin)*PCOEFF, &
                                & virialptavg(1,2,ibin,ibox)*PCOEFF,virialptavg(2,2,ibin,ibox)*PCOEFF, &
                                & (virialptavg(3,2,ibin,ibox)+pkin(ibin))*PCOEFF, virialpnavg(1,ibin,ibox)*PCOEFF, &
                                & virialpnavg(2,ibin,ibox)*PCOEFF, &
                                & (virialpnavg(3,ibin,ibox)+pkin(ibin))*PCOEFF

                        Else 
                            ! Write to file
                            Write(FILE_VIRIALPRESS_SLIT_2, &
                                & '(F8.4,1X,F14.4,3X,F16.4,3X,F16.4,3X,F16.4,3X,F16.4,3X,F16.4,3X,F16.4)')  &
                                & (DBLE(ibin)-0.5d0)*virialpress_slit_dz+steele_position(1), pkin(ibin)*PCOEFF, &
                                & virialptavg(1,2,ibin,ibox)*PCOEFF, virialptavg(2,2,ibin,ibox)*PCOEFF, &
                                & (virialptavg(3,2,ibin,ibox)+pkin(ibin))*PCOEFF, virialpnavg(1,ibin,ibox)*PCOEFF, &
                                & virialpnavg(2,ibin,ibox)*PCOEFF, &
                                & (virialpnavg(3,ibin,ibox)+pkin(ibin))*PCOEFF

                        End if

                        

                    End Do
                
                ! End integral contour definition
                End If


                ! Variation of the Harasima type (H-VR)
                If ((virialpress_ctype .EQ. 4) .OR. (virialpress_ctype .EQ. 6)) Then

                    ! Open file
                    OPEN(UNIT=FILE_VIRIALPRESS_SLIT_3,FILE='press_slit_3.txt', &
                         & STATUS='UNKNOWN',ACCESS='SEQUENTIAL',ACTION='WRITE')

                    ! Write file head
                    Write(FILE_VIRIALPRESS_SLIT_3,*) 'Pressure tensor from virial route using H-VR1 definition'
                    Write(FILE_VIRIALPRESS_SLIT_3,*) 'Unit: pressure in [bar] and length in [Angstrom]'
                    Write(FILE_VIRIALPRESS_SLIT_3,'(A)') &
                                    & '  z              Pkin               Pt(fw)             Pt(ff)             Pt(tot)        &
                                    &    Pn(fw)             Pn(ff)             Pn(tot)'

                    ! Write data to file
                    ! Loop over bins
                    Do ibin = 1, virialpress_slit_bins

                        If (field_type .eq. HARD_SLIT_FINITEX) Then

                            ! Write to file
                            Write(FILE_VIRIALPRESS_SLIT_3, &
                                & '(F8.4,1X,F14.4,3X,F16.4,3X,F16.4,3X,F16.4,3X,F16.4,3X,F16.4,3X,F16.4)') &
                                & (DBLE(ibin)-0.5d0)*virialpress_slit_dz+wall_radius, pkin(ibin)*PCOEFF, &
                                & virialptavg(1,3,ibin,ibox)*PCOEFF,virialptavg(2,3,ibin,ibox)*PCOEFF, &
                                & (virialptavg(3,3,ibin,ibox)+pkin(ibin))*PCOEFF, virialpnavg(1,ibin,ibox)*PCOEFF, &
                                & virialpnavg(2,ibin,ibox)*PCOEFF, &
                                & (virialpnavg(3,ibin,ibox)+pkin(ibin))*PCOEFF

                        ! Else 
                        !   ! Write to file
                        !   Write(FILE_VIRIALPRESS_SLIT_H, &
                        !       & '(F8.4,1X,F14.4,3X,F16.4,3X,F16.4,3X,F16.4,3X,F16.4,3X,F16.4,3X,F16.4)')  &
                        !       & (DBLE(ibin)-0.5d0)*virialpress_slit_dz+steele_position(1), pkin(ibin)*PCOEFF, &
                        !       & virialptavg(1,2,ibin,ibox)*PCOEFF, virialptavg(2,2,ibin,ibox)*PCOEFF, &
                        !       & (virialptavg(3,2,ibin,ibox)+pkin(ibin))*PCOEFF, virialpnavg(1,ibin,ibox)*PCOEFF, &
                        !       & virialpnavg(2,ibin,ibox)*PCOEFF, &
                        !       & (virialpnavg(3,ibin,ibox)+pkin(ibin))*PCOEFF

                        End if

                        

                    End Do
                
                ! End integral contour definition
                End If

                ! IK-VR1: variation of the IK type of the first kind 
                If ((virialpress_ctype .EQ. 5) .OR. (virialpress_ctype .EQ. 6)) Then
                    ! Open file
                    OPEN(UNIT=FILE_VIRIALPRESS_SLIT_4,FILE='press_slit_4.txt',STATUS='UNKNOWN',ACCESS='SEQUENTIAL',ACTION='WRITE')

                    ! Write file head
                    Write(FILE_VIRIALPRESS_SLIT_4,*) 'Pressure tensor from virial route using IK-VR1 definition'
                    Write(FILE_VIRIALPRESS_SLIT_4,*) 'Unit: pressure in [bar] and length in [Angstrom]'
                    Write(FILE_VIRIALPRESS_SLIT_4,'(A)') &
                                    & '  z              Pkin               Pt(fw)             Pt(ff)             Pt(tot)        &
                                    &    Pn(fw)             Pn(ff)             Pn(tot)'

                    ! Write data to file
                    ! Loop over bins
                    Do ibin = 1, virialpress_slit_bins

                        If (field_type .eq. HARD_SLIT_FINITEX) Then

                            ! Write to file
                            Write(FILE_VIRIALPRESS_SLIT_4,'(F8.4,1X,F14.4,3X,F16.4,3X,F16.4,3X,F16.4,3X,F16.4,3X,F16.4,3X,F16.4)')&
                                & (DBLE(ibin)-0.5d0)*virialpress_slit_dz+wall_radius, pkin(ibin)*PCOEFF, &
                                & virialptavg(1,4,ibin,ibox)*PCOEFF,virialptavg(2,4,ibin,ibox)*PCOEFF, &
                                & (virialptavg(3,4,ibin,ibox)+pkin(ibin))*PCOEFF, virialpnavg(1,ibin,ibox)*PCOEFF, &
                                & virialpnavg(2,ibin,ibox)*PCOEFF, &
                                & (virialpnavg(3,ibin,ibox)+pkin(ibin))*PCOEFF


                        End if

                        

                    End Do
                Endif
                

            ! End virial pressure tensor    
            End If
            ! -------------------------------------------
            
            ! Surface excess 
            If (lsurfex) Then

                ! Open file
                OPEN(UNIT=FILE_SURFEX,FILE='surface_excess.txt',STATUS='UNKNOWN',ACCESS='SEQUENTIAL',ACTION='WRITE')
                
                ! Initialize surface excess variables (optional sampling)
                surfexavg = 0.0d0
                surfexstd = 0.0d0

                ! Loop over blocks
                Do iblk = n_blocks_equil+1, n_blocks_tot

                    ! Loop over molecule types
                    Do itype = 1, n_mol_types

                        ! Exclude external structure
                        If (initstyle(itype,1) .eq. 'coords') CYCLE

                        ! Update sum
                        surfexavg(itype) = surfexavg(itype) + avg_surfex(itype,iblk)
                        
                    End Do
                    
                End Do

                ! Average
                Do itype = 1, n_mol_types

                    ! Exclude external structure
                    If (initstyle(itype,1) .eq. 'coords') CYCLE

                    surfexavg(itype) = surfexavg(itype)/DBLE(n_blocks_prod)
                End Do

                ! Start Computing standard deviation
                ! Loop over blocks
                Do iblk = n_blocks_equil+1, n_blocks_tot
                
                    ! Loop over molecule types
                    Do itype = 1, n_mol_types

                        ! Exclude external structure
                        If (initstyle(itype,1) .eq. 'coords') CYCLE

                        ! Update sum
                        surfexstd(itype) = surfexstd(itype) + (avg_surfex(itype,iblk) - surfexavg(itype))**2
                    End Do
                End Do


                Write(FILE_SURFEX,*) ' Units     Type        Average          Standard Deviation'

                Do itype = 1, n_mol_types

                    ! Exclude external structure
                    If (initstyle(itype,1) .eq. 'coords') CYCLE

                    surfexstd(itype) = dSQRT(surfexstd(itype)/DBLE(n_blocks_prod))

                    Write(FILE_SURFEX,'(A,3X,I2,7X,F15.7,3X,E15.7)') &
                    &'mmol/m^2', itype, 1.0d23/Na*surfexavg(itype), 1.0d23/Na*surfexstd(itype)

                End Do


                ! Close file
                Close(UNIT=FILE_SURFEX)

            ! End of surface excess
            End If

            ! Lattice constant
            If (llattconst) Then

                ! Open file
                OPEN(UNIT=FILE_LATTCONST,FILE='lattconst.txt',STATUS='UNKNOWN',ACCESS='SEQUENTIAL',ACTION='WRITE')
                
                ! Initialize surface excess variables (optional sampling)
                lattconst = 0.0d0

                ! Loop over blocks
                Do iblk = n_blocks_equil+1, n_blocks_tot

                    ! Loop over molecule types
                    Do inum = 1, lattconst_n

                        ! Update sum
                        lattconst(inum) = lattconst(inum) + avg_lattconst(inum,iblk)
                        
                    End Do
                    
                End Do

                ! Average
                Do inum = 1, lattconst_n
                    lattconst(inum) = lattconst(inum)/DBLE(n_blocks_prod)
                End Do

                Write(FILE_LATTCONST,*) 'Layer          Range              Average'

                Do inum = 1, lattconst_n

                    Write(FILE_LATTCONST,'(I2,7X,A,F7.3,A,F7.3,A,3X,F15.7)') &
                    & inum,'(',lattconst_cut(1,inum),',',lattconst_cut(2,inum),')',lattconst(inum)

                End Do


                ! Close file
                Close(UNIT=FILE_LATTCONST)

            ! End of lattice constant
            End If


            ! Isosteric heat of adsorption
            If (lqst) Then              

                ! Assume 1 box
                ibox = 1 

                ! Initialize
                qstavg = 0.0d0

                ! Open file
                OPEN(UNIT=FILE_QST,FILE='dump_qst.txt',STATUS='UNKNOWN',ACCESS='SEQUENTIAL',ACTION='WRITE')

                ! Loop over blocks
                Do iblk = n_blocks_equil+1, n_blocks_tot

                    ! Loop over molecule types
                    Do itype = 1, n_mol_types

                        ! Exclude external structure
                        If (initstyle(itype,1) .eq. 'coords') CYCLE

                        ! Update sum 
                        ! <N>
                        qstavg(1) = qstavg(1) + avg_qst(1,iblk,ibox)
                        ! <U>
                        qstavg(2) = qstavg(2) + avg_qst(2,iblk,ibox)
                        ! <N^2>
                        qstavg(3) = qstavg(3) + avg_qst(3,iblk,ibox)
                        ! <U*N>
                        qstavg(4) = qstavg(4) + avg_qst(4,iblk,ibox)
                        
                    End Do
                    
                End Do

                ! Average
                Do itype = 1, n_mol_types

                    ! Exclude external structure
                    If (initstyle(itype,1) .eq. 'coords') CYCLE

                    qstavg(1) = qstavg(1)/DBLE(n_blocks_prod)
                    qstavg(2) = qstavg(2)/DBLE(n_blocks_prod)
                    qstavg(3) = qstavg(3)/DBLE(n_blocks_prod)
                    qstavg(4) = qstavg(4)/DBLE(n_blocks_prod)

                End Do

                ! Calculate the isosteric heat of adsorption in J/mol
                fluc_U = qstavg(1)*qstavg(2) - qstavg(4)
                fluc_N = qstavg(3) - qstavg(1)**2
                qst = R * (temp + fluc_U/fluc_N)

                ! Convert to kJ/mol
                qst = qst*1d-3

                Write(FILE_QST,*) ' Units     Type        Average '

                Do itype = 1, n_mol_types

                    ! Exclude external structure
                    If (initstyle(itype,1) .eq. 'coords') CYCLE

                    Write(FILE_QST,'(A,5X,I2,7X,F15.7)') &
                    &'kJ/mol', itype, qst

                End Do


                ! Close file
                Close(UNIT=FILE_QST)

            ! End of surface excess
            End If





            
        Case default
            
      
      End Select
      
      Return
      
      End Subroutine 















