! =========================================================
!  A Monte Carlo simulation program by Kaihang Shi     
!  Prof. Keith Gubbins' research group                 
!  Department of Chemical & Biomolecular Engineering.  
!  North Carolina State University, NC, USA.           
!                                                      
!  Created on Dec. 4th, 2016                 
! ==========================================================


      Program MC

      Use global

      IMPLICIT NONE

      ! Local
      Integer :: ibox, iblock, istep, nstep
      Double Precision :: cpu_start, cpu_end
      Double Precision :: engtst, engerr
      ! For retriving system time
      Integer, Dimension(8) :: date_time
      Character(len=10), Dimension(3) :: time_value

      ! Initialize local variables
      nstep = 0

      ! Start CPU TIME COUNTING
      Call CPU_TIME(cpu_start)

      ! Print initial banner
      Write(*,'(A)') '====================================================='
      Write(*,'(A)') '                                                     '
      Write(*,'(A)') ' Monte Carlo Simulation Package v2.2, by Kaihang Shi'
      Write(*,'(A)') ' Department of Chemical & Biomolecular Engineering'
      Write(*,'(A)') ' North Carolina State University, Raleigh, NC, USA'
      Write(*,'(A)') ' Created on Dec. 4th, 2016'
      Write(*,'(A)') ' Last modified in June 2020'
      Write(*,'(A)') '                                                     '
      Write(*,'(A)') ' ===================================================='

      ! Initialize the random number generator
      Call seed_random(idum)
      Write(*,'(A,I0)') 'Random number seed: ', idum

      ! Read the input file
      Call read_input

      ! Calculate energy and print initial configuration (iblock = 0)
      Do ibox = 1, n_box
        Call eng_total(ibox,energy(ibox))
        Call print_xyz(ibox,nstep,'sites')

        ! Check if overlap
        If (OVERLAP) Then
            Write(*,'(A,I0)') 'FATAL ERROR: INITIAL CONFIGURATION OVERLAPS FOR BOX: ', ibox
            STOP
        End If

      End Do

      ! Start simulation
      Write(*,*) 'Starting Simulation ...'

      ! Relax system first
      Call relax

      ! Loop over the blocks
      Do iblock = 1, n_blocks_tot

        ! Flush I/O unit (write information in buffer immediately to file)
        ! unit=6 is screen
        Call FLUSH(6)

        ! Reset basic statistics 
        Call stats(1,0,0)

        ! Reset sampling statistics
        ! Only perform during the production stage
        if(lsampling .and. (iblock .gt. n_blocks_equil)) Call sample(1,0,0)

        ! Manual option to adjust probability of each type of MC move
!       Call update_prob(iblock)

        ! Gradually change chemical potential for uVT for dense system (add on 8-13-2017)
        ! First 1000 blocks gradually changing chempot
!       if(iblock .le. 1000) Call update_chempot(iblock)

        ! Loop over the steps per block
        Do istep = 1, block_size

            ! Update nstep 
            nstep = nstep + 1

            ! Recalculate and check the total and updated energy consistency for each box 
            ! Modified on 6-9-2020
            IF(MOD(nstep,check_freq) .eq. 0) THEN
                DO ibox = 1,n_box
                    Call eng_total(ibox,engtst)
                    engerr = dABS((engtst - energy(ibox)))
                    
                    IF (engerr .GE. 1.0d-8) THEN
                        Write(*,*) '!---------------------------------------------------------------!'
                        WRITE(*,*) '!   WARNING: LARGE ERROR BETWEEN THE TOTAL AND UPDATED ENERGY   !'
                        Write(*,*) '!---------------------------------------------------------------!'
                        WRITE(*,'(A6,I0,A8,I0,A6,I0,A10,F20.8,A10,F20.8)')&
                        &'STEP: ',nstep,' BLOCK: ',iblock,' BOX: ',ibox,' UPD ENG: ',energy(ibox),' TOT ENG: ',engtst
                    ENDIF

                    ! Assign accurate energy to the global energy variable
                    energy(ibox) = engtst

                ENDDO
            ENDIF

            ! Perform a trial move
            Call trial_move(iblock, istep)

            ! Accumulate current block statistics
            Call stats(2,iblock,istep)

            ! Update max values (max_trans, max_rotat etc)
            If (MOD(nstep,10) .eq. 0) Call update_max

            ! Do optional sampling
            ! Only perform during the production stage
            if(lsampling .and. (iblock .gt. n_blocks_equil)) Call sample(2,iblock,istep)

            ! Perform volume change to calculate pressure tensor (Thermal route)
!           if(lthermopress_slit .and. (iblock .gt. n_blocks_equil)) Call press_thermo_slit(1,iblock,istep)

            ! Calculate pressure tensor for slit pore geometry (virial route)
            if(lvirialpress_slit .and. (iblock .gt. n_blocks_equil)) Call press_virial_slit(1,iblock,istep)

            ! Calculate pressure tensor for cylindrical system (virial route using IK and/or Harasima definition)
            if(lvirialpress_cylin .and. (iblock .gt. n_blocks_equil)) Call press_virial_cylinder(1,iblock,istep)

            ! Calculate pressure tensor for planar surface (virial route using IK and/or Harasima definition)
            if(lvirialpress .and. (iblock .gt. n_blocks_equil)) Call press_virial(1,iblock,istep)

            ! Perform Widom Insertion Method
            ! Only perform during the production stage
            if(lwdm .and. (iblock .gt. n_blocks_equil)) Call widom(1,iblock,istep)

            ! Dump molecules' coordinates to file for post process
            If (ldumpxyz .and. (dumpxyz_freq .ne. 0) .and. &
                    & (MOD(nstep,dumpxyz_freq) .eq. 0) .and. (iblock .gt. n_blocks_equil)) Then 
                Do ibox = 1, n_box
                    If (dump_moltype .ne. 0) Then
                        Call print_xyz(ibox,nstep,'mol')
                    Else 
                        Call print_xyz(ibox,nstep,'all')
                    End If
                    
                End Do
            End if

            ! Write restart file 
            If (lwriterestart .and. (rst_freq .ne. 0) .and. (MOD(nstep,rst_freq) .eq. 0)) Then
                Do ibox = 1, n_box
                    Call print_xyz(ibox,nstep,'restart')
                End Do
            End If 



        ! End of the block run
        End Do

        ! Average statistics for iblock
        Call stats(3,iblock,0)

        ! Average optional sampling statistics for iblock
        if(lsampling .and. (iblock .gt. n_blocks_equil)) Call sample(3,iblock,0)

        ! Average pressure tensor from virial route for iblock in slit pore geometry
        if(lvirialpress_slit .and. (iblock .gt. n_blocks_equil)) Call press_virial_slit(2,iblock,0)

        ! Average pressure tensor from virial route for iblock in cylindrical coordinates
        if(lvirialpress_cylin .and. (iblock .gt. n_blocks_equil)) Call press_virial_cylinder(2,iblock,0)

        ! Average pressure tensor from virial route for iblock for planar surface version
        if(lvirialpress .and. (iblock .gt. n_blocks_equil)) Call press_virial(2,iblock,0)

        ! Average Widom statistics for iblock
        if(lwdm .and. (iblock .gt. n_blocks_equil)) Call widom(2,iblock,0)

        ! Dump statistics to screen
        Call stats(4,iblock,0)

        

      ! End loop over blocks
      ! Finish simulation runs  
      End Do

      ! Calculate total energy for the final configuration
      ! & print final configurations
      Do ibox = 1, n_box
        Call eng_total(ibox,energy(ibox))
        Call print_xyz(ibox,nstep,'sites')
      End Do

      ! Calculate & Dump final results to screen
      Call stats(5,0,0)

      ! Write optional sampling statistics to files
      ! Including the pressure tensor
      If (lsampling) Call sample(4,0,0)


      ! End CPU TIME COUNTING
      Call CPU_TIME(cpu_end)
      Write(*,*) ' '
      Write(*,'(A,E15.7,A,4X,F7.2,A)') 'Total CPU TIME: ', (cpu_end - cpu_start), ' s', &
                                        & (cpu_end - cpu_start)/3600.0d0, ' hrs'

      ! Print time
      Call date_and_time(time_value(1),time_value(2),time_value(3),date_time)
      Write(*,'(A,I0,A,I0,4X,I0,A,I0,A,I0)') 'Date: ', &
        & date_time(5),':',date_time(6), date_time(2),'/',date_time(3),'/',date_time(1)

      End Program MC





