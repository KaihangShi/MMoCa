! ==========================================================
! This subroutine is used to allocate global variables declared 
! in global.f90.
! Created on Dec. 4th, 2016 by Kaihang Shi
! Last modified in July, 2019
! ==========================================================



      Subroutine global_allocate(selector)

      Use global

      IMPLICIT NONE

      ! Passed
      Character(Len=*) selector

      ! Local
	  Integer :: ierr




	  ! Initialize ierr, '0' represents successful 
	  ierr = 0

	  ! Select array group to allocate
	  Select Case (selector)

	  	! Box dimension
	  	Case('Box')
	  		!Allocate box dimensions
	  		Allocate(box(3,n_box), STAT=ierr)
	  		Allocate(vol(n_box), STAT=ierr)
	  		Allocate(vol_pore(n_box), STAT=ierr)
	  		Allocate(linit(n_box), STAT=ierr)

	  		! Check allocation error
	  		If(ierr .ne. 0) Then
	  			Write(*,*) 'ALLOCATION ERROR FOR BOX DIMENSION IN global_allocate.f90'
	  			STOP
	  		End If

	  	! System arrays
	  	Case('System')

	  		! Molecular topology
	  		Allocate(n_mol(n_mol_types,n_box), STAT=ierr)
		    Allocate(n_mol_tot(n_box), STAT=ierr)
		    Allocate(n_sites(n_mol_types), STAT=ierr)
		    Allocate(n_site_types(n_mol_types), STAT=ierr)
		    Allocate(n_sites_disp(n_mol_types), STAT=ierr)
		    Allocate(n_sites_elec(n_mol_types), STAT=ierr)
		    Allocate(mol_type(2*n_mol_max,n_box), STAT=ierr)
		    Allocate(mol_type_name(n_mol_types), STAT=ierr)
		    Allocate(site_type(n_sites_max,n_mol_types), STAT =ierr)
		    Allocate(site_type_name(n_nonbond_max,n_mol_types), STAT=ierr)  
		    Allocate(disp_sites(n_sites_max,n_mol_types),STAT=ierr)
		    Allocate(elec_sites(n_sites_max,n_mol_types), STAT=ierr)
		    Allocate(initstyle(n_mol_types,n_box), STAT=ierr)
		    Allocate(mol_mass(n_mol_types), STAT=ierr)
		    Allocate(lambda(n_mol_types), STAT=ierr)
		    Allocate(mol_act(n_mol_types), STAT=ierr)
		    Allocate(rx(2*n_mol_max,n_box), STAT=ierr)
		    Allocate(ry(2*n_mol_max,n_box), STAT=ierr)
		    Allocate(rz(2*n_mol_max,n_box), STAT=ierr)
		    Allocate(q1(2*n_mol_max,n_box), STAT=ierr)
		    Allocate(q2(2*n_mol_max,n_box), STAT=ierr)
		    Allocate(q3(2*n_mol_max,n_box), STAT=ierr)
		    Allocate(q4(2*n_mol_max,n_box), STAT=ierr)
		    Allocate(rx_s(n_sites_max,2*n_mol_max,n_box), STAT=ierr)
		    Allocate(ry_s(n_sites_max,2*n_mol_max,n_box), STAT=ierr)
		    Allocate(rz_s(n_sites_max,2*n_mol_max,n_box), STAT=ierr)
		    Allocate(rx_i(n_sites_max,n_mol_types), STAT=ierr)
		    Allocate(ry_i(n_sites_max,n_mol_types), STAT=ierr)
		    Allocate(rz_i(n_sites_max,n_mol_types), STAT=ierr)
		    Allocate(sites_name(n_sites_max,n_mol_types), STAT=ierr)

		    ! External structure flag
		    Allocate(ext_struc(n_box), STAT=ierr)

		    ! External field
		    Allocate(energy_field(n_box), STAT=ierr)
		    ! 10-4-3 Steele potential
		    Allocate(steele_position(2), STAT=ierr)
		    Allocate(steele_site_name(n_nonbond_max), STAT=ierr)
		    Allocate(steele_sigmasf(n_nonbond_max,n_mol_types_max), STAT=ierr)
		    Allocate(steele_sigmasfsq(n_nonbond_max,n_mol_types_max), STAT=ierr)
		    Allocate(steele_epsilonsf(n_nonbond_max,n_mol_types_max), STAT=ierr)
		    Allocate(steele_posx(2), STAT=ierr)
		    Allocate(steele_avgx(2), STAT=ierr)
		    ! Coarse-grained wall potential
		    Allocate(cg_wall_position(2), STAT=ierr)
		    ! Maximum 10,000 points
		    Allocate(cg_wall_egrid(10000), STAT=ierr)
		    Allocate(cg_wall_dgrid(10000), STAT=ierr)
		    Allocate(cg_ff_egrid(10000), STAT=ierr)
		    Allocate(cg_ff_dgrid(10000), STAT=ierr)
		    Allocate(eff_sig_z(10000), STAT=ierr)
		    Allocate(eff_sig(10000), STAT=ierr)
		    ! cg_wall_strc
		    Allocate(cg_strc_dgrid(10000), STAT=ierr)
		    Allocate(cg_strc_dacc(10000), STAT=ierr)
		    ! For modified cg potential
		    Allocate(cg_wall1_egrid(10000), STAT=ierr)
		    Allocate(cg_wall1_dgrid(10000), STAT=ierr)




		    ! Potential parameters
		    ! Note: epsilon(1,1,2,2) means epsilon between site 1 on molecule type 1 and 
		    ! site 2 on moelcule type 2
		    Allocate(epsilon(n_nonbond_max,n_mol_types,n_nonbond_max,n_nonbond_max), STAT=ierr)
		    Allocate(sigma(n_nonbond_max,n_mol_types,n_nonbond_max,n_nonbond_max), STAT=ierr)
		    Allocate(sigmasq(n_nonbond_max,n_mol_types,n_nonbond_max,n_mol_types), STAT=ierr)

			! Ewald sum parameters
		    Allocate(k_max(3,n_box), STAT=ierr)
		    Allocate(q(n_nonbond_max,n_mol_types,n_nonbond_max,n_mol_types), STAT=ierr)
		    Allocate(qsq(n_nonbond_max,n_mol_types,n_nonbond_max,n_mol_types), STAT=ierr)
		    Allocate(qex(2,n_nonbond_max,n_mol_types,n_nonbond_max,n_mol_types), STAT=ierr)
		    Allocate(qsqex(2,n_nonbond_max,n_mol_types,n_nonbond_max,n_mol_types), STAT=ierr)
		   	ALLOCATE(maxk(n_box), STAT=ierr)
	        ALLOCATE(ksq_max(n_box), STAT=ierr)
	        ALLOCATE(ewld_self(n_mol_types), STAT=ierr)
	        Allocate(ewld_intra(n_mol_types), STAT=ierr)
	        Allocate(ewld_real(n_box), STAT=ierr)
	        Allocate(ewld_fourier(n_box), STAT=ierr)
	        Allocate(ewld_tot(n_box), STAT=ierr)
	        Allocate(ewld_slab(n_box), STAT=ierr)
	        ! Initialize
	        ewld_self = 0.0d0
	        ewld_intra = 0.0d0

		    ! Thermodynamic properties
		    Allocate(mu(n_mol_types), STAT=ierr)
		    Allocate(dmu(n_mol_types), STAT=ierr)
		    Allocate(energy(n_box), STAT=ierr)
		    Allocate(energy_disp(n_box), STAT=ierr)

		    ! Tail Correction
		    Allocate(energy_tail(n_box), STAT=ierr)
		    Allocate(vdw_tail(n_mol_types,n_mol_types), STAT=ierr)

		    ! MC move statistical variables
		    Allocate(trans_prob(n_mol_types), STAT=ierr)
		    Allocate(rotat_prob(n_mol_types), STAT=ierr)
		    Allocate(transfer_prob(n_mol_types), STAT=ierr)
		    Allocate(max_trans(n_mol_types,n_box), STAT=ierr)
		    Allocate(max_rotat(n_mol_types,n_box), STAT=ierr)
		    Allocate(max_vol(n_box), STAT=ierr)
		    Allocate(in_stat(2,n_mol_types), STAT=ierr)
		    Allocate(rem_stat(2,n_mol_types), STAT=ierr)
		    Allocate(vol_stat(2,n_box), STAT=ierr)
		    Allocate(swap_ident_stat(2,n_box), STAT=ierr)
		    Allocate(trans_stat(2,n_mol_types,n_box), STAT=ierr)
		    Allocate(rotat_stat(2,n_mol_types,n_box), STAT=ierr)
		    Allocate(transfer_stat(2,n_mol_types,n_box), STAT=ierr)

		    ! Widom insertion 
		    Allocate(muid(n_mol_types), STAT=ierr)
		    Allocate(widom_stat(n_mol_types,n_box), STAT=ierr)
		    Allocate(widom_freq(n_mol_types,n_box), STAT=ierr)

		    ! Surface excess (Sampling options)
		    Allocate(surfex_bulk(n_mol_types), STAT=ierr)
		    Allocate(surfex_n_mol(n_mol_types), STAT=ierr)


		    ! Check for allocation error
	  		If(ierr .ne. 0) Then
	  			Write(*,*) 'ALLOCATION ERROR FOR SYSTEM ARRAYS IN global_allocate.f90'
	  			STOP
	  		End If


	  	! Allocate block average array
	  	! Allocatable variables need block info
	  	Case('Block')

			Allocate(blk_nmol(block_size,n_mol_types,n_box), STAT=ierr)
			Allocate(blk_rho(block_size,n_mol_types,n_box), STAT=ierr)
			Allocate(blk_eng(block_size,n_box), STAT=ierr)
			Allocate(blk_eng_disp(block_size,n_box), STAT=ierr)
			Allocate(blk_vol(block_size,n_box), STAT=ierr)
			Allocate(avg_nmol(n_blocks_tot,n_mol_types,n_box), STAT=ierr)
			Allocate(avg_rho(n_blocks_tot,n_mol_types,n_box), STAT=ierr)
			Allocate(avg_eng(n_blocks_tot,n_box), STAT=ierr)
			Allocate(avg_eng_disp(n_blocks_tot,n_box), STAT=ierr)
			Allocate(avg_vol(n_blocks_tot,n_box), STAT=ierr)

			! Widom total chemical potential
			Allocate(widom_sample(n_blocks_tot,n_mol_types,n_box), STAT=ierr)

			! Check for allocation error
	  		If(ierr .ne. 0) Then
	  			Write(*,*) 'ALLOCATION ERROR FOR BLOCK ARRAYS IN global_allocate.f90'
	  			STOP
	  		End If


	  	! Allocate sampling variables
	  	Case('Sampling')

	  		! Allocate z-density variables
	  		If (lzdensity) Then
	  			Allocate(blk_zden(zden_bins,n_mol_types,block_size), STAT=ierr)
	  			Allocate(avg_zden(zden_bins,n_mol_types,n_blocks_tot), STAT=ierr)
	  		End If

	  		! Allocate r-density variables
	  		If (lrdensity) Then
	  			Allocate(blk_rden(rden_bins,n_mol_types,block_size), STAT=ierr)
	  			Allocate(avg_rden(rden_bins,n_mol_types,n_blocks_tot), STAT=ierr)
	  		End If

	  		! Allocate surface_excess variables
	  		If (lsurfex) Then
	  			Allocate(blk_surfex(n_mol_types,block_size), STAT=ierr)   
	  			Allocate(avg_surfex(n_mol_types,n_blocks_tot), STAT=ierr) 
	  		End If

	  		! Allocate thermopress variables
	  		If (lthermopress_slit) Then
	  			Allocate(thermopress_slit_stat(n_box), STAT=ierr)
	  			Allocate(thermopress_slit_sample(3,n_blocks_tot,n_box), STAT=ierr)
	  			Allocate(thermopress_slit_pt(thermopress_slit_bins,n_blocks_tot,n_box), STAT=ierr)
	  			Allocate(thermopress_slit_pn(thermopress_slit_bins,n_blocks_tot,n_box), STAT=ierr)
	  		End If

	  		! Allocate virialpress variables for slit geometry
	  		If (lvirialpress_slit) Then
	  			Allocate(virialpress_slit_stat(n_box), STAT=ierr)
	  			! first rank: 2 - denotes two contributions; 1 is fluid-wall, 2 is fluid-fluid
	  			! second rank: 5 - denotes different definitions; 1 is IK, 2 is Harasima, 3 is H-VR1, 4 is IK-VR1 ...
	  			Allocate(virialpress_slit_pt(2,5,virialpress_slit_bins,n_blocks_tot,n_box), STAT=ierr)
	  			Allocate(virialpress_slit_pn(2,virialpress_slit_bins,n_blocks_tot,n_box), STAT=ierr)
	  		End If

	  		! Allocate virialpress variables for planar surface
	  		If (lvirialpress) Then
	  			Allocate(virialpress_stat(n_box), STAT=ierr)
	  			! first rank: 2 - denotes two definitions; 1 is IK, 2 is Harasima
	  			Allocate(virialpress_pt(2,zden_bins,n_blocks_tot,n_box), STAT=ierr)
	  			Allocate(virialpress_pn(zden_bins,n_blocks_tot,n_box), STAT=ierr)
	  		End If

	  		! Allocate virialpress variables for cylindrical coordinates
	  		If (lvirialpress_cylin) Then
	  			Allocate(virialpress_cylin_stat(n_box), STAT=ierr)
	  			! first rank: 2 - denotes two contributions; 1 is fluid-wall, 2 is fluid-fluid
	  			! second rank: 2 - denotes two definitions; 1 is IK, 2 is Harasima
	  			Allocate(virialpress_cylin_ptt(2,2,rden_bins,n_blocks_tot,n_box), STAT=ierr)
	  			Allocate(virialpress_cylin_ptz(2,2,rden_bins,n_blocks_tot,n_box), STAT=ierr)
	  			Allocate(virialpress_cylin_pnr(2,2,rden_bins,n_blocks_tot,n_box), STAT=ierr)
	  			
	  		End If

	  		! Allocate lattice constant sampling
	  		If (llattconst) Then
	  			Allocate(lattconst_stat(lattconst_n), STAT=ierr)
	  			Allocate(blk_lattconst(lattconst_n,block_size), STAT=ierr)
	  			Allocate(avg_lattconst(lattconst_n,n_blocks_tot), STAT=ierr)
	  			Allocate(lattconst(lattconst_n), STAT=ierr)
	  		End If

	  		! Allocate variables for isosteric heat of adsorption
	  		If (lqst) Then
	  			Allocate(blk_qst(4,block_size,n_box), STAT=ierr)
	  			Allocate(avg_qst(4,n_blocks_tot,n_box), STAT=ierr)
	  		End If

	  		! Allocate variables for virial
	  		If (ldumpvir) Then
	  			Allocate(vir(n_box), STAT=ierr)
	  			Allocate(vir_hyp(n_box), STAT=ierr)
	  		End If


	  	! Allocate variables for cell list 
	  	! Added on July 15, 2019
	  	Case('Clist')

  			! Assume max 50*50*50 cells 
  			Allocate(clist_hoc(125000), STAT=ierr)
  			! This array for storing head of chain for substrate sites	  			
  			Allocate(clist_hoc_sub(125000), STAT=ierr)	  			
  			Allocate(clist_llist(2*n_mol_max), STAT=ierr)
  			Allocate(clist_llist_sub(n_sites_max), STAT=ierr)
  			! second rank 27 means only loop over the 27 neighbor cells for 3D (including self)
  			Allocate(clist_neigh(125000,27), STAT=ierr)
  			Allocate(clist_loca(50,50,50), STAT=ierr)
	  			

	  	! Default
	  	Case Default

	  		Write(*,*) 'FATAL ERROR: INVALID SELECTOR IN global_allocate.F90'
	  		STOP


	  End Select

	  Return


	  End Subroutine global_allocate







		 











































