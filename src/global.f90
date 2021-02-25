! ==========================================================
! This module is used to set global variables
! Created on Dec. 2, 2016 by Kaihang Shi
! ==========================================================


	  Module global

	  IMPLICIT NONE

	  ! Ensemble 
	  !================================================!
	  Integer, save :: ensmbl


	  ! Simulation box info
	  !================================================!

	  ! Max number of simulation boxes allowed
	  Integer, Parameter :: n_box_max = 3

	  ! Number of boxes
	  Integer, save :: n_box

	  ! Simulation box dimension [Angstrom]
	  Double Precision, Dimension(:,:), Allocatable :: box

	  ! Simulation box volume [Angstrom**3] (Lx*Ly*Lz)
	  Double Precision, Dimension(:), Allocatable :: vol

	  ! Pore volume for slit geometry [Angstrom**3] (Lx*Ly*H) (add 8-18-2017)
	  Double Precision, Dimension(:), Allocatable :: vol_pore

	  ! Flag to determine if read the initial config from external file
	  Logical, Dimension(:), Allocatable :: linit


	  ! Molecular topology
	  !================================================!

	  ! Max number of molecules in each box
	  Integer, Parameter :: n_mol_max = 100000

	  ! Max number of molecule types allowed
	  Integer, Parameter :: n_mol_types_max = 20

	  ! Number of nonbonded sites in the system (equal to Max numebr of site types allowed)
	  Integer, Parameter :: n_nonbond_max = 100

	  ! Number of nonbonded sites in the system 
	  Integer, save :: n_nonbond_types

	  ! Number of molecule types
	  Integer, save :: n_mol_types

	  ! Maximum number of sites on each molecule
	  Integer, save :: n_sites_max

	  ! Number of molecule of each type 
	  Integer, Dimension(:,:), Allocatable :: n_mol

	  ! Total number of molecules in each box
	  Integer, Dimension(:), Allocatable :: n_mol_tot

	  ! Total number of sites on each type of molecule 
	  Integer, Dimension(:), Allocatable :: n_sites

	  ! Number of site types on each molecule
	  Integer, Dimension(:), Allocatable :: n_site_types

	  ! Number of dispersion sites on each type of molecule
	  Integer, Dimension(:), Allocatable :: n_sites_disp

	  ! Number of ion sites on each type of molecule
	  Integer, Dimension(:), Allocatable :: n_sites_elec

	  ! Type of each molecule
	  Integer, Dimension(:,:), Allocatable :: mol_type

	  ! Name of each molecule type 
	  Character(Len=6), Dimension(:), Allocatable :: mol_type_name

	  ! Type of each site
	  Integer, Dimension(:,:), Allocatable :: site_type

	  ! Name of each molecule site
	  Character(Len=6), Dimension(:,:), Allocatable :: site_type_name

	  ! List of dispersion sites on each molecule type
	  Integer, Dimension(:,:), Allocatable :: disp_sites

	  ! List of ion sites on each molecule type
	  Integer, Dimension(:,:), Allocatable :: elec_sites

	  ! Initial style of each molecule type ('simple_cubic', 'coords' etc.)
	  Character(Len=50), Dimension(:,:), Allocatable :: initstyle

	  ! Mass of each molecule type [g/mol]
	  Double Precision, Dimension(:), Allocatable :: mol_mass

	  ! de Broglie wavelength for each molecule type
	  Double Precision, Dimension(:), Allocatable :: lambda

	  ! Activity*volume for each molecule type  z=exp(mu/T)/lambda**3(*vol) 
	  Double Precision, Dimension(:), Allocatable :: mol_act

	  ! Center of mass (com) position of each molecule
	  Double Precision, Dimension(:,:), Allocatable :: rx, ry, rz

	  ! Orientation quaternions of each molecule
	  Double Precision, Dimension(:,:), Allocatable :: q1, q2, q3, q4

	  ! Position of the sites on each molecule 
	  Double Precision, Dimension(:,:,:), Allocatable :: rx_s, ry_s, rz_s

	  ! Internal coordinates of the sites on each molecule (wrt center of mass)
	  Double Precision, Dimension(:,:), Allocatable :: rx_i, ry_i, rz_i

	  ! Site name for all sites in the molecule 
	  Character(Len=6), Dimension(:,:), Allocatable :: sites_name

	  ! External Structure flag 
	  ! Set to .true. if we have molecules loading from coords.in file
	  Logical, Dimension(:), Allocatable :: ext_struc


	  ! Potential parameters 
	  !================================================!

	  Integer, save :: potential
	  Integer, save :: mix_rule
	  Double Precision, Dimension(:,:,:,:), Allocatable :: epsilon, sigma, sigmasq

	  ! Cutoff distance for potential in [Angstrom]
	  Double Precision, save :: r_cut, r_cutsq

	  ! Inner cutoff distance to speed up LJ system in [Angstrom]
	  Double Precision, save :: r_min, r_minsq

	  ! Check overlap of sites (r < r_min)
	  Logical, save :: OVERLAP



	  ! External Field
	  !================================================!

	  ! Flag to turn on/off external field
	  Logical, save :: lfield

	  ! External field types
	  Integer, save :: field_type

	  ! External field energy
	  Double Precision, Dimension(:), Allocatable :: energy_field

	  ! Hard Wall 
	  ! Assume hard wall center is always at z = 0.0d0
	  !-----------------

	  ! Hard Wall thickness (radius)
	  Double Precision, save :: wall_radius

	  ! 10-4-3 Steele Potential 
	  ! Surface alway perpendicular to z
	  !-----------------

	  ! Position of Steele wall [A]
	  ! 1 - down wall; 2 - upper wall
	  Double Precision, Dimension(:), Allocatable :: steele_position

	  ! Cutoff radius for Steele wall [A]
	  Double Precision, save :: steele_cut

	  ! Delta value (the spacing between two adjacent graphene layers) [A]
	  Double Precision, save :: steele_delta

	  ! Rho_s (solid density) [A^-3]
	  Double Precision, save :: steele_rhos

	  ! Interaction Parameters
	  Double Precision, Dimension(:,:), Allocatable :: steele_sigmasf, steele_sigmasfsq, steele_epsilonsf

	  ! Nonbonded type name for Steele potential use 
	  Character(Len=6), Dimension(:), Allocatable :: steele_site_name

	  ! Pore mouth position for the steele_slit_finitex 
	  ! 3-24-2018: The pore length in x-direction is now automatically set to 1/3 of the Lx
	  ! 1 - lower boundary, 2- upper boundary
	  Double Precision, Dimension(:), Allocatable :: steele_posx

	  ! Averaging region position for the steele_slit_finitex
	  ! 3-24-2018: The averaging length in x-direction is now automatically set to 6sig_ff (following Yun Long's method)
	  ! 6-9-2020:  Now the size of the averaing region can be set freely by users
	  Double Precision, Dimension(:), Allocatable :: steele_avgx
	  ! Length of the averaging region in the x-direction 
	  Double Precision, save :: l_avgx

	  ! Coarse-grained FEA solid-fluid potential (cg_wall)
	  ! Surface is always perpendicular to z
	  ! Slit pore model
	  ! ---------------------
	  ! 1 - bottom wall; 2 - top wall
	  ! Wall position is symmetric
	  Double Precision, Dimension(:), Allocatable :: cg_wall_position

	  ! ! Cutoff radius for Steele wall [A]
	  ! Double Precision, save :: cg_wall_cut

	  ! Starting z_position for fluid-wall calculation [A]
	  Double Precision, save :: cg_wall_initz

	  ! Resolution of cg potential [A]
	  Double Precision, save :: cg_wall_res

	  ! Resolution of effective fluid-fluid interaction [A]
	  Double Precision, save :: cg_ff_res

	  ! Potential grid point counter
	  Integer, save :: cg_wall_maxnum

	  ! Potential grid point counter
	  Integer, save :: cg_ff_maxnum

	  ! Cutoff range for piecewise fluid-fluid interaction [A]
	  ! Lower limit and upper limit
	  Double Precision, save :: cg_ff_lob, cg_ff_upb

	  ! Effective sigma for the particle
	  Double Precision, save :: eff_sigma
	  Double Precision, save :: eff_sigmasq

	  ! Coarse-grained potential energy [K]
	  Double Precision, Dimension(:), Allocatable :: cg_wall_egrid
	  ! Vertical distance to the wall [A]
	  Double Precision, Dimension(:), Allocatable :: cg_wall_dgrid

	  ! Set up for modified cg potential use
	  ! Coarse-grained potential energy [K]
	  Double Precision, Dimension(:), Allocatable :: cg_wall1_egrid
	  ! Vertical distance to the wall [A]
	  Double Precision, Dimension(:), Allocatable :: cg_wall1_dgrid

	  ! Coarse-grained effective fluid-fluid energy [K]
	  Double Precision, Dimension(:), Allocatable :: cg_ff_egrid
	  ! Vertical distance to the wall [A]
	  Double Precision, Dimension(:), Allocatable :: cg_ff_dgrid

	  ! z position for effective sigma value [A]
	  Double Precision, Dimension(:), Allocatable :: eff_sig_z
	  ! Effective sigma value [A]
	  Double Precision, Dimension(:), Allocatable :: eff_sig

	  ! Coarse-grained fluid-wall potential (cg_wall_cos)
	  ! Use cosine to mimic the geometric roughness of the surface
	  ! Added on Feb 27 2018
	  ! -------------------------
	  ! Number of periodicity 
	  Double Precision, save :: cg_wall_cosnp
	  ! Periodicity
	  Double Precision, save :: cg_wall_cospe

	  ! Amplitude of the cosine function
	  Double Precision, save :: cg_wall_cosap
	

	  ! For cg_wall_strc
	  ! Use an input file "cg_wall_strc.in" to read in accessible volume for COM
	  ! Added on March 4 2018
	  ! -------------------------
	  ! Accsible length at specified z-direction [A]
	  Double Precision, Dimension(:), Allocatable :: cg_strc_dacc
	  ! Vertical distance to the wall [A]
	  Double Precision, Dimension(:), Allocatable :: cg_strc_dgrid
	  ! Potential grid point counter
	  Integer, save :: cg_strc_maxnum



	  ! Tail Correction (12-6 Lennard Jones potential)
	  ! Applicable to multi-type atoms system
	  !================================================!

	  ! Tail correction flag
	  Logical, save :: ltailc

	  ! Total tail correction energy
	  Double Precision, Dimension(:), Allocatable :: energy_tail

	  ! Tail correction energy for each molecule pair
	  Double Precision, Dimension(:,:), Allocatable :: vdw_tail

	  ! Ewald sum
	  !================================================!

	  ! Flag to turn on/off ewald sum
	  Logical, save :: lewld

	  ! Ewald Style
	  Integer, save :: ewld_style

	  ! Electrostatic parameters
	  Double Precision, Dimension(:,:,:,:), Allocatable :: q, qsq

	  ! Perturbated electrostatic parameters for Explicit Mixing rules
	  ! First rank indicates perturbated parameters for ad-surf or ad-ad
	  ! It can be expanded to multi-interaction mode
	  Double Precision, Dimension(:,:,:,:,:), Allocatable :: qex, qsqex

	  ! Cut-off radius in real space 
	  Double Precision, save :: rcelect, rcelectsq

	  ! Real space convergence parameter
	  Double Precision, save :: kalp, alpha

	  ! Maximum number of reciprocal vector 
	  Integer, Dimension(:), Allocatable :: maxk
	  Integer, Dimension(:,:), Allocatable :: k_max
	  Integer, Dimension(:), Allocatable :: ksq_max

	  ! Reciprocal vector array
	  Double Precision, Dimension(:,:), Allocatable :: k_vec

	  ! Reciprocal component vectors 
	  Complex*16, Dimension(:,:,:,:), Allocatable :: eikx, eiky, eikz

	  ! Reciprocal component vectors for a single molecule
	  Complex*16, Dimension(:,:), Allocatable :: eikx_mol, eiky_mol, eikz_mol

	  ! Reciprocal structure factor
	  Complex*16, Dimension(:,:), Allocatable :: skewld, skewld_new

	  ! Reciprocal structure factor 
	  ! For explicit mixing rule
	  ! First ranking: 1-surface sites, 2-adsorbates sites
	  Complex*16, Dimension(:,:,:), Allocatable :: skewld_ex, skewld_ex_new

	  ! Self Interaction energy
	  Double Precision, Dimension(:), Allocatable :: ewld_self

	  ! Real space intramolecular energy
	  Double Precision, Dimension(:), Allocatable :: ewld_intra

	  ! Real space interaction energy
	  Double Precision, Dimension(:), Allocatable :: ewld_real

	  ! Fourier space energy
	  Double Precision, Dimension(:), Allocatable :: ewld_fourier

	  ! Total Coulomb energy
	  Double Precision, Dimension(:), Allocatable :: ewld_tot

	  ! Flag to turn on/off slab correction
	  Logical, save :: lslabc

	  ! Surface components (for slab correction in slit pore)
	  Double Precision, Dimension(:), Allocatable :: ewld_surfc, ewld_surfc_new, ewld_slab


	  ! Spherical cut-off treatment for adsorbate-surface interaction
	  !================================================!

	  ! Flag to turn on/off this feature
	  Logical, save :: lcoulsc

	  ! Spherical cut-off radius 
	  Double Precision, save :: rscelect, rscelectsq




	  ! Thermodyanmic properties
	  !================================================!

	  ! Temperature [K]
	  Double Precision, save :: temp 

	  ! Pressure [bar]
	  Double Precision, save :: press 

	  ! Pre-defined chemical potential for each type of molecule [K]
	  Double Precision, Dimension(:), Allocatable :: mu

	  ! Update step for chemical potential [K] (add 8-13-2017 for update_chempot subroutine)
	  Double Precision, Dimension(:), Allocatable :: dmu

	  ! Total energy of each box [K]
	  Double Precision, Dimension(:), Allocatable :: energy 

	  ! Total dispersion energy of each box [K]
	  Double Precision, Dimension(:), Allocatable :: energy_disp

	  ! Instant virial [K] of each box
	  Double Precision, Dimension(:), Allocatable :: vir

	  ! Instant hypervirial [K] of each box
	  Double Precision, Dimension(:), Allocatable :: vir_hyp



	  ! Monte Carlo Parameters
	  !================================================!

	  ! Total number of simulation steps
	  Integer, save :: n_steps_tot

	  ! Number of equilibrium blocks
	  Integer, save :: n_blocks_equil

	  ! Number of production blocks
	  Integer, save :: n_blocks_prod

	  ! Total Number of blocks
	  Integer, save :: n_blocks_tot

	  ! Number of relax blocks (added on Aug 12, 2017)
	  Integer, save :: n_blocks_relax

	  ! Number of steps in each block
	  Integer, save :: block_size

	  ! PDB output frequency
	  Integer, save :: pdb_out_freq 

	  ! XYZ output frequency
	  Integer, save :: xyz_out_freq

	  ! Coordinates file output frequency in production stage
	  Integer, save :: coords_out_freq

	  
	  ! Statistical variables 
	  !================================================!

	  ! probability of translational move for each type of molecule
	  Double Precision, Dimension(:), Allocatable :: trans_prob

	  ! Probability of rotational move for each type of molecule
	  Double Precision, Dimension(:), Allocatable :: rotat_prob

	  ! Prabability of transfer move for each type of molecule
	  Double Precision, Dimension(:), Allocatable :: transfer_prob

	  ! Probability for selecting each move type
	  Double Precision, Dimension(10), save :: move_type_prob

	  ! Probability for selecting each move type (use for update_prob subroutine)
	  Double Precision, Dimension(2,10), save :: move_type_prob_update

	  ! Maximum translation and rotation for each molecule in each box
	  Double Precision, Dimension(:,:), Allocatable :: max_trans, max_rotat

	  ! Maximum volume change
	  Double Precision, Dimension(:), Allocatable :: max_vol

	  ! Move acceptance statistics
	  Integer, Dimension(:,:), Allocatable :: in_stat, rem_stat
	  Integer, Dimension(:,:), Allocatable :: vol_stat, swap_ident_stat
	  Integer, Dimension(:,:,:), Allocatable :: trans_stat, rotat_stat, transfer_stat


	  ! Block average statistics
	  !================================================!

	  Integer, Dimension(:,:,:), Allocatable :: blk_nmol
	  Double Precision, Dimension(:,:,:), Allocatable :: blk_rho
	  Double Precision, Dimension(:,:), Allocatable :: blk_vol
	  Double Precision, Dimension(:,:), Allocatable :: blk_eng, blk_eng_disp

	  Double Precision, Dimension(:,:,:), Allocatable :: avg_nmol, avg_rho
	  Double Precision, Dimension(:,:), Allocatable :: avg_vol
	  Double Precision, Dimension(:,:), Allocatable :: avg_eng, avg_eng_disp




	  ! Optional sampling
	  !================================================!

	  ! Optional sampling flag
	  Logical, save :: lsampling


	  ! ------------------ z_density statistics 
	  ! (Density profile of each molecule type along z-axis)
	  ! Only perform during production stage
	  ! ONLY for box 1 right now and Only could be applied to volume-fixed ensemble

	  ! Flag to turn on/off z_density
	  Logical, save :: lzdensity

	  ! Number of bins to produce z-density
	  Integer, save :: zden_bins

	  ! Sampling frequency
	  Integer, save :: zden_freq

	  ! Counter 
	  Double Precision, save :: zden_stat

	  ! z-density volume unit and length in z direction for each bin
	  Double Precision, save :: dvol, dz

	  ! Statistics 
	  Double Precision, Dimension(:,:,:), Allocatable :: blk_zden, avg_zden


	  ! ----------------- Surface excess statistics (Excess adsorption number)
	  ! Only perform during production stage
	  ! ONLY for box 1 right now and Only could be applied to volume-fixed ensemble

	  ! Flag to turn on/off surface_excess
	  Logical, save :: lsurfex

	  ! Accessible volume [A^3]
	  Double Precision, save :: surfex_vol

	  ! Surface area [A^2]
	  Double Precision, save :: surfex_area

	  ! Pre-set bulk number density [1/A^3]
	  Double Precision, Dimension(:), Allocatable :: surfex_bulk

	  ! Number of molecules in accessible volume
	  Double Precision, Dimension(:), Allocatable :: surfex_n_mol

	  ! Statistics
	  Double Precision, Dimension(:,:), Allocatable :: blk_surfex, avg_surfex


	  ! ------------- Thermo route to calculate pressure tensor for slit pore geometry
	  ! P = - dF/dV (Pressure tensor equation for GCMC is the same as that for NVT)
	  ! Use zden_bins for the number of bins
	  ! Flag to turn on/off pressure calculation
	  Logical, save :: lthermopress_slit

	  ! Lx/Ly/Lz ratio change factor
	  ! fac = dL/L
	  Double Precision, save :: thermopress_slit_ratio

	  ! Volume change frequency
	  Integer, save :: thermopress_slit_freq

	  ! Number of bins to produce pressure tensor
	  Integer, save :: thermopress_slit_bins

	  ! dvol and dz for each bin
	  Double Precision, save :: thermopress_slit_dvol, thermopress_slit_dz 

	  ! Volume change counter
	  Double Precision, Dimension(:), Allocatable :: thermopress_slit_stat

	  ! Statistics
	  ! 1 - <exp(-beta*dU)>. 2 - <P>
	  Double Precision, Dimension(:,:,:), Allocatable :: thermopress_slit_sample
	  Double Precision, Dimension(:,:,:), Allocatable :: thermopress_slit_pt, thermopress_slit_pn


	  ! ---------------- Virial route to calculate pressure tensor for slit pore geometry
	  ! Definition of integral contour type. 
	  ! 1-Irving-Kirkwood, 2-Harasima, 3-both, 4-H-VR1 (added on May 29,2019), 5-IK-VR1, 6-ALL (May 30, 2019)
	  Integer, save :: virialpress_ctype

	  ! Definition of the contribution to be accounted
	  ! 1 - fluid-fluid, 2- fluid-wall, 3-both
	  Integer, save :: calc_type

	  ! Flag to turn on/off pressure calculation
	  Logical, save :: lvirialpress_slit

	  ! Calculation frequency
	  Integer, save :: virialpress_slit_freq

	  ! Number of bins to produce pressure tensor
	  Integer, save :: virialpress_slit_bins

	  ! dvol and dz for each bin
	  Double Precision, save :: virialpress_slit_dz 

	  ! Volume change counter
	  Double Precision, Dimension(:), Allocatable :: virialpress_slit_stat

	  ! Statistics
	  Double Precision, Dimension(:,:,:,:), Allocatable :: virialpress_slit_pn
	  Double Precision, Dimension(:,:,:,:,:), Allocatable :: virialpress_slit_pt



	  ! -------------- Virial route to calculate pressure tensor for planar surface geometry (Added on March 1, 2019)
	  ! Definition of integral contour type. 1-Irving-Kirkwood, 2-Harasima, 3-both
	  ! Flag to turn on/off pressure calculation
	  Logical, save :: lvirialpress

	  ! Calculation frequency
	  Integer, save :: virialpress_freq

	  ! Volume change counter
	  Double Precision, Dimension(:), Allocatable :: virialpress_stat

	  ! Statistics
	  ! 1-fluid-wall; 2-fluid-fluid; 3-total
	  Double Precision, Dimension(:,:,:), Allocatable :: virialpress_pn
	  Double Precision, Dimension(:,:,:,:), Allocatable :: virialpress_pt



	  ! --------------- Virial route to calculate pressure tensor for cylindrical geometry (added on Feb 27 2019)
	  ! Flag to turn on/off pressure calculation
	  Logical, save :: lvirialpress_cylin

	  ! Calculation frequency
	  Integer, save :: virialpress_cylin_freq

	  ! delta value in equation for Harasima contour use
	  ! Defined in 'initialize.f90' to be half of rden_dr
	  Double Precision :: delrr

	  ! counter
	  Double Precision, Dimension(:), Allocatable :: virialpress_cylin_stat

	  ! Statistics
	  ! 1-fluid-wall; 2-fluid-fluid; 3-total
	  Double Precision, Dimension(:,:,:,:,:), Allocatable :: virialpress_cylin_pnr
	  Double Precision, Dimension(:,:,:,:,:), Allocatable :: virialpress_cylin_ptt
	  Double Precision, Dimension(:,:,:,:,:), Allocatable :: virialpress_cylin_ptz


	  ! --------------- r-density statistics (Density profile of each molecule type in the cylindrical system)
	  ! Only perform during production stage
	  ! ONLY for box 1 right now and Only could be applied to volume-fixed ensemble
	  ! Added on Feb 27 2019 

	  ! Flag to turn on/off r_density
	  Logical, save :: lrdensity

	  ! Cutoff radius for density and cylindrical pressure tensor calculation
	  Double Precision, save :: rden_cut, rden_lim

	  ! Number of bins to produce r-density
	  Integer, save :: rden_bins

	  ! Sampling frequency
	  Integer, save :: rden_freq

	  ! Counter 
	  Double Precision, save :: rden_stat

	  ! r-density length in r direction for each bin
	  Double Precision, save :: rden_dr, rden_drsq

	  ! Statistics 
	  Double Precision, Dimension(:,:,:), Allocatable :: blk_rden, avg_rden



	  ! ---------------- 2D lattice constant statistics
	  ! Only the lattice constant within the same adsorbed layer will be calculated
	  ! Flag to turn on/off lattice constant calculation
	  Logical, save :: llattconst

	  ! Number of layers to be considered
	  Integer, save :: lattconst_n

	  ! Boundary for the confined layer
	  ! 1 - lower; 2 - upper
	  Double Precision, Dimension(2,100) :: lattconst_cut

	  ! Frequency to do averaging 
	  Integer, save :: lattconst_freq

	  ! Statistics
	  Double Precision, Dimension(:), Allocatable :: lattconst_stat
	  Double Precision, Dimension(:,:), Allocatable :: blk_lattconst
	  Double Precision, Dimension(:,:), Allocatable :: avg_lattconst
	  Double Precision, Dimension(:), Allocatable :: lattconst


	  ! ---------------- Isosteric heat of adsorption
	  ! Now, ONLY works for single adsorbate component
	  ! Needs to be revised when applied to multi-components
	  Logical, save :: lqst

	  ! Isosteric heat 
	  Double Precision :: qst

	  ! Statistics
	  Double Precision, Dimension(:,:,:), Allocatable :: blk_qst, avg_qst


	  ! ------------- Output coordinates of the whole system to file
	  ! Flag to turn on/off dump_xyz option
	  Logical, save :: ldumpxyz

	  ! Specified molecule type to dump
	  Integer, save :: dump_moltype

	  ! Specify dumping modes
	  Character(Len=6), save :: dumpxyz_mode

	  ! Frequency to dump
	  Integer, save :: dumpxyz_freq

	  ! -------------- Output density of the whole system to file
	  ! Flag to turn on/off dump_density option
	  Logical, save :: ldumpdensity

	  ! Frequency to dump
	  Integer, save :: dumpdensity_freq


	  ! -------------- Output total energy of the system to file
	  ! Added on Oct 13, 2019
	  ! Flag to turn on/off dump_energy option
	  Logical, save :: ldumpenergy

	  ! Frequency to dump
	  Integer, save :: dumpenergy_freq

	  ! ------------ Output instant virial, hypervirial and hydrostatic pressure
	  ! Added on Jan, 2020
	  ! Flag to turn on/off dump_virial option
	  Logical, save :: ldumpvir

	  ! Frequency to dump
	  Integer, save :: dumpvir_freq


	  ! ------------ Print instantaneous configuration for restart 
	  ! Added on Oct 13, 2019
	  Logical, save :: lwriterestart

	  ! Frequency to rewrite
	  Integer, save :: rst_freq


	  ! -------------- Turn off Periodic boundary condition
	  ! Flag to turn off PBC option and turn on a hard wall boundary condition
	  ! Added on July 14, 2019
	  Logical, save :: lno_pbc


	  ! -------------- Use the cell list and linked-list algorithm 
	  ! to accelerate dispersive (short-range) energy calculations 
	  ! Added on July 15, 2019
	  ! Flag to turn on/off cell list
	  Logical, save :: lclist

	  ! head of chain in the linked-list algorithm (reference page 552, Frenkel & Smit book)
	  Integer, Dimension(:), Allocatable :: clist_hoc, clist_hoc_sub

	  ! linked list 
	  Integer, Dimension(:), Allocatable :: clist_llist, clist_llist_sub

	  ! Neighbor cells 
	  Integer, Dimension(:,:), Allocatable :: clist_neigh

	  ! Cell locator
	  Integer, Dimension(:,:,:), Allocatable :: clist_loca

	  ! Number of cells in each dimension
	  Integer, save :: clist_nx, clist_ny, clist_nz, clist_ncel

	  ! Size of cell in each dimension
	  Double Precision, save :: clist_dx, clist_dy, clist_dz


	  ! --------------- Rolling energy check ----------------
	  ! Added on 6-9-2020 
	  ! To recalculate and check the total energy to mitigate the error
	  ! check frequency, default value 5000
	  Integer, save :: check_freq = 5000




	  ! Widom Insertion Method 
	  !================================================!
	  ! Note: Allocated in global_allocate.f90, 
	  !		  initialized in initialize.f90 and stat.f90

	  ! Flag turn on/off widom insertion method
	  Logical, save :: lwdm

	  ! Insertion Frequency for each molecule type
	  Integer, Dimension(:,:), Allocatable :: widom_freq

	  ! Ideal chemical potential of each molecule type [K]
	  Double Precision, Dimension(:), Allocatable :: muid

	  ! Total chemical potential [K]
	  Double Precision, Dimension(:,:,:), Allocatable :: widom_sample

	  ! Insertion counter
	  Integer, Dimension(:,:), Allocatable :: widom_stat

	  
	  ! Random number generator 
	  !================================================!

	  Integer :: idum
	  Real, External :: random 

	  ! Utility function declaration 
	  !================================================!

	  ! Check label function for input files
	  Logical, External :: label_check


	  ! Predifined parameters 
	  !================================================!

	  ! Ensemble type
	  Integer, Parameter :: ENS_NVT = 1
      Integer, Parameter :: ENS_NPT = 2
      Integer, Parameter :: ENS_uVT = 3



      ! File number, number setting is arbitary*/
      ! File unit 6 is screen
      ! input.in
      Integer, Parameter :: FILE_INPUT = 1
      ! coords.in
      Integer, Parameter :: FILE_COORDS = 2
      ! mol.in
      Integer, Parameter :: FILE_MOL = 3
      ! config_MOL.xyz file
      Integer, Parameter :: FILE_XYZ_MOL = 4
      ! config_site.xyz file
      Integer, Parameter :: FILE_XYZ_SITES = 5
      ! old_config file
      Integer, Parameter :: FILE_OLDCONFIG = 7
      ! z-density.txt output file
      Integer, Parameter :: FILE_ZDENSITY = 8
      ! surface_excess.txt output file
      Integer, Parameter :: FILE_SURFEX = 9
      ! press_thermo_slit.txt output file
      Integer, Parameter :: FILE_THERMOPRESS_SLIT = 10
      ! Pre-calculated coarse-grained solid-fluid potential input file
      Integer, Parameter :: FILE_CG_WALL = 13
      ! Pre-calculated effective sigma for fluid-fluid potential input file
      Integer, Parameter :: FILE_EFF_SIG = 14
      ! Pre-calculated accessible volume/length from full structure
      Integer, Parameter :: FILE_CG_WALL_STRC = 15
      ! Pre-calculated coarse-grained solid-fluid potential input file for modified potential
      Integer, Parameter :: FILE_CG_WALL1 = 16
      ! lattconst.txt
      Integer, Parameter :: FILE_LATTCONST = 17
      ! dump_qst.txt
      Integer, Parameter :: FILE_QST = 18
      ! dump_xyz.txt
      Integer, Parameter :: FILE_XYZ_ALL = 19
      ! r-density.txt output file
      Integer, Parameter :: FILE_RDENSITY = 20
      ! press_virial_cylinH.txt (default for Harasima definition)
      Integer, Parameter :: FILE_VIRIALPRESS_CYLINH = 21
      ! density.txt 
      Integer, Parameter :: FILE_DENSITY = 22
      ! press_slit_1.txt output file (IK definition)
      Integer, Parameter :: FILE_VIRIALPRESS_SLIT_1 = 23
      ! press_slit_2.txt output file (Harasima definition)
      Integer, Parameter :: FILE_VIRIALPRESS_SLIT_2 = 24
      ! press_slit_3.txt output file (H-VR1 definition)
      Integer, Parameter :: FILE_VIRIALPRESS_SLIT_3 = 25
      ! press_slit_4.txt output file (IK-VR1 definition)
      Integer, Parameter :: FILE_VIRIALPRESS_SLIT_4 = 26
      ! press_1.txt output file (IK definition for planar interface, bulk for now)
      Integer, Parameter :: FILE_VIRIALPRESS_1 = 27
      ! press_2.txt output file (Harasima definition for planar interface, bulk for now)
      Integer, Parameter :: FILE_VIRIALPRESS_2 = 28
      ! energy.txt
      Integer, Parameter :: FILE_ENERGY = 29
      ! rst1.xyz
      Integer, Parameter :: FILE_RST = 30
      ! rst2.xyz
      Integer, Parameter :: FILE_RST_CP = 31
      ! virial.txt
      Integer, Parameter :: FILE_VIR = 32


      ! MC Move type
      Integer, Parameter :: TRANSLATION = 1
      Integer, Parameter :: ROTATION = 2
      Integer, Parameter :: TRANSFER = 3
      Integer, Parameter :: VOLCHANGE = 4

      ! Ewald Style
      Integer, Parameter :: ewald_fix_kmax = 1


      ! CLASSICAL POTENTIAL ALIASES
      ! Maps classical potential names to integers to speed
      ! up inner loops but retain readability.
      ! Changing these settings is NOT advised  

	  Integer, Parameter :: LENNARD_JONES = 1
	  Integer, Parameter :: HARD_SPHERE = 2

	  ! External Field Types
	  Integer, Parameter :: HARD_WALL = 1
	  Integer, Parameter :: STEELE = 2
	  Integer, Parameter :: STEELE_SLIT_PORE = 3
	  Integer, Parameter :: CG_WALL = 4
	  Integer, Parameter :: CG_WALL_FFPW = 5
	  Integer, Parameter :: CG_WALL_COS = 6
	  Integer, Parameter :: CG_WALL_STRC = 7
	  Integer, Parameter :: STEELE_SLIT_FINITEX = 8
	  Integer, Parameter :: HARD_SLIT_FINITEX = 9


	  ! MIXING RULES 
      ! setting numbers are arbitary 

	  Integer, Parameter :: LORENTZ_BERTHELOT = 11
	  Integer, Parameter :: EXPLICIT = 22

	  ! Indicator for qex and qsqex (Explicit mixing rule)
	  ! Parameters for adsorbate-surface interaction
	  Integer, Parameter :: AS = 1
	  ! Parameters for adsorbate-adsorbate interaction
	  Integer, Parameter :: AA = 2



	  ! Constants [CODATA 2010 recommendation]
	  !================================================!
      Double Precision, Parameter :: Pi = 3.141592653589d0
 	  Double Precision, parameter :: two_Pi = 2.0d0*Pi
      Double Precision, Parameter :: Kb = 1.3806488d-23   ! [J/K]; true units should be [J/(K*molecule)]
	  Double Precision, Parameter :: Na = 6.02214129d23  ! [1/mol]; true units should be [molecule/mol]
	  Double Precision, Parameter :: h = 6.62606957d-34   ! [J*s] = [kg*m^2/s]
	  Double Precision, Parameter :: R = Kb*Na            ! [J/(mol*K)] */
	  Double Precision, Parameter :: EETOK = 167100.9558738398d0 ! EETOK = 1/(4*pi*epsilon*kb) * e^2 * 10^10
	  Double Precision, Parameter :: PVTOK = 0.7242971565d-02   ! PVTOK = 1.0d5*1.0d-30/kb , [bar]*[A^3] to [K]
	  Double Precision, Parameter :: PCOEFF = Kb*1.0d25        ! [K]/[A]^3 to [bar]
	  
	  !================================================!
	  



	  End Module global





      



