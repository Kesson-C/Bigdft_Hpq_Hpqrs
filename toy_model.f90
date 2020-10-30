!> @file
!! @author
!!    Copyright (C) 2007-2013 BigDFT group. This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file or http://www.gnu.org/copyleft/gpl.txt.
!!    For the list of contributors, see ~/AUTHORS
!
!===============================================================================
! JP NOTES:
!   * input files: wavefunction files in data/ folder
!                  input.yaml
!                  input-hpqrs (112)
!    output files: toy_model.log (124)
!                  toy_model.out (6)
!
!   * The number of virtual orbitals must be less than occupied orbitals due to 
!     the dimension of virtual orbital is identical as occupied orbitals

!  gp = wp = dp = 8
!===============================================================================

!===============================!
!> Toy program to use BigDFT API!
!===============================!
program toy_model
  use Poisson_Solver, except_dp => dp, except_gp => gp
  use BigDFT_API
  use bigdft_run
  use dynamic_memory
  use compression
  use yaml_output
  use module_input_dicts
  use module_input_keys
  use module_razero
  use module_atoms
  use module_dpbox,       only: denspot_distribution,dpbox_free,dpbox_set
  use rhopotential,       only: full_local_potential
  use locregs_init,       only: lr_set
  use locreg_operations 
  use communications_base,  only: deallocate_comms
  use communications_init,  only: orbitals_communicators
  use communications,       only: transpose_v, untranspose_v
  use module_base
  use module_types
  use psp_projectors
  use pseudopotentials
  use orbitalbasis
  use ao_inguess, only: lmax_ao
  use dictionaries

  implicit none
  type(input_variables)        :: inputs
  type(atoms_data)             :: atoms
  type(local_zone_descriptors) :: Lzd
  type(orbitals_data)          :: orbs, orbsv
  type(comms_cubic)            :: comms
  type(workarr_sumrho)         :: wisf
  type(xc_info)                :: xc
  type(rho_descriptors)        :: rhodsc
  type(denspot_distribution)   :: dpcom
  type(GPU_pointers)           :: GPU
  type(workarr_locham)         :: wrk_lh
  type(coulomb_operator)       :: pkernel
  type(paw_objects)            :: paw_r, paw_i
  type(DFT_PSP_projectors)     :: nlpsp_r, nlpsp_i
  type(ket)                    :: psi_it_r, psi_it_i
  type(orbital_basis)          :: psi_ob_r, psi_ob_i
  type(DFT_PSP_projector_iter) :: psp_it_r, psp_it_i
  type(atomic_proj_matrix)     :: prj
  type(dictionary), pointer :: user_inputs, options, dict

  logical  :: dosome, alive, init_projectors_completely
! logical  :: paw = .false.
  real(dp) :: nrm, epot_sum_r,epot_sum_i, hxh,hyh,hzh
  real(dp) :: psoffset, ekin, ekin_sum, eproj_sum
  real(gp), dimension(3)            :: shift
  real(wp), dimension(:),   pointer :: w
  real(wp), dimension(:,:),   pointer :: psi, psiv, psir, psir_i,psir_j  ! psi(1,ispinor)
  real(dp), dimension(:),   pointer :: rhor, pot_ion, potential,rho_ion
  real(wp), dimension(:)  , pointer :: hpsi_ptr_r, hpsi_r, hpsi_ptr_i, hpsi_i
  real(gp), dimension(:,:), pointer :: rxyz_old
  real(wp), dimension(:,:), pointer :: ovrlp
  real(dp), dimension(:,:), pointer :: rho_p => null() !needs to be nullified
  real(wp), dimension(:,:), allocatable :: pot_tmp
  real(wp), dimension(:,:,:), allocatable :: eprjo1, eprjo2
  real(wp), dimension(:,:,:), allocatable :: eprjv1, eprjv2
  real(wp), dimension(:,:,:), allocatable :: psirr
  real(wp), dimension(:,:)  , allocatable :: tpsi, tpsi_o, tpsi_v
  real(wp), dimension(:,:)  , allocatable :: tho,tho1,tho2, tho_new1,tho_new2
  real(dp), dimension(:,:,:), allocatable  :: E_local, E_nonlocal, E_kin
  real(dp), dimension(:,:,:), allocatable  :: output1, output2, output3
  real(dp)                               :: eproj

  integer(kind=8), dimension(:), allocatable :: indwvl_i,indwvl_f
  integer,         dimension(:,:), allocatable :: thoi1,thoi2
  integer :: ierr, iproc, nproc, nwarnings, ii,jj,kk, nn,nn_new1,nn_new2,nn_new
  integer :: iorb, nvirtu, nvirtd, ispsi, ilr, ilr_orb, npot, mspin, cspin
  integer :: ityp, nc, m, mm, ishift, nwarning, dry_run, orbtot, nv, istate
  integer :: orbdimocc, orbdimvir, norb_start,norb_end,  orbocc,orbvir

  character*32 con, tmp
  integer(kind=8) :: ip ,iq ,ir ,is , ihpqrs, nhpqrs, nkpqrs, nrpqrs, istat, i,j,k, re_index1, re_index2
  integer(kind=8) :: ipt,iqt,irt,ist
  integer(kind=8), dimension(3) :: irr
  real(kind=8), dimension(:), allocatable :: nii, njj, nkk
  real(kind=8)                        :: rx, ry, rz, rr, EH, r_r_r, r_r_i, ee1, ee2_r, ee2_i, ee3_r, ee3_i, cri
  real(kind=8)                        :: hpqrs_r, hpqrs_i, hpqrs_t1, hpqrs_t2, hpqrs_t3, hpqrs_t4, hpqrs_t5, kpqrs, rpqrs, ecos, esin, eikr
  integer(kind=4)                     :: OMP_get_max_threads, OMP_get_thread_num, OMP_get_num_threads
  
  open(124,file="toy_model.log") ; open(  6,file="toy_model.out")

  !-----------------------------------------------
  ! initializes the mpi_environment and bigdft
  !-----------------------------------------------
  call f_lib_initialize() ; nullify(options)
  call bigdft_init(options) ; call dict_free(options)

  iproc = bigdft_mpi%iproc
  nproc = bigdft_mpi%nproc
  call dict_init(user_inputs)
  call user_dict_from_files(user_inputs, 'input', 'posinp', bigdft_mpi)
  call inputs_from_dict(inputs, atoms, user_inputs)

  if (iproc == 0) call print_general_parameters(inputs,atoms,'input')
  call dict_free(user_inputs)
  GPU%OCLconv = .false.

  call system_properties(iproc,nproc,inputs,atoms,orbs)
  Lzd = default_Lzd() ; Lzd%hgrids=(/ inputs%hx, inputs%hy, inputs%hz /)
  call lr_set(Lzd%Glr, iproc, GPU%OCLconv, .true., inputs%crmult, inputs%frmult, &
              Lzd%hgrids,atoms%astruct%rxyz,atoms,.true.,.false.)
  call orbitals_communicators(iproc,nproc,Lzd%Glr,orbs,comms)
  call xc_init(xc, inputs%ixc, XC_LIBXC , inputs%nspin)
  call dpbox_set(dpcom,Lzd%Glr%mesh,xc,iproc,nproc,MPI_COMM_WORLD, inputs%SIC%approach, inputs%nspin)

!-------------------------------------------------------------------------------------------------------------------------------------
!=====================================================================================================================================
!-------------------------------------------------------------------------------------------------------------------------------------

  if(orbs%nspin .eq. 1) cspin = 2 ; if(orbs%nspin .eq. 1) mspin = 1
  if(orbs%nspin .eq. 2) cspin = 1 ; if(orbs%nspin .eq. 2) mspin = 2
  orbocc = orbs%norb
  orbvir = mspin*inputs%norbv
  orbtot = orbs%norb + mspin*inputs%norbv
  allocate( E_kin(      orbtot,orbtot,orbs%nspinor) ) ; E_kin=0.d0
  allocate( E_local(    orbtot,orbtot,orbs%nspinor) ) ; E_local=0.d0
  allocate( E_nonlocal( orbtot,orbtot,orbs%nspinor) ) ; E_nonlocal=0.d0

  allocate(orbs%eval(orbs%norb*orbs%nkpts)) ; call f_zero(orbs%eval)
  allocate(psi(  max(orbs%npsidim_orbs, orbs%npsidim_comp)+1 ,orbs%nspinor)) ; psi=0._gp
  allocate(psiv( 2*inputs%norbv * (Lzd%Glr%wfd%nvctr_c + 7*Lzd%Glr%wfd%nvctr_f)+1 ,orbs%nspinor)) ; psiv=0._gp
  allocate(rxyz_old(3, atoms%astruct%nat)) ; rxyz_old=0._gp

  !--------------------------------------------------------------------!
  ! Read occupied state wavefunctions from disk and store them in psi. !
  ! orbs%norb - number of occupied orbitals                            !
  !--------------------------------------------------------------------!
  call check_linear_and_create_Lzd(iproc,nproc,inputs%linear,Lzd,atoms,orbs,inputs%nspin,atoms%astruct%rxyz)
  call readmywaves(iproc, "data/wavefunction", WF_FORMAT_PLAIN, orbs, Lzd%Glr%d%n1, Lzd%Glr%d%n2, Lzd%Glr%d%n3, &
                   inputs%hx, inputs%hy, inputs%hz, atoms, rxyz_old, atoms%astruct%rxyz, Lzd%Glr%wfd, psi)
  if(nproc>1) call fmpi_allreduce(orbs%eval(1), orbs%norb*orbs%nkpts, op=FMPI_SUM)
! do i=1,orbs%norb ; write(1000+i,'("# orb: ",i4)') i ; write(1000+i,'(f20.12)') psi ; end do

  !--------------------------------------------------------------------!
  ! Read virtual  state wavefunctions from disk and store them in psi. !
  ! inputs%norbv - number of virtual orbitals                          !
  ! orbsv%norb - number of total virtual orbitals (2*inputs%norbv)     !
  !--------------------------------------------------------------------!
  nullify(orbsv%eval)
  orbsv%eval = f_malloc_ptr(orbsv%norb*orbsv%nkpts,id='orbsv%eval')
  nvirtu = abs(inputs%norbv) ; nvirtd = nvirtu
  call orbitals_descriptors(iproc, nproc, nvirtu+nvirtd, nvirtu, nvirtd, &
                            orbs%nspin, orbs%nspinor, orbs%nkpts, orbs%kpts, orbs%kwgts, orbsv, LINEAR_PARTITION_NONE)
  call check_linear_and_create_Lzd(iproc,nproc,inputs%linear,Lzd,atoms,orbsv,inputs%nspin,atoms%astruct%rxyz)
  call readmywaves(iproc, "data/virtuals",     WF_FORMAT_PLAIN, orbsv, Lzd%Glr%d%n1, Lzd%Glr%d%n2, Lzd%Glr%d%n3, &
                   inputs%hx, inputs%hy, inputs%hz, atoms, rxyz_old, atoms%astruct%rxyz, Lzd%Glr%wfd, psiv)
  if(nproc>1) call fmpi_allreduce(orbsv%eval(1), orbsv%norb*orbsv%nkpts, op=FMPI_SUM)
! do i=1,mspin*inputs%norbv ; write(2000+i,'("# orb: ",i4)') i ;  write(2000+i,'(f20.12)') psiv ; end do

  orbdimocc = (Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)*orbs%norb
  orbdimvir = (Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)*orbsv%norb

!-------------------------------------------------------------------------------------------------------------------------------------
!=====================================================================================================================================
!-------------------------------------------------------------------------------------------------------------------------------------

    call system('echo "kinetic calculation"')
    !--------------------------------------
    ! kinetic energy calculation
    !--------------------------------------
    write(124,*) "================"
    write(124,*) " kinetic energy "
    write(124,*) "================"
    norb_start = orbs%norb+1
    norb_end   = orbs%norb+mspin*inputs%norbv
    allocate( tpsi_o(orbdimocc+1,orbs%nspinor), indwvl_i(orbtot) ) ; tpsi_o = 0._dp
    allocate( tpsi_v(orbdimvir+1,orbs%nspinor), indwvl_f(orbtot) ) ; tpsi_v = 0._dp
  
    !------------------------------------------------------
    ! loop on the localisation regions (occupied orbitals)
    !------------------------------------------------------
    allocate(tpsi(orbdimocc,orbs%nspinor)) ; tpsi = 0._dp
    write(124,20) orbs%norb, orbs%nspinor, mspin, orbdimocc, Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f
    call initialize_work_arrays_locham(Lzd%nlr,Lzd%Llr,orbs%nspinor,.true.,wrk_lh)
    ekin = 0.d0 ; ekin_sum = 0.0_gp 
    loop_lr_kin: do ilr=1,Lzd%nlr
      dosome=.false.   ! check if this localisation region is used by one of the orbitals
      do iorb=1,orbs%norb ; dosome = (orbs%inwhichlocreg(iorb+orbs%isorb) == ilr) ; if(dosome) exit ; end do
      if (.not. dosome) cycle loop_lr_kin
      call initialize_work_arrays_locham(Lzd%Llr(ilr),orbs%nspinor,.false.,wrk_lh)
      ispsi = 1
      loop_orbs: do iorb=1,orbs%norb
        indwvl_i(iorb) = ispsi
        ilr_orb = orbs%inwhichlocreg(iorb+orbs%isorb)
        ilr_orb = 1
        if (ilr_orb /= ilr) then
          ispsi = ispsi + (Lzd%Llr(ilr_orb)%wfd%nvctr_c+7*Lzd%Llr(ilr_orb)%wfd%nvctr_f)*orbs%nspinor
          cycle loop_orbs
        end if
        call psi_to_tpsi(Lzd%hgrids, orbs%kpts(1,orbs%iokpt(iorb)),orbs%nspinor,Lzd%Llr(ilr), psi(ispsi,:), wrk_lh, tpsi(ispsi,:), ekin)
        ekin_sum = ekin_sum + orbs%kwgts(orbs%iokpt(iorb))*orbs%occup(iorb+orbs%isorb)*ekin
        ispsi  = ispsi + (Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*orbs%nspinor
        indwvl_f(iorb) = ispsi
        write(124,22) iorb, ekin, indwvl_i(iorb),indwvl_f(iorb)
      end do loop_orbs
    end do loop_lr_kin
    write(124,*) " occupied orbitals ... DONE"
    write(124,'("total kinetic energy = ",f15.8)') ekin_sum ; write(124,*)
    tpsi_o=tpsi
    call deallocate_work_arrays_locham(wrk_lh)
    deallocate(tpsi)

    !------------------------------------------------------
    ! loop on the localisation regions ( virtual orbitals)
    !------------------------------------------------------
    allocate(tpsi(orbdimvir,orbs%nspinor)) ; tpsi = 0._dp
    write(124,21) inputs%norbv*mspin, orbs%nspinor, mspin, orbdimvir, Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f
    call initialize_work_arrays_locham(Lzd%nlr,Lzd%Llr,orbs%nspinor,.true.,wrk_lh) 
    ekin = 0.d0 ; ekin_sum = 0.0_gp ; tpsi=0._dp
    loop_lr_kinv: do ilr=1,Lzd%nlr
      dosome=.false.   ! check if this localisation region is used by one of the orbitals
      do iorb=1,inputs%norbv ; dosome = (orbs%inwhichlocreg(iorb+orbs%isorb) == ilr) ; if(dosome) exit ; end do
      if (.not. dosome) cycle loop_lr_kinv
      call initialize_work_arrays_locham(Lzd%Llr(ilr),orbs%nspinor,.false.,wrk_lh)
      ispsi = 1
      loop_orbsv: do iorb=1,mspin*inputs%norbv
        indwvl_i(iorb+orbs%norb) = ispsi
        ilr_orb = orbs%inwhichlocreg(iorb+orbs%isorb)
        ilr_orb = 1
        if (ilr_orb /= ilr) then
          ispsi = ispsi + (Lzd%Llr(ilr_orb)%wfd%nvctr_c+7*Lzd%Llr(ilr_orb)%wfd%nvctr_f)*orbs%nspinor
          cycle loop_orbsv
        end if
        call psi_to_tpsi(Lzd%hgrids, orbs%kpts(1,orbs%iokpt(iorb)),orbs%nspinor,Lzd%Llr(ilr), psiv(ispsi,:), wrk_lh, tpsi(ispsi,:), ekin)
        ekin_sum = ekin_sum + ekin
        ispsi  = ispsi + (Lzd%Llr(ilr)%wfd%nvctr_c+7*Lzd%Llr(ilr)%wfd%nvctr_f)*orbs%nspinor
        indwvl_f(iorb+orbs%norb) = ispsi
        write(124,22) iorb+orbs%norb, ekin, indwvl_i(iorb+orbs%norb),indwvl_f(iorb+orbs%norb)
      end do loop_orbsv
    end do loop_lr_kinv
    write(124,*) " virtual  orbitals ... DONE"
    write(124,'("total kinetic energy = ",f15.8)') ekin_sum ; write(124,*)
    tpsi_v=tpsi
    call deallocate_work_arrays_locham(wrk_lh)
    deallocate(tpsi)

  20 format("   # occupied orbitals:",i3,", nspinor:",i3,", mspin: ",i4,", psi_dim:",3i9)
  21 format("   # virtual  orbitals:",i3,", nspinor:",i3,", mspin: ",i4,", psi_dim:",3i9)
  22 format(" orbital ",i3,"  ekin = ",f15.8,2x,"index of wavelet:",2i10)
        !hij
    do i=1,   orbs%norb       ; do j=1,  orbs%norb
      E_kin(i,j,1) = sum( psi( indwvl_i(i):indwvl_f(i)-1,1) * tpsi_o(indwvl_i(j):indwvl_f(j)-1,1) ) + sum( psi( indwvl_i(i):indwvl_f(i)-1,2) * tpsi_o(indwvl_i(j):indwvl_f(j)-1,2) )
      E_kin(i,j,2) = -sum( psi( indwvl_i(i):indwvl_f(i)-1,2) * tpsi_o(indwvl_i(j):indwvl_f(j)-1,1) ) + sum( psi( indwvl_i(i):indwvl_f(i)-1,1) * tpsi_o(indwvl_i(j):indwvl_f(j)-1,2) )
    end do ; end do
    do i=1,   orbs%norb       ; do j=norb_start,norb_end
      E_kin(i,j,1) = sum( psi( indwvl_i(i):indwvl_f(i)-1,1) * tpsi_v(indwvl_i(j):indwvl_f(j)-1,1) ) + sum( psi( indwvl_i(i):indwvl_f(i)-1,2) * tpsi_v(indwvl_i(j):indwvl_f(j)-1,2) )
      E_kin(i,j,2) = -sum( psi( indwvl_i(i):indwvl_f(i)-1,2) * tpsi_v(indwvl_i(j):indwvl_f(j)-1,1) ) + sum( psi( indwvl_i(i):indwvl_f(i)-1,1) * tpsi_v(indwvl_i(j):indwvl_f(j)-1,2) )
    end do ; end do
    do i=norb_start,norb_end  ; do j=1,  orbs%norb
      E_kin(i,j,1) = sum( psiv(indwvl_i(i):indwvl_f(i)-1,1) * tpsi_o(indwvl_i(j):indwvl_f(j)-1,1) ) + sum( psiv(indwvl_i(i):indwvl_f(i)-1,2) * tpsi_o(indwvl_i(j):indwvl_f(j)-1,2) )
      E_kin(i,j,2) = -sum( psiv(indwvl_i(i):indwvl_f(i)-1,2) * tpsi_o(indwvl_i(j):indwvl_f(j)-1,1) ) + sum( psiv(indwvl_i(i):indwvl_f(i)-1,1) * tpsi_o(indwvl_i(j):indwvl_f(j)-1,2) )
    end do ; end do
    do i=norb_start,norb_end  ; do j=norb_start,norb_end 
      E_kin(i,j,1) = sum( psiv(indwvl_i(i):indwvl_f(i)-1,1) * tpsi_v(indwvl_i(j):indwvl_f(j)-1,1) ) + sum( psiv(indwvl_i(i):indwvl_f(i)-1,2) * tpsi_v(indwvl_i(j):indwvl_f(j)-1,2) )
      E_kin(i,j,2) = -sum( psiv(indwvl_i(i):indwvl_f(i)-1,2) * tpsi_v(indwvl_i(j):indwvl_f(j)-1,1) ) + sum( psiv(indwvl_i(i):indwvl_f(i)-1,1) * tpsi_v(indwvl_i(j):indwvl_f(j)-1,2) )
    end do ; end do
 
    deallocate(tpsi_o, tpsi_v, indwvl_i, indwvl_f)
    call system('echo "kinetic calculation ... DONE"')

!-------------------------------------------------------------------------------------------------------------------------------------
!=====================================================================================================================================
!-------------------------------------------------------------------------------------------------------------------------------------

  call system('echo "potential calculation (local)"')
  !---------------------------------------------------------------------
  !          energy of the    local potential, unit in Hartree         -
  !---------------------------------------------------------------------
  write(124,*) 
  write(124,*) "================================================"
  write(124,*) " energy of the LOCAL potential, unit in Hartree "
  write(124,*) "       -- psir(i)*potential*psir(j) --       "
  write(124,*) "================================================"
  allocate(psir(  Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i, orbs%nspinor)) ; psir = 0._wp
  allocate(psir_i(Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i, orbs%nspinor)) ; psir_i = 0._wp
  allocate(psir_j(Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i, orbs%nspinor)) ; psir_j = 0._wp
  allocate(psirr( Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i,0:orbs%norb+orbsv%norb-1, orbs%nspinor)) ; psirr=0._wp

  call local_potential_dimensions(iproc, Lzd, orbs, xc, dpcom%ngatherarr(0,1))
  dict => dict_new()
  pkernel = pkernel_init(iproc,nproc,dict, atoms%astruct%geocode, &
            (/Lzd%Glr%d%n1i, Lzd%Glr%d%n2i, Lzd%Glr%d%n3i/), (/inputs%hx/2._gp,inputs%hy/2._gp,inputs%hz/2._gp/) )
  call dict_free(dict)
  call pkernel_set(pkernel,verbose=.false.)
  nullify(pot_ion) ; allocate(pot_ion(Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*dpcom%n3p)) ; pot_ion=0._gp
  nullify(rho_ion) ; allocate(rho_ion(Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*dpcom%n3p)) ; rho_ion=0._gp

  write(124,*) " call createIonicPotential"
  call createIonicPotential(iproc,(iproc==0),atoms, atoms%astruct%rxyz, inputs%elecfield, dpcom, pkernel, pot_ion, rho_ion, psoffset)

  write(124,*) " call full_local_potential"
  call full_local_potential(iproc,nproc, orbs,Lzd, 0,dpcom,xc, pot_ion, potential)       ! allocate the potential in the full box
! write(566,'(f15.8)') potential

  write(124,*) " call initialize_work_arrays_sumrho"
  call initialize_work_arrays_sumrho(    Lzd%Glr , .true., wisf)
  epot_sum_r = 0._dp
  epot_sum_i = 0._dp
  do i=1,orbtot
    do j=1,orbtot
      if(i .le. orbs%norb) then
        call daub_to_isf(Lzd%Glr,wisf,  psi( (i-1)          *(Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)+1,1) ,psir_i(:,1))
        call daub_to_isf(Lzd%Glr,wisf,  psi( (i-1)          *(Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)+1,2) ,psir_i(:,2))
      else
        call daub_to_isf(Lzd%Glr,wisf, psiv( (i-1-orbs%norb)*(Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)+1,1) ,psir_i(:,1))
        call daub_to_isf(Lzd%Glr,wisf, psiv( (i-1-orbs%norb)*(Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)+1,2) ,psir_i(:,2))
      end if

      if(j .le. orbs%norb) then
        call daub_to_isf(Lzd%Glr,wisf,  psi( (j-1)          *(Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)+1,1) ,psir_j(:,1))
        call daub_to_isf(Lzd%Glr,wisf,  psi( (j-1)          *(Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)+1,2) ,psir_j(:,2))
        if(i .eq. 1) psirr(:,j-1,1) = psir_j(:,1); psirr(:,j-1,2) = psir_j(:,2)
      else
        call daub_to_isf(Lzd%Glr,wisf, psiv( (j-1-orbs%norb)*(Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)+1,1) ,psir_j(:,1))
        call daub_to_isf(Lzd%Glr,wisf, psiv( (j-1-orbs%norb)*(Lzd%Glr%wfd%nvctr_c+7*Lzd%Glr%wfd%nvctr_f)+1,2) ,psir_j(:,2))
        if(i .eq. 1) psirr(:,j-1,1) = psir_j(:,1); psirr(:,j-1,2) = psir_j(:,2)
      end if

      do k=1,Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i
         epot_sum_r = epot_sum_r + psir_i(k,1)*potential(k)*psir_j(k,1) + psir_i(k,2)*potential(k)*psir_j(k,2)
         epot_sum_i = epot_sum_i + psir_i(k,1)*potential(k)*psir_j(k,2) - psir_i(k,2)*potential(k)*psir_j(k,1)
      end do
    ! write(124,"(f19.9)",advance = "NO") epot_sum  
      E_local(i,j,1) = epot_sum_r
      E_local(i,j,2) = epot_sum_i
      epot_sum_r = 0.0_gp
      epot_sum_i = 0.0_gp
    end do
  ! write(124,"(A)") " "
  end do ! write(124,*)

! do i=0,orbs%norb+orbsv%norb-1  
!   write(3000+i,*)   Lzd%Glr%d%n1i,   Lzd%Glr%d%n2i,   Lzd%Glr%d%n3i
!   write(3000+i,*) inputs%hx/2._gp, inputs%hy/2._gp, inputs%hz/2._gp
!   write(3000+i,'(f15.8)') psirr(:,i)  
! end do

  call free_full_potential(dpcom%mpi_env%nproc,0,xc,potential)
  if (nproc>1) call fmpi_allreduce(epot_sum_r,1,op=FMPI_SUM); if (nproc>1) call fmpi_allreduce(epot_sum_i,1,op=FMPI_SUM)
  deallocate(pot_ion, rho_ion, psir)
  call deallocate_work_arrays_sumrho(wisf)
  call system('echo "potential calculation (local) ... DONE"')

!-------------------------------------------------------------------------------------------------------------------------------------
!=====================================================================================================================================
!-------------------------------------------------------------------------------------------------------------------------------------

    call system('echo "potential calculation (nonlocal)"')
  !---------------------------------------------------------------------
  !          energy of the nonlocal potential, unit in Hartree         -
  !---------------------------------------------------------------------

    allocate(eprjo1(4,orbs%norb,orbs%nspinor),eprjv1(4,norb_start:norb_end,orbs%nspinor)) ; eprjo1 = 0._dp ; eprjv1 = 0._dp
    allocate(eprjo2(4,orbs%norb,orbs%nspinor),eprjv2(4,norb_start:norb_end,orbs%nspinor)) ; eprjo2 = 0._dp ; eprjv2 = 0._dp
    init_projectors_completely = .true. ; dry_run = 0 ; nwarnings = 0
    write(124,'("eprjo1 ",2i3," eprjv1 ",2i3,1h:,i3)') 4,orbs%norb, 4,norb_start,norb_end
    write(124,'("eprjo2 ",2i3," eprjv2 ",2i3,1h:,i3)') 4,orbs%norb, 4,norb_start,norb_end

    write(124,*) 
    write(124,*) "==============="
    write(124,*) " nonlocal term "
    write(124,*) "==============="
    !-------------------
    !  occupied orbitals
    !-------------------
    allocate(hpsi_r(orbdimocc)) ;allocate(hpsi_i(orbdimocc)); nn=0 
    hpsi_ptr_r => ob_ket_map(hpsi_r,psi_it_r) ; hpsi_ptr_i => ob_ket_map(hpsi_i,psi_it_i) ; 
    if(orbdimocc > 0) call f_zero(orbdimocc,hpsi_r(1)); if(orbdimocc > 0) call f_zero(orbdimocc,hpsi_i(1))
    call createProjectorsArrays(Lzd%Glr, atoms%astruct%rxyz,atoms,orbs,&
                                inputs%frmult,inputs%frmult,inputs%projection,dry_run,nlpsp_r,init_projectors_completely)
    call createProjectorsArrays(Lzd%Glr, atoms%astruct%rxyz,atoms,orbs,&
                                inputs%frmult,inputs%frmult,inputs%projection,dry_run,nlpsp_i,init_projectors_completely)
    call orbital_basis_associate(psi_ob_r, orbs=orbs, phis_wvl=psi(:,1), Lzd=Lzd, id='nonlocalham')
    call orbital_basis_associate(psi_ob_i, orbs=orbs, phis_wvl=psi(:,2), Lzd=Lzd, id='nonlocalham')
    if(associated(nlpsp_r%iagamma)) call f_zero(nlpsp_r%gamma_mmp) ; if(associated(nlpsp_i%iagamma)) call f_zero(nlpsp_i%gamma_mmp)
    psi_it_r = orbital_basis_iterator(psi_ob_r)
    psi_it_i = orbital_basis_iterator(psi_ob_i)
    loop_kpt: do while(ket_next_kpt(psi_it_r) .AND. ket_next_kpt(psi_it_i))
       loop_lr: do while(ket_next_locreg(psi_it_r, ikpt=psi_it_r%ikpt) .AND. ket_next_locreg(psi_it_i, ikpt=psi_it_i%ikpt))
          call DFT_PSP_projectors_iter_new(psp_it_r, nlpsp_r); call DFT_PSP_projector_iter_new(psp_it_i, nlpsp_i)
          loop_proj: do while (DFT_PSP_projectors_iter_next(psp_it_r, ilr=psi_it_r%ilr, lr=psi_it_r%lr, glr=Lzd%glr) .AND. DFT_PSP_projectors_iter_next(psp_it_i, ilr=psi_it_i%ilr, lr=psi_it_i%lr, glr=Lzd%glr))
             call DFT_PSP_projectors_iter_ensure(psp_it_r, psi_it_r%kpoint, 0, nwarnings, Lzd%Glr)
             call DFT_PSP_projectors_iter_ensure(psp_it_i, psi_it_i%kpoint, 0, nwarnings, Lzd%Glr)
             loop_psi_kpt: do while(ket_next(psi_it_r, ikpt=psi_it_r%ikpt, ilr=psi_it_r%ilr) .AND. ket_next(psi_it_i, ikpt=psi_it_i%ikpt, ilr=psi_it_i%ilr))
                call DFT_PSP_projectors_iter_apply(psp_it_r, psi_it_r, atoms, eproj, hpsi=hpsi_r, paw=paw_r)
                call DFT_PSP_projectors_iter_apply(psp_it_i, psi_it_i, atoms, eproj, hpsi=hpsi_i, paw=paw_i)
                eproj_sum = eproj_sum + psi_it_r%kwgt * psi_it_r%occup * eproj
                nn = nn + 1
                eprjo1(:,nn,1) = psp_it_r%parent%cproj
                eprjo1(:,nn,2) = psp_it_i%parent%cproj
                eprjo2(:,nn,1) = psp_it_r%parent%hcproj
                eprjo2(:,nn,2) = psp_it_i%parent%hcproj
                write(124,'(i2,3x,2(4f7.3,1x))',advance="no" )  nn, eprjo1(:,nn,1)! psp_it%parent%cproj
                write(124,'(   5x,2(4f7.3,1x))',advance="yes")      eprjo2(:,nn,1)! psp_it%parent%hcproj
             end do loop_psi_kpt
          end do loop_proj
       end do loop_lr
    end do loop_kpt
    call orbital_basis_release(psi_ob_r) ; deallocate(hpsi_r)
    call orbital_basis_release(psi_ob_i) ; deallocate(hpsi_i)
    !-------------------
    !  virtual orbitals
    !-------------------
    allocate(hpsi_r(orbdimvir)) ; allocate(hpsi_i(orbdimvir))! nn=0
    hpsi_ptr_r => ob_ket_map(hpsi_r,psi_it_r) ; if(orbdimvir > 0) call f_zero(orbdimvir,hpsi_r(1))
    hpsi_ptr_i => ob_ket_map(hpsi_i,psi_it_i) ; if(orbdimvir > 0) call f_zero(orbdimvir,hpsi_i(1))
    call createProjectorsArrays(Lzd%Glr, atoms%astruct%rxyz,atoms,orbsv,&
                                inputs%frmult,inputs%frmult,inputs%projection,dry_run,nlpsp_r,init_projectors_completely)
    call createProjectorsArrays(Lzd%Glr, atoms%astruct%rxyz,atoms,orbsv,&
                                inputs%frmult,inputs%frmult,inputs%projection,dry_run,nlpsp_i,init_projectors_completely)
    call orbital_basis_associate(psi_ob_r, orbs=orbsv, phis_wvl=psiv(:,1), Lzd=Lzd, id='nonlocalham')
    call orbital_basis_associate(psi_ob_i, orbs=orbsv, phis_wvl=psiv(:,2), Lzd=Lzd, id='nonlocalham')
    if(associated(nlpsp_r%iagamma)) call f_zero(nlpsp_r%gamma_mmp) ; if(associated(nlpsp_i%iagamma)) call f_zero(nlpsp_i%gamma_mmp)
    psi_it_r = orbital_basis_iterator(psi_ob_r)
    psi_it_i = orbital_basis_iterator(psi_ob_i)
    loop_kpt_v: do while(ket_next_kpt(psi_it_r) .AND. ket_next_kpt(psi_it_i))
       loop_lr_v: do while(ket_next_locreg(psi_it_r, ikpt=psi_it_r%ikpt) .AND. ket_next_locreg(psi_it_i, ikpt=psi_it_i%ikpt))
          call DFT_PSP_projectors_iter_new(psp_it_r, nlpsp_r)
          call DFT_PSP_projectors_iter_new(psp_it_i, nlpsp_i)
          loop_proj_v: do while (DFT_PSP_projectors_iter_next(psp_it_r, ilr=psi_it_r%ilr, lr=psi_it_r%lr, glr=Lzd%glr) .AND. DFT_PSP_projectors_iter_next(psp_it_i, ilr=psi_it_i%ilr, lr=psi_it_i%lr, glr=Lzd%glr))
             call DFT_PSP_projectors_iter_ensure(psp_it_r, psi_it_r%kpoint, 0, nwarnings, Lzd%Glr)
             call DFT_PSP_projectors_iter_ensure(psp_it_i, psi_it_i%kpoint, 0, nwarnings, Lzd%Glr)
             loop_psi_kpt_v: do while(ket_next(psi_it_r, ikpt=psi_it_r%ikpt, ilr=psi_it_r%ilr) .AND. ket_next(psi_it_i, ikpt=psi_it_i%ikpt, ilr=psi_it_i%ilr))
                call DFT_PSP_projectors_iter_apply(psp_it_r, psi_it_r, atoms, eproj, hpsi=hpsi_r, paw=paw_r)
                call DFT_PSP_projectors_iter_apply(psp_it_i, psi_it_i, atoms, eproj, hpsi=hpsi_i, paw=paw_i)
                eproj_sum = eproj_sum + psi_it_r%kwgt * psi_it_r%occup * eproj
                nn = nn + 1
                eprjv1(:,nn,1) = psp_it_r%parent%cproj
                eprjv2(:,nn,1) = psp_it_r%parent%hcproj
                eprjv1(:,nn,2) = psp_it_i%parent%cproj
                eprjv2(:,nn,2) = psp_it_i%parent%hcproj
                write(124,'(i2,3x,2(4f7.3,1x))',advance="no" )  nn, eprjv1(:,nn,1)! psp_it%parent%cproj
                write(124,'(   5x,2(4f7.3,1x))',advance="yes")      eprjv2(:,nn,1)! psp_it%parent%hcproj
          !     if(nn .eq. norb_end) go to 99
             end do loop_psi_kpt_v
          end do loop_proj_v
       end do loop_lr_v
    end do loop_kpt_v
99  call orbital_basis_release(psi_ob_r) ; deallocate(hpsi_r) ; call orbital_basis_release(psi_ob_i) ; deallocate(hpsi_i)
 
    write(124,*) ; write(124,*) "---- eproj1 (psp_it%parent%cproj)  ----"
    do i=1,orbs%norb ; write(124,'(i3)',advance="no") i  
      do j=1,4 ; write(124,"(f12.7),2x,(f12.7)",advance="no") eprjo1(j,i,1),eprjo1(j,i,2) ; end do ; write(124,*) ; end do
    do i=norb_start,norb_end ; write(124,'(i3)',advance="no") i  
      do j=1,4 ; write(124,"(f12.7),2x,(f12.7)",advance="no") eprjv1(j,i,1),eprjo1(j,i,2) ; end do ; write(124,*) ; end do

    write(124,*) ; write(124,*) "---- eproj2 (psp_it%parent%hcproj) ----" 
    do i=1,orbs%norb ; write(124,'(i3)',advance="no") i
     do j=1,4 ; write(124,"(f12.7),2x,(f12.7)",advance="no") eprjo2(j,i,1),eprjo2(j,i,2) ; end do ; write(124,*) ; end do
    do i=norb_start,norb_end ; write(124,'(i3)',advance="no") i  
     do j=1,4 ; write(124,"(f12.7),2x,(f12.7)",advance="no") eprjv2(j,i,1),eprjo2(j,i,2) ; end do ; write(124,*) ; end do
    write(124,*)

    do i=1,orbs%norb          ; do j=1,orbs%norb 
      E_nonlocal(i,j,1) = sum(eprjo1(:,i,1)*eprjo2(:,j,1))+sum(eprjo1(:,i,2)*eprjo2(:,j,2)) ; E_nonlocal(i,j,2) = sum(eprjo1(:,i,1)*eprjo2(:,j,2))-sum(eprjo1(:,i,2)*eprjo2(:,j,1))
      write(124,787) i,eprjo1(:,i,1), j,eprjo2(:,j,1) ; write(124,787) i,eprjo1(:,i,2),j,eprjo2(:,j,2)
    end do ; end do ; write(124,*) 
    do i=1,orbs%norb          ; do j=norb_start,norb_end 
      E_nonlocal(i,j,1) = sum(eprjo1(:,i,1)*eprjv2(:,j,1))+sum(eprjo1(:,i,2)*eprjv2(:,j,2)) ; E_nonlocal(i,j,2) = sum(eprjo1(:,i,1)*eprjv2(:,j,2))-sum(eprjo1(:,i,2)*eprjv2(:,j,1))
      write(124,787) i,eprjo1(:,i,1), j,eprjv2(:,j,1) ; write(124,787) i,eprjo1(:,i,2),j,eprjv2(:,j,2)
    end do ; end do ; write(124,*)
    do i=norb_start,norb_end  ; do j=1,orbs%norb    
      E_nonlocal(i,j,1) = sum(eprjv1(:,i,1)*eprjo2(:,j,1))+sum(eprjv1(:,i,2)*eprjo2(:,j,2)) ; E_nonlocal(i,j,2) = sum(eprjv1(:,i,1)*eprjo2(:,j,2))-sum(eprjv1(:,i,2)*eprjo2(:,j,1))
      write(124,787) i,eprjv1(:,i,1), j,eprjo2(:,j,1) ; write(124,787) i,eprjv1(:,i,2),j,eprjo2(:,j,2)
    end do ; end do ; write(124,*)
    do i=norb_start,norb_end  ; do j=norb_start,norb_end 
      E_nonlocal(i,j,1) = sum(eprjv1(:,i,1)*eprjv2(:,j,1))+sum(eprjv1(:,i,2)*eprjv2(:,j,2)) ; E_nonlocal(i,j,2) = sum(eprjv1(:,i,1)*eprjv2(:,j,2))-sum(eprjv1(:,i,2)*eprjv2(:,j,1))
      write(124,787) i,eprjv1(:,i,1), j,eprjv2(:,j,1) ; write(124,787) i,eprjv1(:,i,2),j,eprjv2(:,j,2)
    end do ; end do ; write(124,*)
    do i=1,orbtot ; write(124,34) (E_nonlocal(j,i,1),j=1,orbtot) ; end do 
    do i=1,orbtot ; write(124,34) (E_nonlocal(j,i,2),j=1,orbtot) ; end do 
787 format(i2,1x,4f10.6,2x,i2,1x,4f10.6)
34  format(20f10.6)
    
  deallocate(eprjo1, eprjo2, eprjv1, eprjv2)
  deallocate(rxyz_old, psi, psir_i, psir_j)
  close(6)
  call system('echo "potential calculation (nonlocal) ... DONE"')

!-------------------------------------------------------------------------------------------------------------------------------------
!=====================================================================================================================================
!-------------------------------------------------------------------------------------------------------------------------------------

  !----------------
  !   data output 
  !----------------
  write(6,*) 
  write(6,'("number of occupied orbitals - orbs%norb  : ",5i4)') orbs%norb
  write(6,'("number of  virtual orbitals - orbsv%norb : ",5i4)') orbsv%norb
  write(6,'("    input parameter norbv - inputs%norbv : ",5i4)') inputs%norbv
  write(6,'("number of    total orbitals -     orbtot : ",5i4)') orbtot
  write(6,'("                            orbs%nspinor : ",5i4)') orbs%nspinor
  write(6,'("                                   mspin : ",5i4)') mspin
  write(6,'("                              OMP number : ",5i4)') OMP_get_thread_num(), OMP_get_num_threads(), OMP_get_max_threads()
  write(6,'("                          orbs%npsidim_orbs : ", i8)') orbs%npsidim_orbs
  write(6,'(" dimension of occupied orbitals (orbdimocc) : ", i8)') orbdimocc
  write(6,'(" dimension of  virtual orbitals (orbdimvir) : ", i8)') orbdimvir

  allocate( output1(0:orbtot-1,0:orbtot-1,orbs%nspinor) ) ; output1=0._dp
  allocate( output2(0:orbtot-1,0:orbtot-1,orbs%nspinor) ) ; output2=0._dp
  allocate( output3(0:orbtot-1,0:orbtot-1,orbs%nspinor) ) ; output3=0._dp
  if(orbs%nspin .eq. 2) then
    do i=1,orbtot ; do j=1,orbtot
     ii = re_index1(i-1,orbtot,orbocc)
     jj = re_index1(j-1,orbtot,orbocc)
     output1(ii,jj,:) =      E_kin(i,j,:)
     output2(ii,jj,:) =    E_local(i,j,:)
     output3(ii,jj,:) = E_nonlocal(i,j,:)
!    write(4567,'(2i3,2x,2i3,2x,f15.8)') i-1,j-1, ii,jj, E_kin(ii,jj)+E_local(ii,jj)+E_nonlocal(ii,jj)
    end do ; end do
   
    write(6,*)  
    write(6,*) ; write(6,*) " =====  kinetic energy: psi(i)*tpsi(j) =====",orbtot
    do j=0,orbtot-1; do i=0,orbtot-1
      write(6,1,advance="NO") output1(i,j,1) ; write(6,1,advance="NO") output1(i,j,2) ; if(i .eq. orbocc/2) write(6,"(A)",advance="NO") "  " ; end do
      write(6,"(A)") " " ; if(j .eq. orbocc/2) write(6,"(A)") "  " ; end do

    write(6,*) ; write(6,*) " =====  local potential: psir(i)*potential(ij)*psir(j)  ===== "
    do j=0,orbtot-1; do i=0,orbtot-1
      write(6,1,advance="NO") output2(i,j,1) ; write(6,1,advance="NO") output2(i,j,2) ; if(i .eq. orbocc/2) write(6,"(A)",advance="NO") "  " ; end do
      write(6,"(A)") " " ; if(j .eq. orbocc/2) write(6,"(A)") "  " ; end do
    write(6,*) ; write(6,*) " =====  nonlocal potential: psir(i)*hpsir(j)  ===== "

    do j=0,orbtot-1; do i=0,orbtot-1 
      write(6,1,advance="NO") output3(i,j,1) ; write(6,1,advance="NO") output3(i,j,2) ; if(i .eq. orbocc/2) write(6,"(A)",advance="NO") "  " ; end do
      write(6,"(A)") " " ; if(j .eq. orbocc/2) write(6,"(A)") "  " ; end do

    write(6,*) ; write(6,*) " =====  hpq:  sum(kinetic, V_local, V_nonlocal)  ===== "
    do j=0,orbtot-1; do i=0,orbtot-1 
      write(6,1,advance="NO") output1(i,j,1) + output2(i,j,1) + output3(i,j,1)  
      write(6,1,advance="NO") output1(i,j,2) + output2(i,j,2) + output3(i,j,2)  
      if(i .eq. orbocc/2) write(6,"(A)",advance="NO") "  "
      end do ; write(6,"(A)") " " ; if(j .eq. orbocc/2) write(6,"(A)") "  " ; end do

  else

    do i=1,orbtot ; do j=1,orbtot
     output1(i,j,:) =      E_kin(i,j,:)
     output2(i,j,:) =    E_local(i,j,:)
     output3(i,j,:) = E_nonlocal(i,j,:)
    end do ; end do
   
    write(6,*)  
    write(6,*) ; write(6,*) " =====  kinetic energy: psi(i)*tpsi(j) ====="
    do j=1,orbtot; do i=1,orbtot; write(6,1,advance="NO") output1(i,j,1);write(6,1,advance="NO") output1(i,j,2); if(i .eq. orbs%norb) write(6,"(A)",advance="NO") "  " ; end do
     write(6,"(A)") " " ; if(j .eq. orbs%norb) write(6,"(A)") "  " ; end do 

    write(6,*) ; write(6,*) " =====  local potential: psir(i)*potential(ij)*psir(j)  ===== "
    do j=1,orbtot; do i=1,orbtot; write(6,1,advance="NO") output2(i,j,1); write(6,1,advance="NO") output2(i,j,2); if(i .eq. orbs%norb) write(6,"(A)",advance="NO") "  " ; end do
      write(6,"(A)") " " ; if(j .eq. orbs%norb) write(6,"(A)") "  " ; end do

    write(6,*) ; write(6,*) " =====  nonlocal potential: psir(i)*hpsir(j)  ===== "
    do j=1,orbtot; do i=1,orbtot; write(6,1,advance="NO") output3(i,j,1); write(6,1,advance="NO") output3(i,j,2); if(i .eq. orbs%norb) write(6,"(A)",advance="NO") "  " ; end do
      write(6,"(A)") " " ; if(j .eq. orbs%norb) write(6,"(A)") "  " ; end do

    write(6,*) ; write(6,*) " =====  hpq:  sum(kinetic, V_local, V_nonlocal)  ===== "
    do j=1,orbtot; do i=1,orbtot; write(6,1,advance="NO") output1(i,j,1) + output2(i,j,1) + output3(i,j,1) ;  write(6,1,advance="NO") output1(i,j,2) + output2(i,j,2) + output3(i,j,2) ; 
      if(i .eq. orbs%norb) write(6,"(A)",advance="NO") "  "
      end do ; write(6,"(A)") " " ; if(j .eq. orbs%norb) write(6,"(A)") "  " ; end do
    
  end if
! -------------------------------------------------------------------------------------------

  write(6,*)

! ----------------------------  H2, UHF  ----------------------------
! i=0 ; j=0 ;write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i=0 ; j=2 ;write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i=0 ; j=4 ;write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i=0 ; j=6 ;write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i=1 ; j=1 ;write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i=1 ; j=3 ;write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i=1 ; j=5 ;write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i=1 ; j=7 ;write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i=2 ; j=2 ;write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i=2 ; j=4 ;write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i=2 ; j=6 ;write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i=3 ; j=3 ;write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i=3 ; j=5 ;write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i=3 ; j=7 ;write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i=4 ; j=4 ;write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i=4 ; j=6 ;write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i=5 ; j=5 ;write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i=5 ; j=7 ;write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i=6 ; j=6 ;write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i=7 ; j=7 ;write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)

! ----------------------------  H2, RHF  ----------------------------
! i=1 ; j=1 ;write(555,2) (i-1)*2,(j-1)*2, output1(i,j)+output2(i,j)+output3(i,j)
! i=2 ; j=2 ;write(555,2) (i-1)*2,(j-1)*2, output1(i,j)+output2(i,j)+output3(i,j)

! ----------------------------  LiH, RHF  ----------------------------
! i=1 ; j=1 ;write(555,2) (i-1)*2,(j-1)*2, output1(i,j)+output2(i,j)+output3(i,j)
! i=3 ; j=3 ;write(555,2) (i-1)*2,(j-1)*2, output1(i,j)+output2(i,j)+output3(i,j)
! i=2 ; j=2 ;write(555,2) (i-1)*2,(j-1)*2, output1(i,j)+output2(i,j)+output3(i,j)
! i=5 ; j=5 ;write(555,2) (i-1)*2,(j-1)*2, output1(i,j)+output2(i,j)+output3(i,j)
! i=5 ; j=2 ;write(555,2) (i-1)*2,(j-1)*2, output1(i,j)+output2(i,j)+output3(i,j)
! i=1 ; j=2 ;write(555,2) (i-1)*2,(j-1)*2, output1(i,j)+output2(i,j)+output3(i,j)
! i=1 ; j=5 ;write(555,2) (i-1)*2,(j-1)*2, output1(i,j)+output2(i,j)+output3(i,j)

! ----------------------------  LiH, UHF, 1e+1e  ----------------------------
! i=0 ;j=0 ; write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i=0 ;j=2 ; write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i=0 ;j=4 ; write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i=0 ;j=6 ; write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i=0 ;j=8 ; write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i=1 ;j=1 ; write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i=1 ;j=3 ; write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i=1 ;j=5 ; write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i=1 ;j=7 ; write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i=1 ;j=9 ; write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i=2 ;j=2 ; write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i=2 ;j=4 ; write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i=2 ;j=6 ; write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i=2 ;j=8 ; write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i=3 ;j=3 ; write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i=3 ;j=5 ; write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i=3 ;j=7 ; write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i=3 ;j=9 ; write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i=4 ;j=4 ; write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i=4 ;j=6 ; write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i=4 ;j=8 ; write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i=5 ;j=5 ; write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i=5 ;j=7 ; write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i=5 ;j=9 ; write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i=6 ;j=6 ; write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i=6 ;j=8 ; write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i=7 ;j=7 ; write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i=7 ;j=9 ; write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i=8 ;j=8 ; write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i=9 ;j=9 ; write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)

! ----------------------------  LiH, UHF, 3e+1e  ----------------------------
! i= 0 ; j= 0 ;write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i= 0 ; j= 2 ;write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i= 0 ; j= 4 ;write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i= 0 ; j= 6 ;write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i= 0 ; j= 8 ;write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i= 0 ; j=10 ;write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i= 1 ; j= 1 ;write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i= 1 ; j= 3 ;write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i= 1 ; j= 5 ;write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i= 1 ; j= 7 ;write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i= 1 ; j= 9 ;write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i= 1 ; j=11 ;write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i= 2 ; j= 2 ;write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i= 2 ; j= 4 ;write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i= 2 ; j= 6 ;write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i= 2 ; j= 8 ;write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i= 2 ; j=10 ;write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i= 3 ; j= 3 ;write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i= 3 ; j= 5 ;write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i= 3 ; j= 7 ;write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i= 3 ; j= 9 ;write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i= 3 ; j=11 ;write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i= 4 ; j= 4 ;write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i= 4 ; j= 6 ;write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i= 4 ; j= 8 ;write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i= 4 ; j=10 ;write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i= 5 ; j= 5 ;write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i= 5 ; j= 7 ;write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i= 5 ; j= 9 ;write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i= 5 ; j=11 ;write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i= 6 ; j= 6 ;write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i= 6 ; j= 8 ;write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i= 6 ; j=10 ;write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i= 7 ; j= 7 ;write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i= 7 ; j= 9 ;write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i= 7 ; j=11 ;write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i= 8 ; j= 8 ;write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i= 8 ; j=10 ;write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i= 9 ; j= 9 ;write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i= 9 ; j=11 ;write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i=10 ; j=10 ;write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)
! i=11 ; j=11 ;write(555,2) i,j, output1(i,j)+output2(i,j)+output3(i,j)

! ----------------------------   H2O, RHF, new  ----------------------------
  i=1 ; j=1 ;write(555,2) (i-1)*2,(j-1)*2, output1(i,j,1)+output2(i,j,1)+output3(i,j,1)
  i=1 ; j=2 ;write(555,2) (i-1)*2,(j-1)*2, output1(i,j,1)+output2(i,j,1)+output3(i,j,1)
  i=1 ; j=3 ;write(555,2) (i-1)*2,(j-1)*2, output1(i,j,1)+output2(i,j,1)+output3(i,j,1)
  i=1 ; j=4 ;write(555,2) (i-1)*2,(j-1)*2, output1(i,j,1)+output2(i,j,1)+output3(i,j,1)
  i=1 ; j=5 ;write(555,2) (i-1)*2,(j-1)*2, output1(i,j,1)+output2(i,j,1)+output3(i,j,1)
  i=1 ; j=6 ;write(555,2) (i-1)*2,(j-1)*2, output1(i,j,1)+output2(i,j,1)+output3(i,j,1)
  i=2 ; j=2 ;write(555,2) (i-1)*2,(j-1)*2, output1(i,j,1)+output2(i,j,1)+output3(i,j,1)
  i=2 ; j=3 ;write(555,2) (i-1)*2,(j-1)*2, output1(i,j,1)+output2(i,j,1)+output3(i,j,1)
  i=2 ; j=4 ;write(555,2) (i-1)*2,(j-1)*2, output1(i,j,1)+output2(i,j,1)+output3(i,j,1)
  i=2 ; j=5 ;write(555,2) (i-1)*2,(j-1)*2, output1(i,j,1)+output2(i,j,1)+output3(i,j,1)
  i=2 ; j=6 ;write(555,2) (i-1)*2,(j-1)*2, output1(i,j,1)+output2(i,j,1)+output3(i,j,1)
  i=3 ; j=3 ;write(555,2) (i-1)*2,(j-1)*2, output1(i,j,1)+output2(i,j,1)+output3(i,j,1)
  i=3 ; j=4 ;write(555,2) (i-1)*2,(j-1)*2, output1(i,j,1)+output2(i,j,1)+output3(i,j,1)
  i=3 ; j=5 ;write(555,2) (i-1)*2,(j-1)*2, output1(i,j,1)+output2(i,j,1)+output3(i,j,1)
  i=3 ; j=6 ;write(555,2) (i-1)*2,(j-1)*2, output1(i,j,1)+output2(i,j,1)+output3(i,j,1)
  i=4 ; j=4 ;write(555,2) (i-1)*2,(j-1)*2, output1(i,j,1)+output2(i,j,1)+output3(i,j,1)
  i=4 ; j=5 ;write(555,2) (i-1)*2,(j-1)*2, output1(i,j,1)+output2(i,j,1)+output3(i,j,1)
  i=4 ; j=6 ;write(555,2) (i-1)*2,(j-1)*2, output1(i,j,1)+output2(i,j,1)+output3(i,j,1)
  i=5 ; j=5 ;write(555,2) (i-1)*2,(j-1)*2, output1(i,j,1)+output2(i,j,1)+output3(i,j,1)
  i=5 ; j=6 ;write(555,2) (i-1)*2,(j-1)*2, output1(i,j,1)+output2(i,j,1)+output3(i,j,1)
  i=6 ; j=6 ;write(555,2) (i-1)*2,(j-1)*2, output1(i,j,1)+output2(i,j,1)+output3(i,j,1)
  i=1 ; j=1 ;write(555,2) (i-1)*2,(j-1)*2, output1(i,j,2)+output2(i,j,2)+output3(i,j,2)
  i=1 ; j=2 ;write(555,2) (i-1)*2,(j-1)*2, output1(i,j,2)+output2(i,j,2)+output3(i,j,2)
  i=1 ; j=3 ;write(555,2) (i-1)*2,(j-1)*2, output1(i,j,2)+output2(i,j,2)+output3(i,j,2)
  i=1 ; j=4 ;write(555,2) (i-1)*2,(j-1)*2, output1(i,j,2)+output2(i,j,2)+output3(i,j,2)
  i=1 ; j=5 ;write(555,2) (i-1)*2,(j-1)*2, output1(i,j,2)+output2(i,j,2)+output3(i,j,2)
  i=1 ; j=6 ;write(555,2) (i-1)*2,(j-1)*2, output1(i,j,2)+output2(i,j,2)+output3(i,j,2)
  i=2 ; j=2 ;write(555,2) (i-1)*2,(j-1)*2, output1(i,j,2)+output2(i,j,2)+output3(i,j,2)
  i=2 ; j=3 ;write(555,2) (i-1)*2,(j-1)*2, output1(i,j,2)+output2(i,j,2)+output3(i,j,2)
  i=2 ; j=4 ;write(555,2) (i-1)*2,(j-1)*2, output1(i,j,2)+output2(i,j,2)+output3(i,j,2)
  i=2 ; j=5 ;write(555,2) (i-1)*2,(j-1)*2, output1(i,j,2)+output2(i,j,2)+output3(i,j,2)
  i=2 ; j=6 ;write(555,2) (i-1)*2,(j-1)*2, output1(i,j,2)+output2(i,j,2)+output3(i,j,2)
  i=3 ; j=3 ;write(555,2) (i-1)*2,(j-1)*2, output1(i,j,2)+output2(i,j,2)+output3(i,j,2)
  i=3 ; j=4 ;write(555,2) (i-1)*2,(j-1)*2, output1(i,j,2)+output2(i,j,2)+output3(i,j,2)
  i=3 ; j=5 ;write(555,2) (i-1)*2,(j-1)*2, output1(i,j,2)+output2(i,j,2)+output3(i,j,2)
  i=3 ; j=6 ;write(555,2) (i-1)*2,(j-1)*2, output1(i,j,2)+output2(i,j,2)+output3(i,j,2)
  i=4 ; j=4 ;write(555,2) (i-1)*2,(j-1)*2, output1(i,j,2)+output2(i,j,2)+output3(i,j,2)
  i=4 ; j=5 ;write(555,2) (i-1)*2,(j-1)*2, output1(i,j,2)+output2(i,j,2)+output3(i,j,2)
  i=4 ; j=6 ;write(555,2) (i-1)*2,(j-1)*2, output1(i,j,2)+output2(i,j,2)+output3(i,j,2)
  i=5 ; j=5 ;write(555,2) (i-1)*2,(j-1)*2, output1(i,j,2)+output2(i,j,2)+output3(i,j,2)
  i=5 ; j=6 ;write(555,2) (i-1)*2,(j-1)*2, output1(i,j,2)+output2(i,j,2)+output3(i,j,2)
  i=6 ; j=6 ;write(555,2) (i-1)*2,(j-1)*2, output1(i,j,2)+output2(i,j,2)+output3(i,j,2)

 1 format(f10.6)
 2 format(5x,2i3,3x,f20.16,6x,2i3)
  deallocate(E_local, E_nonlocal, E_kin)
  deallocate(output1, output2, output3)

!-------------------------------------------------------------------------------------------------------------------------------------
!=====================================================================================================================================
!-------------------------------------------------------------------------------------------------------------------------------------

! write(6,*) "continue ............... ?" ; read(5,*) con
  con="y"
! con="n"

  if(con .eq. "y" .or. con .eq. "Y") then
   continue
  else
    open(  6,file="toy_model.out",status="old")
    !----------------------
    ! Free allocated space.
    !----------------------
    call deallocate_comms(comms)
    call deallocate_locreg_descriptors(Lzd%Glr)
    call deallocate_Lzd_except_Glr(Lzd)
    call deallocate_orbs(orbs)
    call deallocate_atoms_data(atoms)
    call xc_end(xc)
    call dpbox_free(dpcom)
    call pkernel_free(pkernel)
    call free_input_variables(inputs)
    
    !--------------------------------------
    !wait all processes before finalisation
    !--------------------------------------
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_FINALIZE(ierr)
    call f_lib_finalize()
    close(6)
    write(6,*) "program stop!" ; stop
  end if

!-------------------------------------------------------------------------------------------------------------------------------------

  mm = Lzd%Glr%d%n1 *Lzd%Glr%d%n2 *Lzd%Glr%d%n3
  nn = Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i
  hxh = inputs%hx/2._gp
  hyh = inputs%hy/2._gp
  hzh = inputs%hz/2._gp
  write(6,*) " ---------------  parameters check  ---------------------"
  write(6,'("    Number of orbitals: ",i10)') orbs%norb
  write(6,'("    Number of virtuals: ",i10)') orbsv%norb
  write(6,'("            n1, n2, n3: ",3i7,i10)')  Lzd%Glr%d%n1,  Lzd%Glr%d%n2, Lzd%Glr%d%n3,  mm
  write(6,'("            hx, hy, hz: ",3f7.3  )') inputs%hx, inputs%hy, inputs%hz
  write(6,'("         n1i, n2i, n3i: ",3i7,i10)') Lzd%Glr%d%n1i, Lzd%Glr%d%n2i, Lzd%Glr%d%n3i, nn
  write(6,'("            hx, hy, hz: ",3f7.3  )') hxh, hyh, hzh
  write(6,*) " --------------------------------------------------------"

!-------------------------------------------------------------------------------------------------------------------------------------
 
  allocate(nii(nn),njj(nn),nkk(nn)) ; mm=0 
  do k=1,Lzd%Glr%d%n3i ; do j=1,Lzd%Glr%d%n2i ; do i=1,Lzd%Glr%d%n1i
    mm = mm+1
    nii(mm) = real(i)*hxh
    njj(mm) = real(j)*hyh
    nkk(mm) = real(k)*hzh
  end do ; end do ; end do
 
!-------------------------------------------------------------------------------------------------------------------------------------
 
     call system("grep --color=no EH: log-*.yaml|tail -n 1|wc|awk '{print $2}' > .11")
     open(99,file='.11') ; read(99,*) nv ; close(99)
     open(99,file='.command') ; EH=0._gp  
     if(nv .eq. 6) write(99,33)
     if(nv .eq. 7) write(99,32)
     call system("chmod 755 ./.command") ; call system("sh .command") ; close(99)
     open(99,file='.11') ; read(99,*) EH ; close(99) ; ee1 = 2._gp * EH
     call system('rm -f .11 .command') 
  33 format(73hgrep --color=no EH: log-*.yaml|tail -n 1|awk '{printf "%15.9f",$2}' > .11)
  32 format(73hgrep --color=no EH: log-*.yaml|tail -n 1|awk '{printf "%15.9f",$7}' > .11)
     if(EH .eq. 0._gp) then
       write(6,*) "EH is missing! program stop ..." ; stop
     end if
 
     allocate(tho(Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i,orbs%nspinor)); ee3_r = 0._gp ; ee3_i= 0._gp ; tho = 0._gp
     do i=0,orbs%norb-1 ; tho(:,1) = tho(:,1) + (psirr(:,i,1)*psirr(:,i,1)+psirr(:,i,2)*psirr(:,i,2)) ;tho(:,2) = tho(:,2) + (psirr(:,i,1)*psirr(:,i,2)-psirr(:,i,2)*psirr(:,i,1)) ; end do
     ee3_r = (cspin**2) * sum(tho(:,1)*tho(:,1))+sum(tho(:,2)*tho(:,2))
     ee3_i = (cspin**2) * sum(tho(:,1)*tho(:,2))-sum(tho(:,2)*tho(:,2))
     
     cri=1e-10 ; nn_new=0
     call tho_diet1(tho(:,1),nn,cri,cspin, nn_new)
     call tho_diet1(tho(:,2),nn,cri,cspin, nn_new)
     allocate(tho_new1(nn_new,orbs%nspinor),thoi1(nn_new,orbs%nspinor)) ; tho_new1=0.d0 ; thoi1=0
     call tho_diet2(tho(:,1),nn,cri,nn_new,cspin, tho_new1(:,1),thoi1(:,1))
     call tho_diet2(tho(:,2),nn,cri,nn_new,cspin, tho_new1(:,2),thoi1(:,2))
     
     !real part(11)
     !$OMP PARALLEL DO PRIVATE(rx,ry,rz,rr,ii,jj, hpqrs_t1,hpqrs_t2,hpqrs_t3,hpqrs_t4) &
     !$OMP REDUCTION(+:hpqrs_r)
     do i=1,nn_new,1
       hpqrs_t1 = cspin * tho_new1(i,1) ; ii = thoi1(i,1)
       do j=1,nn_new,1
         hpqrs_t3 = cspin * tho_new1(j,1) ; jj = thoi1(j,1)
         if(ii .eq. jj) cycle
         rx = nii(ii) - nii(jj)
         ry = njj(ii) - njj(jj)
         rz = nkk(ii) - nkk(jj)
         rr = sqrt(rx**2 + ry**2 + rz**2)
         hpqrs_t2 = 1.d0/rr
         hpqrs_t4 = hpqrs_t1 * hpqrs_t2 * hpqrs_t3
         hpqrs_r = hpqrs_r + hpqrs_t4
       end do
     end do
     !$OMP END PARALLEL DO

     !real part(22)
     !$OMP PARALLEL DO PRIVATE(rx,ry,rz,rr,ii,jj, hpqrs_t1,hpqrs_t2,hpqrs_t3,hpqrs_t4) &
     !$OMP REDUCTION(+:hpqrs_r)
     do i=1,nn_new,1
       hpqrs_t1 = cspin * tho_new1(i,2) ; ii = thoi1(i,2)
       do j=1,nn_new,1
         hpqrs_t3 = cspin * tho_new1(j,2) ; jj = thoi1(j,2)
         if(ii .eq. jj) cycle
         rx = nii(ii) - nii(jj)
         ry = njj(ii) - njj(jj)
         rz = nkk(ii) - nkk(jj)
         rr = sqrt(rx**2 + ry**2 + rz**2)
         hpqrs_t2 = 1.d0/rr
         hpqrs_t4 = hpqrs_t1 * hpqrs_t2 * hpqrs_t3
         hpqrs_r = hpqrs_r + hpqrs_t4
       end do
     end do
     !$OMP END PARALLEL DO


     !imaginary part(12)
     !$OMP PARALLEL DO PRIVATE(rx,ry,rz,rr,ii,jj, hpqrs_t1,hpqrs_t2,hpqrs_t3,hpqrs_t4) &
     !$OMP REDUCTION(+:hpqrs_i)
     do i=1,nn_new,1
       hpqrs_t1 = cspin * tho_new1(i,1) ; ii = thoi1(i,1)
       do j=1,nn_new,1
         hpqrs_t3 = cspin * tho_new1(j,2) ; jj = thoi1(j,2)
         if(ii .eq. jj) cycle
         rx = nii(ii) - nii(jj)
         ry = njj(ii) - njj(jj)
         rz = nkk(ii) - nkk(jj)
         rr = sqrt(rx**2 + ry**2 + rz**2)
         hpqrs_t2 = 1.d0/rr
         hpqrs_t4 = hpqrs_t1 * hpqrs_t2 * hpqrs_t3
         hpqrs_i = hpqrs_i + hpqrs_t4
       end do
     end do
     !$OMP END PARALLEL DO

     !imaginary part(21)
     !$OMP PARALLEL DO PRIVATE(rx,ry,rz,rr,ii,jj,hpqrs_t1,hpqrs_t2,hpqrs_t3,hpqrs_t4) &
     !$OMP REDUCTION(+:hpqrs_i)
     do i=1,nn_new,1
       hpqrs_t1 = cspin * tho_new1(i,2) ; ii = thoi1(i,2)
       do j=1,nn_new,1
         hpqrs_t3 = cspin * tho_new1(j,1) ; jj = thoi1(j,1)
         if(ii .eq. jj) cycle
         rx = nii(ii) - nii(jj)
         ry = njj(ii) - njj(jj)
         rz = nkk(ii) - nkk(jj)
         rr = sqrt(rx**2 + ry**2 + rz**2)
         hpqrs_t2 = 1.d0/rr
         hpqrs_t4 = hpqrs_t1 * hpqrs_t2 * hpqrs_t3
         hpqrs_i = hpqrs_i - hpqrs_t4
       end do
     end do
     !$OMP END PARALLEL DO

     deallocate(tho_new1,thoi1)
     ee2_r = hpqrs_r*1.d0
     ee2_i = hpqrs_i*1.d0

     r_r_r = (ee1 - ee2_r) / ee3_r
     r_r_i = (ee1 - ee2_i) / ee3_i
     write( 6 ,'("EH: ",5f12.6,"    1/(ra-rb) at ra=rb : ",2f12.5)') ee1, ee2_r, ee2_i, ee3_r, ee3_i, r_r_r, r_r_i
     write(124,'("EH: ",5f12.6,"    1/(ra-rb) at ra=rb : ",2f12.5)') ee1, ee2_r, ee2_i, ee3_r, ee3_i, r_r_r, r_r_i
     deallocate(tho)
     end do
!-------------------------------------------------------------------------------------------------------------------------------------
    !exp(ik'r')
    nrpqrs=0 ; open(111,file="input-r")
    do ; read(111,*,iostat=istat) tmp ; nrpqrs = nrpqrs + 1 ; if(istat .ne. 0) exit ; end do
    nrpqrs = nrpqrs - 1 ; close(111) ; open(111,file="input-r")
    !unit cell data
    hpqrs_r=0._gp ; hpqrs_i ;  nhpqrs=0 ; open(112,file="input-hpqrs")
    do ; read(112,*,iostat=istat) tmp ; nhpqrs = nhpqrs + 1 ; if(istat .ne. 0) exit ; end do
    nhpqrs = nhpqrs - 1 ; close(112) ; open(112,file="input-hpqrs")
 
    allocate(tho1(Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i,orbs%nspinor))
    allocate(tho2(Lzd%Glr%d%n1i*Lzd%Glr%d%n2i*Lzd%Glr%d%n3i,orbs%nspinor))
      do rpqrs=1,nrpqrs
      read(111,*) irr(0), irr(1), irr(2)
      do ihpqrs=1,nhpqrs
      read(112,*) ip,iq,ir,is ; hpqrs_r=0._gp ;hpqrs_i=0._gp;
      ipt = ip  ;  irt = ir
      iqt = iq  ;  ist = is
      if(mspin .eq. 1) then
        ip = int(ip/2) ; ir = int(ir/2)  
        iq = int(iq/2) ; is = int(is/2)
      else
        ip = re_index2(ip,orbtot,orbocc)
        iq = re_index2(iq,orbtot,orbocc)
        ir = re_index2(ir,orbtot,orbocc)
        is = re_index2(is,orbtot,orbocc)
      end if

      tho1 = 0._gp ; tho1(:,1) = psirr(:,ip,1)*psirr(:,is,1) + psirr(:,ip,2)*psirr(:,is,2) ; tho1(:,2) = psirr(:,ip,1)*psirr(:,is,2) - psirr(:,ip,2)*psirr(:,is,1) 
      tho2 = 0._gp ; tho2(:,2) = psirr(:,iq,1)*psirr(:,ir,1) + psirr(:,iq,2)*psirr(:,ir,2) ; tho1(:,2) = psirr(:,iq,1)*psirr(:,ir,2) - psirr(:,iq,2)*psirr(:,ir,1)
      call tho_diet1(tho1(1,:),nn,cri,cspin, nn_new1) ; call tho_diet1(tho1(:,2),nn,cri,cspin, nn_new1)
      call tho_diet1(tho2(1,:),nn,cri,cspin, nn_new2) ; call tho_diet1(tho2(:,2),nn,cri,cspin, nn_new2) ; nn_new = nn_new2 
      if(nn_new1 .gt. nn_new2) nn_new = nn_new1
      allocate(tho_new1(nn_new,orbs%nspinor),thoi1(nn_new,orbs%nspinor)) ; tho_new1=0.d0 ; thoi1=0
      allocate(tho_new2(nn_new,orbs%nspinor),thoi2(nn_new,orbs%nspinor)) ; tho_new2=0.d0 ; thoi2=0
      
      call tho_diet2(tho1(:,1),nn,cri,nn_new,cspin, tho_new1(:,1),thoi1(:,1)) ; call tho_diet2(tho1(:,2),nn,cri,nn_new,cspin, tho_new1(:,2),thoi1(:,2))

      call tho_diet2(tho2(:,1),nn,cri,nn_new,cspin, tho_new2(:,1),thoi2(:,1)) ; call tho_diet2(tho2(:,2),nn,cri,nn_new,cspin, tho_new2(:,2),thoi2(:,2))

      !real(11)
      !$OMP PARALLEL DO PRIVATE(rx,ry,rz,rr,ii,jj,ecos,esin,eikr,hpqrs_t1,hpqrs_t2,hpqrs_t3) &
      !$OMP& REDUCTION(+:hpqrs_r,hpqrs_i)
      do i=1,nn_new,1 ; ii = thoi1(i,1)
        do j=1,nn_new,1 ; jj = thoi2(j,1)
          hpqrs_t1 = tho_new1(i) * tho_new2(j)
          if(ii .eq. jj) then
            hpqrs_t2 = r_r
          else
            rx = nii(ii) - nii(jj)
            ry = njj(ii) - njj(jj)
            rz = nkk(ii) - nkk(jj)
            rr = sqrt(rx**2 + ry**2 + rz**2)
            !real & imaginary part separate
            ecos = cos(ik(0)*rx+ik(1)*ry+ik(2)*rz)
            esin = sin(ik(0)*rx+ik(1)*ry+ik(2)*rz)
            hpqrs_t2 = ecos/rr
            hpqrs_t3 = esin/rr
          end if
          hpqrs_r = hpqrs_r + hpqrs_t1 * hpqrs_t2
          hpqrs_i = hpqrs_i + hpqrs_t1 * hpqrs_t3
        end do
      end do
      !$OMP END PARALLEL DO
  
      !real(22)
      !$OMP PARALLEL DO PRIVATE(rx,ry,rz,rr,ii,jj,ecos,esin,eikr,hpqrs_t1,hpqrs_t2,hpqrs_t3) &
      !$OMP& REDUCTION(+:hpqrs_r,hpqrs_i)
      do i=1,nn_new,1 ; ii = thoi1(i)
        do j=1,nn_new,1 ; jj = thoi2(j)
          hpqrs_t1 = tho_new1(i) * tho_new2(j)
          if(ii .eq. jj) then
            hpqrs_t2 = r_r
          else
            rx = nii(ii) - nii(jj)
            ry = njj(ii) - njj(jj)
            rz = nkk(ii) - nkk(jj)
            rr = sqrt(rx**2 + ry**2 + rz**2)
            !real & imaginary part separate
            ecos = cos(ik(0)*rx+ik(1)*ry+ik(2)*rz)
            esin = sin(ik(0)*rx+ik(1)*ry+ik(2)*rz)
            hpqrs_t2 = ecos/rr
            hpqrs_t3 = esin/rr
          end if
          hpqrs_r = hpqrs_r + hpqrs_t1 * hpqrs_t2
          hpqrs_i = hpqrs_i + hpqrs_t1 * hpqrs_t3
        end do
      end do
      !$OMP END PARALLEL DO

      write(555,35) ipt,iqt,irt,ist, hpqrs_r, hpqrs_i
      write( 6 ,36) ihpqrs,nhpqrs, ipt,iqt,irt,ist, hpqrs_r, hpqrs_i
      deallocate(tho_new1,thoi1,tho_new2,thoi2)
    end do
    end do
 35 format(2x,4i4,2f24.10)
 36 format(i3,2h /,i3,2x,4i4,2f24.10)
    close(112)
    deallocate(psirr,nii,njj,nkk, tho1,tho2)

!-------------------------------------------------------------------------------------------------------------------------------------

    open(6,file='hpqrs.log')
    !----------------------
    ! Free allocated space.
    !----------------------
    call deallocate_comms(comms)
    call deallocate_locreg_descriptors(Lzd%Glr)
    call deallocate_Lzd_except_Glr(Lzd)
    call deallocate_orbs(orbs)
    call deallocate_atoms_data(atoms)
    call xc_end(xc)
    call dpbox_free(dpcom)
    call pkernel_free(pkernel)
    call free_input_variables(inputs)
    
    !--------------------------------------
    !wait all processes before finalisation
    !--------------------------------------
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_FINALIZE(ierr)
    call f_lib_finalize()
    close(6) ; close(124)

stop
end program toy_model

!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
subroutine tho_diet1(tho,nn,cri,cspin, nn_new)
implicit none
  real*8 cri
  integer(kind=8) i,j
  integer  nn, nn_new, ncorb, cspin
  real*8,  dimension(nn) :: tho
  ncorb=0 
  do i=1,nn  
   if(real(cspin)*abs(tho(i)) .le. cri) cycle
   ncorb = ncorb + 1  
  end do
  nn_new = ncorb
! write(6,*) "number of new grids: ",nn_new
return
end
subroutine tho_diet2(tho,nn,cri,nn_new,cspin, tho_new,thoi)
implicit none
  real*8 cri
  integer(kind=8) i,j,k
  integer nn, nn_new, ncount, cspin
  integer, dimension(nn_new) :: thoi
  real*8,  dimension(nn_new) :: tho_new
  real*8,  dimension(nn    ) :: tho
  thoi=0 ; tho_new=0.d0 ; ncount=0

  do i=1,nn
    if(real(cspin)*abs(tho(i)) .le. cri) cycle
    ncount = ncount + 1
    tho_new(ncount) = tho(i)
       thoi(ncount) = i
  end do
return
end
!-------------------------------------------------------
! -- index in this code    0 1 2 3 4 5    6 7 8 9  -----
!                                 transfer to           
! --   chemical   index    0 2 4 1 3 5    6 8 7 9  -----
!-------------------------------------------------------
integer(kind=8) function re_index1(i,itot,iocc)
implicit none
  integer(kind=8) :: i
  integer         :: itot,iocc
  if(i .lt. iocc) then
    if(i .lt. iocc/2) then
      re_index1 =  i*2
    else
      re_index1 = (i-(iocc/2))*2 + 1
    end if
  else
    if(i .lt. (iocc+(itot-iocc)/2)) then
      re_index1 = iocc + (i-iocc)*2
    else
      re_index1 = 2*i - itot + 1
    end if
  end if
return
end
!-------------------------------------------------------
! --   chemical   index    0 2 4 1 3 5    6 8 7 9  -----
!                                 transfer to           
! -- index in this code    0 1 2 3 4 5    6 7 8 9  -----
!-------------------------------------------------------
integer(kind=8) function re_index2(i,itot,iocc)
implicit none
  integer(kind=8) :: i
  integer         :: itot,iocc
  if(i .lt. iocc) then
     if(mod(i,2) .eq. 0) then
        re_index2 = i/2
     else
        re_index2 = (i + iocc - 1)/2
     end if
  else
    if(mod(i,2) .eq. 0) then
      re_index2 = (i + iocc)/2 
    else
      re_index2 = (i + itot - 1)/2
    end if
  end if
return
end
