program ir2tb

    use comms
    use init
    use wave_data
    implicit none 

    ! command argument
    integer            :: narg, iarg, lens, stat
    character(len=100) :: arg, cmd 

    integer            :: i, kkk

    integer,     allocatable :: rot(:,:,:)
    real(dp),    allocatable :: SO3(:,:,:)
    complex(dp), allocatable :: SU2(:,:,:)
    complex(dp), allocatable :: kphase(:)
    complex(dp), allocatable :: tilte_mat(:,:,:)
    complex(dp), allocatable :: wfuncs(:,:,:)
    real(dp),    allocatable :: eig(:,:)
    integer :: ik,ib,nk,nb,n3atm

interface
subroutine irrep_bcs(sgn, num_sym, &
                      rot_input, tau_input, SO3_input, SU2_input, &
                      KKK, WK, kphase, &
                      num_bands, m, n, ene_bands, &
                      spinor, dim_basis, num_basis, &
                      coeffa, coeffb, &
                      G_phase_pw, rot_vec_pw, rot_mat_tb)

    integer, parameter :: dp = 8
    integer,     intent(in) :: sgn
    integer,     intent(in) :: num_sym 
    integer,     intent(in) :: rot_input(3,3,num_sym)
    real(dp),    intent(in) :: tau_input(3,num_sym)
    real(dp),    intent(in) :: SO3_input(3,3,num_sym)
    complex(dp), intent(in) :: SU2_input(2,2,num_sym)

    integer,     intent(in) :: KKK 
    real(dp),    intent(in) :: WK(3)
    complex(dp), intent(in) :: kphase(num_sym)

    integer,     intent(in) :: num_bands, m, n
    real(dp),    intent(in) :: ene_bands(num_bands) 

    logical,     intent(in) :: spinor
    integer,     intent(in) :: dim_basis, num_basis  
    complex(dp), intent(in) :: coeffa(dim_basis, num_bands)
    complex(dp), intent(in) :: coeffb(dim_basis, num_bands)
    complex(dp), intent(in), optional :: G_phase_pw(dim_basis, num_sym)
    complex(dp), intent(in), optional :: rot_vec_pw(dim_basis, num_sym)
    complex(dp), intent(in), optional :: rot_mat_tb(dim_basis, dim_basis, num_sym)
end subroutine irrep_bcs
end interface 

    sgn = 0
    bot_band = 0
    top_band = 0
    call get_command(cmd)
    write(*,*) 'Current command : ', trim(cmd) 
    narg = command_argument_count()
    write(*,*) 'Argument count : ', narg 

    if (narg == 0) then
        write(*,"('Please input the correct space group number (sgn) &
                   by the command below:')")
        write(*,"('###$: irvsp -sg $sgn')")
        stop
    else
        iarg = 1
        do while (.true.)
            call get_command_argument(iarg, arg, lens, stat)
            if (len_trim(arg) == 0) exit 
            if (trim(arg) == '-sg') then 
                iarg = iarg + 1
                call get_command_argument(iarg, arg, lens, stat)
                read(arg, *) sgn 
            elseif (trim(arg) == '-nb') then 
                iarg = iarg + 1
                call get_command_argument(iarg, arg, lens, stat)
                read(arg, *) bot_band
                iarg = iarg + 1
                call get_command_argument(iarg, arg, lens, stat)
                read(arg, *) top_band
            else
                write(*,"('Please input the correct space group number (sgn) &
                           by the command below : ')")
                write(*, "('###$: irvsp -sg $sgn')")
                stop
            endif 

            iarg = iarg + 1
        enddo 
    endif 

    if (sgn == 0 .or. sgn > 230) then
        write(*,"('Please input the correct space group number (sgn) &
                   by the command below:')")
        write(*,"('###$: irvsp -sg $sgn')")
        stop
    endif 

    call read_tbbox()

    call read_HmnR()

    call setarray()
    !
    open(55,file='modes.dat',status='old')
   !read(55,*) wzjnk,wzjnb!; print*,"nk,nb=",nk,nb
   !allocate (wfuncs(wzjnb,wzjnb,wzjnk)); wfuncs=0._8
    read(55,*) nk, nb, n3atm !; print*,"nk,nb=",nk,nb
   !if(nk_ .neq. nk) STOP ""
   !if(nb_ .neq. nb) STOP ""
    allocate (wfuncs(n3atm,nb,nk)); wfuncs=0._8
    allocate (eig(nb,nk)); eig=0._8
    do ik=1,nk
    do ib=1,nb
    read(55,"(F13.6)") eig(ib,ik)
    read(55,"(10F12.6)") wfuncs(:,ib,ik)
    !write(*,"(10F12.6)") wfuncs(:,ib,ik)
    enddo
    enddo
    close(55)
    !
    do kkk = 1, num_kpoints

        allocate(rot(3,3,num_sym))                         ; rot = 0
        allocate(SO3(3,3,num_sym))                         ; SO3 = 0.d0
        allocate(SU2(2,2,num_sym))                         ; SU2 = 0.d0
        allocate(kphase(num_sym))                          ; kphase = 0.d0
        allocate(tilte_mat(num_sumorb,num_sumorb,num_sym)) ; tilte_mat = 0.d0

        WK(:) = k(:,kkk)
        do i = 1, 3
            if (abs(3.d0*WK(i)-1.d0).lt.0.0002d0) WK(i) = 1.d0/3.d0
        enddo 

        call get_WF()
        
        !tmp
        coeffb(:,:)   =0._dp
        ! coeffa(:,:)   =0._dp
        !write(*,"(2f12.6)") coeffa(:,1)
        !print*,"wzj"
        !write(*,"(2f12.6)") coeffa(1,:)
        !write(*,"(2f12.6)") wfuncs(:,1,kkk)
        !print*,"wzj"
        !write(*,"(2f12.6)") wfuncs(:,2,kkk)
        coeffa(:,:)   =wfuncs(:,:,kkk)
        EE(:)         =eig(:,kkk)
        !tmp

        call tb_setup(WK, br2, &
                      num_sym, det_input, angle_input, axis_input, tau_input, &
                      ncenter, pos_center, &
                      num_sumorb, num_sumorb, norb_center, orbt, &
                      rot, SO3, SU2, &
                      kphase, tilte_mat)

        call irrep_bcs(sgn, num_sym, &
                        rot, tau_input, SO3, SU2, &
                        kkk, WK, kphase, &
                        num_wann, bot_band, top_band, EE, &
                        isSpinor, num_sumorb, num_sumorb, &
                        coeffa, coeffb, &
                        rot_mat_tb=tilte_mat)

        deallocate(rot)
        deallocate(SO3)
        deallocate(SU2)
        deallocate(kphase)
        deallocate(tilte_mat)

    enddo 
    deallocate (wfuncs)
    call downarray()

end program ir2tb 
