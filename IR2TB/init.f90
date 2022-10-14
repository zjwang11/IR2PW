module init

    use comms 
    implicit none 

contains

subroutine setarray()

    integer :: irot 

    allocate(EE(num_wann));                               EE = 0.d0 
    allocate(coeffa(num_sumorb, num_wann));               coeffa=czero 
    allocate(coeffb(num_sumorb, num_wann));               coeffb=czero 

    allocate(det_input(num_sym))
    allocate(angle_input(num_sym))
    allocate(axis_input(3,num_sym))
    allocate(tau_input(3,num_sym))

    do irot = 1, num_sym 
        det_input(irot) = det_read(irot)
        angle_input(irot) = angle_read(irot)
        axis_input(:,irot) = axis_read(:,irot)
        tau_input(:,irot) = tau_read(:,irot)
    enddo 

end subroutine setarray 


subroutine downarray()
    if (allocated(EE))        deallocate(EE)
    if (allocated(coeffa))    deallocate(coeffa)
    if (allocated(coeffb))    deallocate(coeffb)
    if (allocated(HmnR))      deallocate(HmnR)
end subroutine downarray 


subroutine read_HmnR()

    integer :: i, j, ir 
    integer :: itmp, jtmp 

    real(dp) :: r1tmp, r2tmp

    integer, parameter :: hr = 1000

    open(unit=hr, file=hrfile, status='old')
    read(hr, *)
    read(hr, *) num_wann
    read(hr, *) num_block
    num_sumorb = num_wann/nspin
    if (bot_band == 0) bot_band = 1
    if (top_band == 0) top_band = num_wann

    allocate(deg_block(num_block));                deg_block = 0
    allocate(vec_block(3, num_block));             vec_block = 0
    allocate(HmnR(num_wann, num_wann, num_block)); HmnR = czero

    read(hr, "(15I5)") deg_block(:)
    do ir = 1, num_block
        do j = 1, num_wann
        do i = 1, num_wann
            read(hr, *) vec_block(:,ir), itmp, jtmp, r1tmp, r2tmp
            HmnR(i,j,ir) = cmplx(r1tmp, r2tmp, dp)
        enddo 
        enddo

        if (deg_block(ir) == 1) cycle
        HmnR(:,:,ir) = HmnR(:,:,ir)/deg_block(ir)
        deg_block(ir) = 1
    enddo 

    close(hr)

end subroutine read_HmnR


subroutine read_tbbox()

    character(len=15)  :: chtp15

    integer            :: i, j, i1, i2, irot 
    integer            :: ierr 
    
    integer            :: kmesh, Nk 
    real(dp), allocatable, save :: node_kpoints(:,:)
    
    integer            :: itmp
    real(dp)           :: rtmp(8), rotmt(3,3)
    real(dp)           :: tmp(3), tmpr(3)
    complex(dp)        :: crotmt(3,3), protmt(3,3), drotmt(5,5), srotmt(2,2), frotmt(7,7)
    complex(dp)        :: cmat3(3,3), cmat5(5,5), cmat7(7,7)

    integer, parameter :: tbbox = 1001

    open(unit=tbbox, file='tbbox.in',form='formatted',status='old')

    chtp15 = 'case'
    call get_key_para_cht(chtp15, tbbox, casename)
    hrfile = trim(casename)//'_hr.dat'
    nspin = 1
    isSpinor = .false.
    if (casename == 'soc') then 
        nspin = 2
        isSpinor = .true.
    endif 

    chtp15 = 'proj'
    call get_key_para_loc(chtp15, tbbox)
    chtp15 = 'orbt'
    call get_key_para_intct(chtp15, tbbox, orbt)
    if (orbt == 1) then 
        write(6, *) "Orbital convention 1: &
                     s, px, py, pz, xy, yz, zx, x2-y2, 3z2-r2 !!!"
    elseif (orbt == 2) then 
        write(6, *) "Orbital convention 2: &
                     s, pz, px, py, 3z2-r2, xz, yz, x2-y2, xy(Wannier90) !!!"
    else
        write(*,*) "Error: orbt !!!!"
        stop 
    endif 
    chtp15 = 'ntau'
    call get_key_para_intct(chtp15, tbbox, ncenter) 
    allocate(type_center(ncenter), norb_center(ncenter), startorb_center(ncenter))
    allocate(pos_center(3, ncenter))
    startorb_center = 0
    do i = 1, ncenter 
        read(tbbox, *) pos_center(:, i), type_center(i), norb_center(i)
        if (i /= 1) startorb_center(i) = startorb_center(i-1) + norb_center(i-1)
    enddo 

    isInv         = .false.
    isSymmorphic  = .false. 
    chtp15 = 'unit_cell'
    call get_key_para_loc(chtp15, tbbox)
    do i = 1, 3
        read(tbbox, *) br2(:,i), (br4(i,j), j=1,3)
    enddo 
    DO irot = 1, MAXSYM
        read(tbbox, *, iostat=ierr) itmp, rtmp 
        if (ierr /= 0) exit 
        num_sym = irot 
        det_read(irot) = rtmp(1)
        angle_read(irot) = -rtmp(2)
        axis_read(:,irot) = rtmp(3:5)
        tau_read(:, irot) = rtmp(6:8)
        if (abs(rtmp(1)+1.d0+rtmp(2)) < 1.d-5) isInv = .true.
        if (abs(rtmp(6))+abs(rtmp(7))+abs(rtmp(8)) < 1.d-5) isSymmorphic = .true.
    ENDDO  ! irot

    ! kpath 
    chtp15 = 'kpoint'
    call get_key_para_loc(chtp15, tbbox)
    chtp15 = 'kmesh'
    call get_key_para_intct(chtp15, tbbox, kmesh)
    chtp15 = 'Nk'
    call get_key_para_intct(chtp15, tbbox, Nk)
    allocate(node_kpoints(3,Nk)); node_kpoints = 0.d0
    call get_key_para_nvec(tbbox, Nk, node_kpoints)

    allocate(len_k((Nk-1)*kmesh+1)); len_k = 0.d0 
    allocate(k(3,(Nk-1)*kmesh+1)); k = 0.d0 
    do i = 1, Nk - 1
        tmp = (node_kpoints(:,i+1)-node_kpoints(:,i))/kmesh
        do j = 1, kmesh
            k(:,(i-1)*kmesh+j) = node_kpoints(:,i)+tmp*(j-1)
            tmpr = matmul(tmp(:), br4(:,:))
            if (i/=1 .or. j/=1) then 
                len_k((i-1)*kmesh+j) = len_k((i-1)*kmesh+j-1) + &
                                  dsqrt(tmpr(1)**2+tmpr(2)**2+tmpr(3)**2)
            endif 
        enddo 
    enddo 
    k(:,(Nk-1)*kmesh+1) = node_kpoints(:,Nk)
    len_k = 2.d0*PI*len_k 
    num_kpoints = (Nk-1)*kmesh + 1

    close(tbbox)    

    write(6, 529) trim(casename)//'wann_TB' 

    isComplexWF = .true. 
    if (.not.isSpinor .and. isInv) isComplexWF = .false. 
    if (     isSymmorphic) write(6,'(A19)',advance='NO') ' Symmorphic crystal'
    if (.not.isSymmorphic) write(6,'(A23)',advance='NO') ' Non-symmorphic crystal'

    if (     isInv) write(6,'(A24)') ' with inversion symmetry'
    if (.not.isInv) write(6,'(A27)') ' without inversion symmetry'

    if (     isComplexWF) write(6, '(A23)') ' Complex eigenfunctions'
    if (.not.isComplexWF) write(6, '(A20)') ' Real eigenfunctions'

    if (     isSpinor) write(6, '(A45)') ' Spin-orbit eigenfunctions (->time inversion)'
    if (.not.isSpinor) write(6, '(A29)') ' No spin-orbit eigenfunctions'

      WRITE(6,590)
      WRITE(6,592)
      WRITE(6,595) ((br2(I1,I2),I2=1,3),I1=1,3)
      WRITE(6,591)
      WRITE(6,596)  (br4( 1,I2),I2=1,3),' : g1/2pi'
      WRITE(6,596)  (br4( 2,I2),I2=1,3),' : g2/2pi'
      WRITE(6,596)  (br4( 3,I2),I2=1,3),' : g3/2pi'
   
      RETURN
!529  FORMAT(1X,A10,/,1X,A4,' lattice')
 529  FORMAT(1X,A10,/)
 590  FORMAT(//,' Transformations:',/, &
             ' Direct lattice vectors in Cartesian coord. system (BR2)')
 591  FORMAT(' Reciprocal lattice vectors in Cartesian coord. system (BR4)')
 592  FORMAT('        t1 :            t2 :            t3 : ')
 595  FORMAT(3(3F16.8,/))
 596  FORMAT((3F16.8),A9)
end subroutine read_tbbox

end module init
