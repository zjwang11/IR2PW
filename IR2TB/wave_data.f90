module wave_data

    use comms
    implicit none

    complex(dp), allocatable, save :: Hk(:,:)
    complex(dp), allocatable, save :: evec(:,:)

contains

subroutine get_WF()

    integer :: i, info 

    if (.not.allocated(Hk))   allocate(  Hk(num_wann,num_wann))
    if (.not.allocated(evec)) allocate(evec(num_wann,num_wann))

    info = 0
    call get_Hk()
    evec = czero 
    call myzheev(num_wann, num_wann, Hk, EE, evec ,info)
    coeffa = czero; coeffb = czero
    coeffa(1:num_sumorb, 1:num_wann) = evec(1:num_sumorb, 1:num_wann)
    if (nspin == 2) coeffb(1:num_sumorb, 1:num_wann) = evec(num_sumorb+1:2*num_sumorb, 1:num_wann)

end subroutine get_WF 

subroutine get_Hk()

    integer  :: i
    real(dp) :: kdotr 

    Hk = czero 
    do i = 1, num_block
        kdotr = dot_product(WK(:), vec_block(:,i))
        Hk(1:num_wann,1:num_wann) = Hk(1:num_wann,1:num_wann) + &
                HmnR(1:num_wann, 1:num_wann, i)*exp(CI*kdotr*2.d0*PI)
    enddo 

end subroutine get_Hk 


subroutine myzheev(lda, N, Hmn, eval, evec, info)

    integer,     intent(in)  :: lda
    integer,     intent(in)  :: N
    complex(dp), intent(in)  :: Hmn(lda, N)
    real(dp),    intent(out) :: eval(N)
    complex(dp), intent(out) :: evec(lda, N)
    integer,     intent(out) :: info 

    ! dummy integer variables for lapack call
    integer     :: LWORK
    real(dp)    :: RWORK(3*N)
    COMPLEX(dp) :: WORK(2*N)
    complex(dp) :: A(N,N)

    LWORK = 2*N

    A      = czero
    A(:,:) = Hmn(1:N,1:N)
    call ZHEEV('V','U',N,A,N,eval,WORK,LWORK,RWORK,INFO)
    evec(1:N,1:N) = A(:,:)

    return 

end subroutine myzheev 

end module wave_data 
