subroutine get_key_para_int(keyword, nfile, intn)

    character(len=15), intent(inout) :: keyword
    integer,           intent(in)    :: nfile
    integer,           intent(out)   :: intn

    character(len=15) :: keyword2 
    character(len=50) :: chartmp
    integer           :: length
    integer           :: ierr 
    integer           :: i, j

    rewind(nfile)
    do while (.true.)
        read(nfile, "(A)", iostat=ierr) chartmp 
        if (ierr /= 0) exit 
        length = len_trim(chartmp)
        do i = 1, length 
            if (chartmp(i:i) == '=' .or. chartmp(i:i) == ':') exit 
        enddo 
        if (i == length+1) cycle 
        read(chartmp(1:i-1), *) keyword2 
        if (trim(adjustl(keyword2)) == keyword) then 
            read(chartmp(i+1:length), *) intn
            exit 
        endif 
    enddo 

end subroutine get_key_para_int


subroutine get_key_para_cht(keyword, nfile, chat)
    
    character(len=15), intent(inout) :: keyword
    integer,           intent(in)    :: nfile 
    character(len=15), intent(out)   :: chat 

    character(len=15) :: keyword2 
    character(len=50) :: chartmp 

    integer           :: length 
    integer           :: ierr
    integer           :: i, j

    rewind(nfile)
    do while (.true.)
        read(nfile, "(A)", iostat=ierr) chartmp 
        if (ierr /= 0) exit 
        length = len_trim(chartmp)
        do i = 1, length 
            if (chartmp(i:i) == '=' .or. chartmp(i:i) == ':') exit
        enddo 
        if (i == length+1) cycle 
        read(chartmp(1:i-1), *) keyword2
        if (trim(adjustl(keyword2)) == keyword) then 
            read(chartmp(i+1:length), *) chat 
            exit
        endif 
    enddo 

end subroutine get_key_para_cht 


subroutine get_key_para_rel(keyword, nfile, relt)

    character(len=15), intent(inout) :: keyword
    integer,           intent(in)    :: nfile 
    real(8),           intent(out)   :: relt

    character(len=15) :: keyword2
    character(len=50) :: chartmp

    integer           :: length 
    integer           :: ierr 
    integer           :: i, j

    rewind(nfile)
    do while (.true.)
        read(nfile, "(A)", iostat=ierr) chartmp
        if (ierr /= 0) exit 
        length = len_trim(chartmp)
        do i = 1, length 
            if (chartmp(i:i) == '=' .or. chartmp(i:i) == ':') exit 
        enddo 
        if (i == length + 1) cycle 
        read(chartmp(1:i-1), *) keyword2
        if (trim(adjustl(keyword2)) == keyword) then 
            read(chartmp(i+1:length), *) relt
            exit
        endif 
    enddo 

end subroutine get_key_para_rel 


subroutine get_key_para_vec(keyword, nfile, nv, vect)

    character(len=15), intent(inout) :: keyword 
    integer,           intent(in)    :: nfile 
    integer,           intent(in)    :: nv 
    real(8),           intent(out)   :: vect(nv)

    character(len=15) :: keyword2
    character(len=50) :: chartmp

    integer           :: length
    integer           :: ierr
    integer           :: i, j

    rewind(nfile)
    do while (.true.)
        read(nfile, "(A)", iostat=ierr) chartmp
        if (ierr /= 0) exit
        length = len_trim(chartmp)
        do i = 1, length 
            if (chartmp(i:i) == '=' .or. chartmp == ':') exit 
        enddo 
        if (i == length + 1) cycle 
        read(chartmp(1:i-1), *) keyword2
        if (trim(adjustl(keyword2)) == keyword) then 
            read(chartmp(i+1:length), *) vect(:)
            exit
        endif 
    enddo 

end subroutine get_key_para_vec


subroutine get_key_para_loc(keyword, nfile)

    character(len=15), intent(inout) :: keyword
    integer,           intent(in)    :: nfile 

    character(len=15) :: keyword2 
    character(len=50) :: chartmp 
    
    integer           :: length 
    integer           :: ierr
    integer           :: i, j

    rewind(nfile)
    do while (.true.)
        read(nfile, "(A)", iostat=ierr) chartmp 
        length = len_trim(chartmp)
        do i = 1, length
            if (chartmp(i:i) == '=' .or. chartmp(i:i) == ':') exit  
        enddo 
        if (i == length + 1) cycle 
        read(chartmp(1:i-1), *) keyword2 
        if (trim(adjustl(keyword2)) == keyword) exit 
    enddo 

end subroutine get_key_para_loc 


subroutine get_key_para_intct(keyword, nfile, intn)

    character(len=15), intent(inout) :: keyword
    integer,           intent(in)    :: nfile 
    integer,           intent(out)   :: intn 

    character(len=15) :: keyword2
    character(len=50) :: chartmp 

    integer           :: length 
    integer           :: ierr
    integer           :: i, j

    do while (.true.)
        read(nfile, "(A)", iostat=ierr) chartmp 
        if (ierr /= 0) exit 
        length = len_trim(chartmp)
        do i = 1, length 
            if (chartmp(i:i) == '=' .or. chartmp(i:i) == ':') exit 
        enddo 
        if (i == length + 1) cycle 
        read(chartmp(1:i-1), *) keyword2 
        if (trim(adjustl(keyword2)) == keyword) then 
            read(chartmp(i+1:length), *) intn 
            exit 
        endif 
    enddo 

end subroutine get_key_para_intct


subroutine get_key_para_nvec(nfile, nk, kpoints)

    integer, intent(in)  :: nfile 
    integer, intent(in)  :: nk 
    real(8), intent(out) :: kpoints(3, 0:nk-1)

    integer :: ierr
    integer :: i, j

    do j = 0, nk-1
        read(nfile, *, iostat=ierr) kpoints(:,j)
        if (ierr /= 0) exit 
    enddo 

end subroutine get_key_para_nvec
