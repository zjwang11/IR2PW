module comms

    integer,           public, parameter :: dp = 8
    
    real(dp),          public, parameter :: PI = 3.141592653589793238462643383279d0 

    complex(dp),       public, parameter :: czero = (0._dp, 0._dp)
    complex(dp),       public, parameter :: cone  = (1._dp, 0._dp)
    complex(dp),       public, parameter :: CI    = (0._dp, 1._dp)

    integer,           public, parameter :: MAXSYM = 96

    integer,           public, parameter :: mix2l = 12

    integer,           public, save      :: sgn
    integer,           public, save      :: bot_band
    integer,           public, save      :: top_band

    character(len=15), public, save      :: casename 
    character(len=22), public, save      :: hrfile 

    logical,           public, save      :: isSymmorphic
    logical,           public, save      :: isComplexWF
    logical,           public, save      :: isSpinor
    logical,           public, save      :: isInv
    logical,           public, save      :: isSpinPola
    
    integer,           public, save      :: nspin
    integer,           public, save      :: nrot
    integer,           public ,save      :: ncenter
    integer,           public, save      :: orbt 

    integer,           public, save      :: num_kpoints
    integer,           public, save      :: num_sym
    integer,           public, save      :: num_wann
    integer,           public, save      :: num_sumorb
    integer,           public, save      :: num_block

    real(dp),          public, save      :: br2(3,3), br4(3,3)

    real(dp),              public, save  :: det_read(MAXSYM)
    real(dp),              public, save  :: angle_read(MAXSYM)
    real(dp),              public, save  :: axis_read(3,MAXSYM)
    real(dp),              public, save  :: tau_read(3,MAXSYM)
    real(dp), allocatable, public, save  :: det_input(:)    
    real(dp), allocatable, public, save  :: angle_input(:)
    real(dp), allocatable, public, save  :: axis_input(:,:)
    real(dp), allocatable, public, save  :: tau_input(:,:)

    real(dp),          public, save      :: WK(3)
    integer,           public, save      :: num_litt_group
    integer,           public, save      :: litt_group(MAXSYM)

    integer,  allocatable, public, save  :: type_center(:)
    integer,  allocatable, public, save  :: norb_center(:)
    integer,  allocatable, public, save  :: startorb_center(:)
    real(dp), allocatable, public, save  :: pos_center(:,:)
    real(dp), allocatable, public, save  :: rot_orb(:,:,:,:)

    real(dp), allocatable, public, save  :: len_k(:)
    real(dp), allocatable, public, save  :: k(:,:)

    integer,  allocatable, public, save  :: deg_block(:)
    integer,  allocatable, public, save  :: vec_block(:,:)
    complex(dp), allocatable, public, save :: HmnR(:,:,:)

    real(dp),    allocatable, public, save  :: EE(:)
    complex(dp), allocatable, public, save  :: coeffa(:,:)
    complex(dp), allocatable, public, save  :: coeffb(:,:)

end module comms
