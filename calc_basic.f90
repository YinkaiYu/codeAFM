module CalcBasic ! Global parameters
    implicit none

! constants
    real(kind=8),           parameter           :: Zero = 1.0d-10
    real(kind=8),           parameter           :: PI = acos(-1.d0)
    real(kind=8),           parameter           :: upbound = 1.0d+200
! lattice parameters
    integer,                parameter           :: Norb = 2 ! orbital x or y in; band index
    integer,                parameter           :: Nspin = 2 ! fermion spin \up or \down
    integer,                parameter           :: Nsub = 2 ! sublattice; square lattice
    integer,                parameter           :: Nbond = 2 ! bonds per site; lattice direction
    integer,                parameter           :: Nboson = 3 ! O(3) direction; 3-component boson field 
    integer,                public              :: Nlx, Nly, NlxTherm, NlyTherm
    integer,                public              :: Lq, LqTherm
    integer,                public              :: Ndim, Nsite
    real(kind=8),           public              :: Dtau
    real(kind=8),           public              :: Beta
    integer,                public              :: Ltrot, LtrotTherm
    integer,                public              :: Lfam
! Hamiltonian parameters
    real(kind=8),           public,     save    :: RT(Norb, Nbond)
    real(kind=8),           public,     save    :: quartic = 0.d0 ! quartic coupling u
    real(kind=8),           public              :: mass ! quadratic mass term r
    real(kind=8),           public              :: lambda ! Yukawa coupling lambda
    real(kind=8),           public              :: hopping1 ! parallel hopping t_hx = t_vy
    real(kind=8),           public              :: hopping2 ! vertical hopping t~_vx = - t~_hy
    real(kind=8),           public              :: mu1 ! chemical potential mu_x
    real(kind=8),           public              :: mu2 ! chemical potential mu_y
    real(kind=8),           public,     save    :: velocity = 3.d0 ! boson velocitiy c
    integer,                public              :: Nflavor ! fermion flavor N_f
! update parameters
    real(kind=8),           public              :: valr0
    logical,                public              :: is_global ! global update switch: Wolff-shift joint update
    integer,                public              :: Nglobal
    real(kind=8),           public              :: valrs(Nbond)
! initial state parameters
    real(kind=8),           public              :: NB_field ! perpendicular flux quanta
    real(kind=8),           public              :: inconf ! Gaussian amplitude of initial bosonic fields
    integer,                public              :: absolute ! type of initial bosonic field configuration
    real(kind=8),           public              :: init(Nboson) ! initial bosonic field balance position
! process control parameters
    logical,                public              :: is_tau ! whether to calculate time-sliced Green function
    integer,                public              :: Nthermal ! calculate time-sliced Green function from the (Nthermal + 1)-th bin
    logical,                public              :: is_warm ! bosonic warm-up switch
    integer,                public              :: Nwarm
    real(kind=8),           public              :: valrt(Nboson)
    integer,                public              :: Nst
    integer,                public              :: Nwrap
    integer,                public              :: Nbin
    integer,                public              :: Nsweep
    integer,                public              :: ISIZE, IRANK, IERR
    
contains
    subroutine read_input()
        include 'mpif.h'
        if (IRANK == 0) then
            open(unit=20, file='paramC_sets.txt', status='unknown')
            read(20,*) mass, lambda, hopping1, hopping2, mu1, mu2, Nflavor
            read(20,*) Nlx, Nly, Ltrot, Beta
            read(20,*) NlxTherm, NlyTherm, LtrotTherm
            read(20,*) Nwrap, Nbin, Nsweep, valr0
            read(20,*) is_tau, Nthermal
            read(20,*) is_warm, Nwarm, valrt(1), valrt(2), valrt(3)
            read(20,*) is_global, Nglobal, valrs(1), valrs(2), valrs(3)
            read(20,*) inconf, absolute, init(1), init(2), init(3)
            read(20,*) NB_field
            close(20)
        endif 
!   MPI process: parallelization
        call MPI_BCAST(Beta, 1, MPI_Real8, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(mass, 1, MPI_Real8, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(lambda, 1, MPI_Real8, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(hopping1, 1, MPI_Real8, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(hopping2, 1, MPI_Real8, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(mu1, 1, MPI_Real8, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(mu2, 1, MPI_Real8, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(valr0, 1, MPI_Real8, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(valrs, Nboson, MPI_Real8, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(valrt, Nboson, MPI_Real8, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(inconf, 1, MPI_Real8, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(init, Nboson, MPI_Real8, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(NB_field, 1, MPI_Real8, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(Nflavor, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(absolute, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(Nlx, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(Nly, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(Ltrot, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(NlxTherm, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(NlyTherm, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(LtrotTherm, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(Nwrap, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(Nbin, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(Nwarm, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(Nglobal, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(Nthermal, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(Nsweep, 1, MPI_Integer, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(is_tau, 1, MPI_Logical, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(is_warm, 1, MPI_Logical, 0, MPI_COMM_WORLD, IERR)
        call MPI_BCAST(is_global, 1, MPI_Logical, 0, MPI_COMM_WORLD, IERR)
        return
    end subroutine read_input
    
    subroutine Params_set()
        RT(1,1) = hopping1 ! (x\up) orbital; lattice horizontal direction
        RT(1,2) = hopping2 ! (x\up) orbital; lattice vertical direction
        RT(2,1) = hopping2 ! (y\down) orbital; lattice horizontal direction
        RT(2,2) = hopping1 ! (y\down) orbital; lattice vertical direction
        Lq = Nlx * Nly
        LqTherm = NlxTherm * NlyTherm
        Dtau = Beta / dble(Ltrot)
        Ndim = Lq * Norb * Nspin
        Nsite = Lq * Norb
        Lfam = int(Lq / Nsub)
        if (mod(Ltrot, Nwrap) == 0) then
            Nst = Ltrot / Nwrap
        else
            write(6,*) "Ltrot is not a multiple of Nwrap"; stop
        endif
        return
    end subroutine Params_set
    
    integer function nranf(iseed, N)
        integer, intent(inout) :: iseed
        integer, intent(in) :: N
        real (kind=8), external :: ranf
        nranf  = nint(ranf(iseed)*dble(N) + 0.5)
        if (nranf .lt. 1 ) nranf = 1
        if (nranf .gt. N ) nranf = N 
        return
    end function nranf
    
    ! define the periodic boundary condition in x,y-directions.
    integer function npbc(nr, L)
        integer, intent(in) :: nr, L
        npbc = nr
        if (nr .gt. L) npbc = nr - L
        if (nr .lt. 1) npbc = nr + L
        return
    end function npbc
    
    pure function sqr_vec(vec)
        real(kind=8) :: sqr_vec
        real(kind=8), dimension(Nboson), intent(in) :: vec
        integer :: ns
        sqr_vec = 0.d0
        do ns = 1, Nboson
            sqr_vec = sqr_vec + vec(ns) * vec(ns)
        enddo
        return
    end function sqr_vec
    
    subroutine write_info()
        if (IRANK == 0) then
            open (unit=50, file='info.txt', status='unknown', action="write")
            write(50,*) '========================='
            write(50,*) 'SU(2) fermi-surface coupling to O(3) SDW boson'
            write(50,*) 'Linear lengh Lx                                :', Nlx
            write(50,*) 'Linear lengh Ly                                :', Nly
            write(50,*) 'Hopping x-orbital horizontal                   :', RT(1,1)
            write(50,*) 'Hopping x-orbital vertical                     :', RT(1,2)
            write(50,*) 'Hopping y-orbital horizontal                   :', RT(2,1)
            write(50,*) 'Hopping y-orbital vertical                     :', RT(2,2)
            write(50,*) 'Bosonic mass r                                 :', mass
            write(50,*) 'Yukawa coupling lambda=g/sqrt(N_f)             :', lambda
            write(50,*) 'Boson velocity c                               :', velocity
            write(50,*) 'Boson quartic term u                           :', quartic
            write(50,*) 'x-orbital chemical potential mu_x              :', mu1
            write(50,*) 'y-orbital chemical potential mu_y              :', mu2
            write(50,*) 'fermion flavor N_f                             :', Nflavor
            write(50,*) 'Local update Boson field magnitude Shift       :', valr0
            if (is_global) then
                write(50,*) '# Global                                   :', Nglobal ! global Metropolis algorithm
                write(50,*) 'Global shift Boson field magnitude x direction :', valrs(1)
                write(50,*) 'Global shift Boson field magnitude y direction :', valrs(2)
            endif
            if (is_warm) then
                write(50,*) '# Warm                                         :', Nwarm
                write(50,*) 'Thermalize Boson field magnitude x direction   :', valrt(1)
                write(50,*) 'Thermalize Boson field magnitude y direction   :', valrt(2)
            endif
            write(50,*) 'Magnitude of Gaussian distribution             :', inconf
            write(50,*) 'Choosing initial distribution                  :', absolute
            write(50,*) 'Theta                                          :', Beta
            write(50,*) 'Trotter number                                 :', Ltrot
            write(50,*) '=>Dtau                                         :', Dtau
            write(50,*) 'N_Ortho                                        :', Nwrap
            write(50,*) '# Bins                                         :', Nbin
            write(50,*) '# Local                                        :', Nsweep
            write(50,*) 'Length of fam.                                 :', Lfam
            write(50,*) 'Magnetic field (in units of Flux quanta )      :', NB_Field
            write(50,*) '# Nodes                                        :', ISIZE
            if (mod(Ltrot, Nwrap) .NE. 0) then
                write(50,*) 'Ltrot is not a multiple of Nwrap'; stop
            endif
        endif
        return
    end subroutine write_info
    
end module CalcBasic