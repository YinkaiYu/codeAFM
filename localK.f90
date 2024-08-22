module LocalK_mod
    ! 该程序已经全面修改过了！！！
    use Multiply_mod
    implicit none
    
    public :: determinant_cmplx, inverse_cmplx
    private :: LocalK_metro, phi_new
    
    type(AccCounter) :: Acc_Kl, Acc_Kt
    real(kind=8), dimension(:,:,:), allocatable :: phi_new
    
contains
    subroutine LocalK_init()
        call Acc_Kl%init()
        call Acc_Kt%init()
        allocate(phi_new(Nboson, Lq, Ltrot))
        return
    end subroutine LocalK_init
    
    subroutine LocalK_clear()
        deallocate(phi_new)
        return
    end subroutine LocalK_clear
    
    subroutine LocalK_reset()
        call Acc_Kl%reset()
        call Acc_Kt%reset()
        phi_new = NsigL_K%phi
        return
    end subroutine LocalK_reset
    
    subroutine LocalK_metro(Gr, iseed, ii, ntau)
        ! 计算用于更新的ratio
        use MyMats
! Arguments:
	    complex(kind=8), dimension(Ndim, Ndim), intent(inout) :: Gr
        integer, intent(inout) :: iseed
        integer, intent(in) :: ii, ntau
!   Local: 
        real(kind=8), external :: ranf
        complex(kind=8) :: Proddet
        complex(kind=8), dimension(Norb*Nspin, Norb*Nspin) :: Prod, Prodinv, Gr_local, mat_tmp
        complex(kind=8) :: Vhlp(Norb*Nspin, Ndim), Uhlp(Ndim, Norb*Nspin), temp(Ndim, Norb*Nspin), Diff(Ndim, Ndim) ! 为Vhlp, Uhlp, temp 增加了2维
        real(kind=8) :: ratio_fermion, ratio_boson, ratio_re, ratio_re_abs
        real(kind=8) :: random, Xdif, xflip
        integer :: ns, P(Norb*Nspin), j, no, sign, nl, nr, nn
        real(kind=8), dimension(Nboson) :: vec_new, vec_old

! Local update on a two-component spin vector on space-time (ii, ntau)
        do ns = 1, Nboson
            xflip = ranf(iseed)
            Xdif = dble((xflip - 0.5) * abs(valr0))
            phi_new(ns, ii, ntau) = NsigL_K%phi(ns, ii, ntau) + Xdif
        enddo
        vec_new(:) = phi_new(:, ii, ntau)
        vec_old(:) = NsigL_K%phi(:, ii, ntau)
        ! 计算P数组，用于计算Gr_local
        nn = 0
        do ns = 1, Nspin
            do no = 1, Norb
                nn = nn + 1
                P(nn) = Latt%inv_dim_list(ii, no, ns)
            enddo
        enddo
        if (Latt%b_list(ii, 2) == 1) sign = 1
        if (Latt%b_list(ii, 2) == 2) sign = -1
! Calculate fermionic Metropolis ratio within 4*4 matrix space
        call Op_K%get_delta(vec_old, vec_new, sign) ! update Delta matrix in Op_K
        Prod = dcmplx(0.d0, 0.d0)
        do nr = 1, Norb * Nspin
            do nl = 1, Norb * Nspin
                Gr_local(nl, nr) = ZKRON(nl, nr) - Gr(P(nl), P(nr))
            enddo
        enddo
        call mmult(mat_tmp, Op_K%Delta, Gr_local) ! 4*4 matrix multiplication
        do nr = 1, Norb * Nspin
            do nl = 1, Norb * Nspin
                Prod(nl, nr) = ZKRON(nl, nr) + mat_tmp(nl, nr)
            enddo
        enddo
        ! 手动计算行列式
        Proddet = determinant_cmplx(Prod)
        ! 费米子部分的ratio，已经进行了平方
        ratio_fermion = real(Proddet * dconjg(Proddet))
! Calculate total Metropolis ratio
        ratio_boson = NsigL_K%bosonratio(phi_new, ii, ntau, Latt)
        ratio_re = dble(ratio_fermion * ratio_boson)
        ratio_re_abs = abs(ratio_re)
        random = ranf(iseed)
        write(6,*) 'ntau', ntau
        write(6,*) 'site', ii
        write(6,*) 'fermion ratio when update yukawa: ', Proddet
        write(6,*) 'boson ratio when update yukawa: ', ratio_boson
        write(6,*) 'update V ratio: ', ratio_re_abs
! Upgrade Green's function
        ! 接受更新
        if (ratio_re_abs .gt. random) then
            call Acc_Kl%count(.true.)
            Prodinv = inverse_cmplx(Prod) ! 计算逆矩阵
            Uhlp = dcmplx(0.d0, 0.d0); Vhlp = dcmplx(0.d0, 0.d0)
            temp = dcmplx(0.d0, 0.d0); Diff = dcmplx(0.d0, 0.d0)
! Vhlp(1:2, 1:Ndim) = Del(1:2) * (1 - Grup)(P(1):P(2), 1:Ndim); Uhlp(1:Ndim, 1:2) = Grup(1:Ndim, P(1):P(2))
            do no = 1, Norb * Nspin
                do j = 1, Ndim
                    Uhlp(j, no) = Gr(j, P(no))
                    Vhlp(no, j) = - Op_K%Delta(no, 1) * Gr(P(1), j) - Op_K%Delta(no, 2) * Gr(P(2), j) &
                        & - Op_K%Delta(no, 3) * Gr(P(3), j) - Op_K%Delta(no, 4) * Gr(P(4), j)
                enddo
                Vhlp(no, P(1)) = Vhlp(no, P(1)) + Op_K%Delta(no, 1)
                Vhlp(no, P(2)) = Vhlp(no, P(2)) + Op_K%Delta(no, 2)
                Vhlp(no, P(3)) = Vhlp(no, P(3)) + Op_K%Delta(no, 3)
                Vhlp(no, P(4)) = Vhlp(no, P(4)) + Op_K%Delta(no, 4)
            enddo
            call mmult(temp, Uhlp, Prodinv)
            call mmult(Diff, temp, Vhlp)
            Gr = Gr - Diff ! output Gr in each spin-orbital sector
! Flip: 
            NsigL_K%phi(:, ii, ntau) = phi_new(:, ii, ntau)
        else
            call Acc_Kl%count(.false.)
            phi_new(:, ii, ntau) = NsigL_K%phi(:, ii, ntau)
        endif
        return
    end subroutine LocalK_metro
    
    subroutine LocalK_prop_L(Prop, iseed, nt)
        class(Propagator), intent(inout) :: Prop
        integer, intent(inout) :: iseed
        integer, intent(in) :: nt
        integer :: ii
        do ii = Lq, 1, -1
            call LocalK_metro(Prop%Gr, iseed, ii, nt)
        enddo
        call Op_K%mmult_L(Prop%Gr, Latt, NsigL_K%phi, nt, 1)
        call Op_K%mmult_R(Prop%Gr, Latt, NsigL_K%phi, nt, -1)
        call Op_K%mmult_L(Prop%UUL, Latt, NsigL_K%phi, nt, 1)
        return
    end subroutine LocalK_prop_L
    
    subroutine LocalK_prop_R(Prop, iseed, nt)
        class(Propagator), intent(inout) :: Prop
        integer, intent(inout) :: iseed
        integer, intent(in) :: nt
        integer :: ii
        call Op_K%mmult_R(Prop%Gr, Latt, NsigL_K%phi, nt, 1)
        call Op_K%mmult_L(Prop%Gr, Latt, NsigL_K%phi, nt, -1)
        do ii = 1, Lq
            call LocalK_metro(Prop%Gr, iseed, ii, nt)
        enddo
        call Op_K%mmult_R(Prop%UUR, Latt, NsigL_K%phi, nt, 1)
        return
    end subroutine LocalK_prop_R
    
    subroutine LocalK_therm(ii, ntau, iseed)
! Arguments: 
        integer, intent(inout) :: iseed
        integer, intent(in) :: ii, ntau
! Local: 
        real(kind=8), external :: ranf
        real(kind=8) :: xflip, random, Xdif, ratio_boson
        integer :: ns

        do ns = 1, Nboson
            xflip = ranf(iseed)
            Xdif = dble((xflip - 0.5) * abs(valrt(ns)))
            phi_new(ns, ii, ntau) = NsigL_K%phi(ns, ii, ntau) + Xdif
        enddo
        ratio_boson = NsigL_K%bosonratio(phi_new, ii, ntau, Latt)
        random = ranf(iseed)
        if (ratio_boson .gt. random) then
            call Acc_Kt%count(.true.)
            NsigL_K%phi(:, ii, ntau) = phi_new(:, ii, ntau)
        else
            call Acc_Kt%count(.false.)
            phi_new(:, ii, ntau) = NsigL_K%phi(:, ii, ntau)
        endif
        return
    end subroutine LocalK_therm

    ! 计算复矩阵的行列式
    function determinant_cmplx(A) result(det)
        use MyMats
        complex(kind=8), dimension(:,:), intent(in) :: A
        complex(kind=8) :: det
        complex(kind=8), dimension(size(A,1), size(A,2)) :: LU
        integer, dimension(size(A,1)) :: ipiv
        integer :: info, i
        integer :: n

        n = size(A, 1)
        LU = A

        ! 使用zgetrf进行LU分解
        call zgetrf(n, n, LU, n, ipiv, info)
        
        ! 初始化行列式为1
        det = (1.0d0_8, 0.0d0_8)
        
        ! 计算行列式，乘以对角线元素
        do i = 1, n
            det = det * LU(i,i)
            if (ipiv(i) /= i) det = -det
        end do
    end function determinant_cmplx

    ! 计算复矩阵的逆
    function inverse_cmplx(A) result(invA)
        use MyMats
        complex(kind=8), dimension(:,:), intent(in) :: A
        complex(kind=8), dimension(size(A,1), size(A,2)) :: invA
        integer, dimension(size(A,1)) :: ipiv
        integer :: info
        integer :: n

        n = size(A, 1)
        invA = A

        ! 使用zgetrf进行LU分解
        call zgetrf(n, n, invA, n, ipiv, info)

        ! 使用zgetri计算矩阵的逆
        call zgetri(n, invA, n, ipiv, info)
    end function inverse_cmplx

end module LocalK_mod