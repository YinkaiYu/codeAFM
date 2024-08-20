module LocalK_mod
    ! 该程序已经全面修改过了！！！
    use Multiply_mod
    implicit none
    
    public
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
        complex(kind=8), dimension(Norb, Norb) :: Prod, Prodinv, Gr_local, mat_tmp
        complex(kind=8) :: Vhlp(4, Ndim), Uhlp(Ndim, 4), temp(Ndim, 4), Diff(Ndim, Ndim) ! 为Vhlp, Uhlp, temp 增加了2维
        real(kind=8) :: ratio_fermion, ratio_boson, ratio_re, ratio_re_abs
        real(kind=8) :: random, Xdif, xflip
        integer :: ns, P(Norb), j, no, sign, nl, nr, nn
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
! Calculate fermionic Metropolis ratio within 2*2 matrix space
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
        Proddet = determinant_4x4(Prod)
        ! Proddet = Prod(1,1) * Prod(2,2) - Prod(1,2) * Prod(2,1)
        ! 费米子部分的ratio，已经进行了平方
        ratio_fermion = real(Proddet * dconjg(Proddet))
! Calculate total Metropolis ratio
        ratio_boson = NsigL_K%bosonratio(phi_new, ii, ntau, Latt)
        ratio_re = dble(ratio_fermion * ratio_boson)
        ! 防止符号问题
        ratio_re_abs = abs(ratio_re)
        random = ranf(iseed)
! Upgrade Green's function
        ! 接受更新
        if (ratio_re_abs .gt. random) then
            call Acc_Kl%count(.true.)
            ! Prodinv(1,1) = Prod(2,2)
            ! Prodinv(2,2) = Prod(1,1)
            ! Prodinv(1,2) = - Prod(1,2)
            ! Prodinv(2,1) = - Prod(2,1)
            Prodinv = adjugate_4x4(Prod) ! 计算伴随矩阵
            Prodinv = Prodinv / Proddet
            Uhlp = dcmplx(0.d0, 0.d0); Vhlp = dcmplx(0.d0, 0.d0)
            temp = dcmplx(0.d0, 0.d0); Diff = dcmplx(0.d0, 0.d0)
! Vhlp(1:2, 1:Ndim) = Del(1:2) * (1 - Grup)(P(1):P(2), 1:Ndim); Uhlp(1:Ndim, 1:2) = Grup(1:Ndim, P(1):P(2))
            do no = 1, Norb * Nspin
                do j = 1, Ndim
                    Uhlp(j, no) = Gr(j, P(no))
                    Vhlp(no, j) = - Op_K%Delta(no, 1) * Gr(P(1), j) - Op_K%Delta(no, 2) * Gr(P(2), j) &
                        - Op_K%Delta(no, 3) * Gr(P(3), j) - Op_K%Delta(no, 4) * Gr(P(4), j)
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

    ! 计算4x4矩阵的行列式
    real(8) function determinant_4x4(A)
        real(8), dimension(4, 4), intent(in) :: A
        determinant_4x4 = A(1,1) * cofactor(A, 1, 1) - A(1,2) * cofactor(A, 1, 2) &
                    + A(1,3) * cofactor(A, 1, 3) - A(1,4) * cofactor(A, 1, 4)
    end function determinant_4x4

    ! 计算4x4矩阵的伴随矩阵
    real(8) function adjugate_4x4(A)
        real(8), dimension(4, 4), intent(in) :: A
        real(8), dimension(4, 4) :: adj
        integer :: i, j
        do i = 1, 4
            do j = 1, 4
                adj(j, i) = cofactor(A, i, j) ! 计算余子式并转置
            end do
        end do
        adjugate_4x4 = adj
    end function adjugate_4x4

    ! 计算给定位置 (i, j) 的余子式
    real(8) function cofactor(A, i, j)
        real(8), dimension(4, 4), intent(in) :: A
        integer, intent(in) :: i, j
        real(8), dimension(3, 3) :: minor
        integer :: m, n, p, q
        m = 1
        do p = 1, 4
            if (p == i) cycle
            n = 1
            do q = 1, 4
                if (q == j) cycle
                minor(m, n) = A(p, q)
                n = n + 1
            end do
            m = m + 1
        end do
        cofactor = (-1.0d0)**(i + j) * determinant_3x3(minor)
    end function cofactor

    ! 计算3x3矩阵的行列式
    real(8) function determinant_3x3(minor)
        real(8), dimension(3, 3), intent(in) :: minor
        determinant_3x3 = minor(1,1) * (minor(2,2) * minor(3,3) - minor(2,3) * minor(3,2)) &
                - minor(1,2) * (minor(2,1) * minor(3,3) - minor(2,3) * minor(3,1)) &
                + minor(1,3) * (minor(2,1) * minor(3,2) - minor(2,2) * minor(3,1))
    end function determinant_3x3

end module LocalK_mod