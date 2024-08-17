module OperatorK_mod
    ! 该程序已经全面修改过了！！！
    use MyLattice
    implicit none
    public
    
    type :: OperatorPhonon
        real(kind=8), private :: alpha ! = Dtau*phg
        real(kind=8), private :: entryC
        complex(kind=8), private :: entryS
        real(kind=8), private :: entryA ! 由于O(3)对称性新加的A项
        
        complex(kind=8), public :: Delta(4,4) ! 用于更新Δ=exp(-alpha*V')*exp(alpha*V)-1(对角)
    contains
        procedure :: set => opK_set
        procedure, private :: get_exp => opK_get_exp
        procedure :: get_delta => opK_get_delta
        procedure :: mmult_R => opK_mmult_R
        procedure :: mmult_L => opK_mmult_L
    end type OperatorPhonon
    
contains
    subroutine opK_set(this)
        class(OperatorPhonon), intent(inout) :: this
        this%alpha = Dtau * lambda
        return
    end subroutine opK_set
    
    subroutine opK_get_exp(this, vec, nflag)
        ! 计算exp(D(NF))
        class(OperatorPhonon), intent(inout) :: this
        integer, intent(in) :: nflag ! +1 or -1; propagating direction
        real(kind=8), dimension(Nboson), intent(in) :: vec ! two-component order parameter with space-time coordinate (ii, ntau)
        real(kind=8) :: magnitude
        ! sqr_vec为求vec的模长的平方
        magnitude = sqrt(sqr_vec(vec)) ! 模长
        this%entryC = cosh( nflag * this%alpha * magnitude )
        this%entryS = sinh( nflag * this%alpha * magnitude ) * dcmplx(-vec(1), vec(2)) / magnitude
        this%entryA = sinh( nflag * this%alpha * magnitude ) * vec(3) / magnitude
        return
    end subroutine opK_get_exp
    
    subroutine opK_get_delta(this, vec_old, vec_new, sign)
        ! 计算用于更新的Δ=exp(-alpha*V')*exp(alpha*V)-1(对角)
        class(OperatorPhonon), intent(inout) :: this
        ! Nboson为玻色场自由度
        real(kind=8), dimension(Nboson), intent(in) :: vec_old, vec_new
        integer, intent(in) :: sign ! = \pm 1
        real(kind=8) :: C_new, C_old, A_new, A_old, mag_new, mag_old
        complex(kind=8) :: S_new, S_old
        ! 手动计算Δ=exp(-alpha*V')*exp(alpha*V)-1(对角)
        mag_new = sqrt(sqr_vec(vec_new))
        C_new = cosh( this%alpha * mag_new )
        S_new = sinh( this%alpha * mag_new ) * dcmplx(- vec_new(1), vec_new(2)) / mag_new
        A_new = sinh( this%alpha * mag_new ) * vec_new(3) / mag_new
        mag_old = sqrt(sqr_vec(vec_old))
        C_old = cosh( -this%alpha * mag_old ) ! 此处有没有正负一样，因为cosh()是偶函数
        S_old = sinh( -this%alpha * mag_old ) * dcmplx(- vec_old(1), vec_old(2)) / mag_old
        A_old = sinh( -this%alpha * mag_old ) * vec_old(3) / mag_old
        this%Delta(1,1) = C_new * C_old + S_new * dconjg(S_old) + A_new * A_old - 1.d0
        this%Delta(1,2) = (C_new * S_old + S_new * C_old ) * dble(sign)
        this%Delta(1,3) = S_new * A_old - A_new * S_old
        this%Delta(1,4) = (- C_new * A_old - A_new * C_old ) * dble(sign)
        this%Delta(2,1) = (dconjg(S_new) * C_old + C_new * dconjg(S_old) ) * dble(sign)
        this%Delta(2,2) = dconjg(S_new) * S_old + C_new * C_old + A_new * A_old - 1.d0
        this%Delta(2,3) = (A_new * C_old + C_new * A_old) * dble(sign)
        this%Delta(2,4) = - dconjg(S_new) * A_old + A_new * dconjg(S_old)
        this%Delta(3,1) = - dconjg(S_new) * A_old + A_new * dconjg(S_old)
        this%Delta(3,2) = (C_new * A_old + A_new * C_old) * dble(sign)
        this%Delta(3,3) = A_new * A_old + C_new * C_old + dconjg(S_new) * S_old - 1.d0
        this%Delta(3,4) = (C_new * dconjg(S_old) + dconjg(S_new) * C_old ) * dble(sign)
        this%Delta(4,1) = (- A_new * C_old - C_new * A_old ) * dble(sign)
        this%Delta(4,2) = - A_new * S_old + S_new * A_old
        this%Delta(4,3) = (S_new * C_old + C_new * S_old) * dble(sign)
        this%Delta(4,4) = A_new * A_old + S_new * dconjg(S_old) + C_new * C_old - 1.d0
        ! this%Delta(1,1) = C_new * C_old + S_new * dconjg(S_old) - 1.d0
        ! this%Delta(1,2) = (C_new * S_old + S_new * C_old) * dble(sign)
        ! this%Delta(2,1) = (C_new * dconjg(S_old) + dconjg(S_new) * C_old) * dble(sign)
        ! this%Delta(2,2) = C_new * C_old + dconjg(S_new) * S_old - 1.d0
        return
    end subroutine opK_get_delta
    
    subroutine opK_mmult_R(this, Mat, Latt, phi, ntau, nflag)
        ! 计算 Mat = EXP(D(NF)) * Mat
! In Mat Out U(NF) * EXP(D(NF)) * Mat
! Arguments: 
        class(OperatorPhonon), intent(inout) :: this
        complex(kind=8), dimension(Ndim, Ndim), intent(inout) :: Mat
        class(SquareLattice), intent(in) :: Latt
        real(kind=8), dimension(Nboson, Lq, Ltrot), intent(in) :: phi
        integer, intent(in) :: ntau, nflag
! Local: 
        integer :: P(Norb + 2), ii, no, j, sign ! 新加了2个P，用于遍历不同自旋的轨道
        real(kind=8), dimension(Nboson) :: vec
        complex(kind=8), dimension(2 + 2, Ndim) :: Vhlp ! 2 + 2 是因为新加了2个P，用于遍历不同自旋的轨道

        do ii = 1, Lq
            vec(:) = phi(:, ii, ntau)
            call this%get_exp(vec, nflag) ! output entryC and entryS 和 entryA
            if (Latt%b_list(ii, 2) == 1) sign = 1
            if (Latt%b_list(ii, 2) == 2) sign = -1
            do no = 1, Norb
                P(no) = Latt%inv_o_list(ii, no)
                P(no+2) = Latt%inv_o_list(ii, no) + 2 * Lq ! 新加的P，用于遍历不同自旋的轨道
            enddo
            Vhlp = dcmplx(0.d0, 0.d0)
            ! 计算MAT中不同行的结果（左乘了矩阵exp(-alpha*V)）
            do j = 1, Ndim
                Vhlp(1, j) = this%entryC * Mat(P(1), j) + sign * this%entryS * Mat(P(2), j) - sign * this%entryA * Mat(P(4), j)
                Vhlp(2, j) = sign * dconjg(this%entryS) * Mat(P(1), j) + this%entryC * Mat(P(2), j) + sign * this%entryA * Mat(P(3), j)
                Vhlp(3, j) = sign * this%entryA * Mat(P(2), j) + this%entryC * Mat(P(3), j) + sign * dconjg(this%entryS) * Mat(P(4), j)
                Vhlp(4, j) = - sign * this%entryA * Mat(P(1), j) + sign * this%entryS * Mat(P(3), j) + this%entryC * Mat(P(4), j)
            enddo
            do j = 1, Ndim
                Mat(P(1), j) = Vhlp(1, j)
                Mat(P(2), j) = Vhlp(2, j)
                Mat(P(3), j) = Vhlp(3, j)
                Mat(P(4), j) = Vhlp(4, j)
            enddo
        enddo
        return
    end subroutine opK_mmult_R
    
    subroutine opK_mmult_L(this, Mat, Latt, phi, ntau, nflag) 
        ! 计算 Mat = EXP(D(NF)) * Mat
! In Mat Out Mat* EXP(D(NF)) * UT(NF)
! Arguments:
        class(OperatorPhonon), intent(inout) :: this
        complex(kind=8), dimension(Ndim, Ndim), intent(inout) :: Mat
        class(SquareLattice), intent(in) :: Latt
        real(kind=8), dimension(Nboson, Lq, Ltrot), intent(in) :: phi
        integer, intent(in) :: ntau, nflag
! Local: 
        integer :: P(Norb + 2), j, ii, no, sign ! 新加了2个P，用于遍历不同自旋的轨道
        real(kind=8), dimension(Nboson) :: vec
        complex(kind=8), dimension(Ndim, 2 + 2) :: Uhlp ! 2 + 2 是因为新加了2个P，用于遍历不同自旋的轨道

        do ii = 1, Lq
            vec(:) = phi(:, ii, ntau)
            call this%get_exp(vec, nflag) ! output entryC and entryS 和 entryA
            if (Latt%b_list(ii, 2) == 1) sign = 1
            if (Latt%b_list(ii, 2) == 2) sign = -1
            do no = 1, Norb
                P(no) = Latt%inv_o_list(ii, no)
                P(no+2) = Latt%inv_o_list(ii, no) + 2 * Lq ! 新加的P，用于遍历不同自旋的轨道
            enddo
            Uhlp = dcmplx(0.d0, 0.d0)
            do j = 1, Ndim
                Uhlp(j, 1) = Mat(j, P(1)) * this%entryC + sign * Mat(j, P(2)) * dconjg(this%entryS) - sign * Mat(j, P(4)) * this%entryA
                Uhlp(j, 2) = sign * Mat(j, P(1)) * this%entryS + Mat(j, P(2)) * this%entryC + sign * Mat(j, P(3)) * this%entryA
                Uhlp(j, 3) = sign * Mat(j, P(2)) * this%entryA + Mat(j, P(3)) * this%entryC + sign * Mat(j, P(4)) * this%entryS
                Uhlp(j, 4) = - sign * Mat(j, P(1)) * this%entryA + sign * Mat(j, P(3)) * dconjg(this%entryS) + Mat(j, P(4)) * this%entryC
            enddo
            do j = 1, Ndim
                Mat(j, P(1)) = Uhlp(j, 1)
                Mat(j, P(2)) = Uhlp(j, 2)
                Mat(j, P(3)) = Uhlp(j, 3)
                Mat(j, P(4)) = Uhlp(j, 4)
            enddo
        enddo
        return
    end subroutine opK_mmult_L
end module OperatorK_mod