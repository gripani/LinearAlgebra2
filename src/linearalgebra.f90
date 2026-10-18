module linalg_mod
    use iso_fortran_env, only: real64
    private :: calc_cofactor, swap, outerand, outerprod
    public :: determinant, gaussj, tred2, tqli, least_square_fit
    
    abstract interface
        function integrand(r1) result(funval)
            import :: real64
            real(real64), intent(in) :: r1
            real(real64) :: funval
        end function
    end interface
    
    contains

    recursive function determinant (mat) result(laplace_det)
        real(real64), dimension(:, :) :: mat
        real(real64) :: det, laplace_det
        real(real64), allocatable, dimension(:, :) :: cf

        n = size(mat, 1)
        if (n == 1) then
            det = mat(1, 1)
        else
            det = 0
            do i=1, n
                allocate(cf(n-1, n-1))
                call calc_cofactor(mat, i, 1, cf)
                det = det + ((-1)**(i+1)) * mat(i, 1) * determinant(cf)
                deallocate(cf)
            end do
        end if

        laplace_det = det

    end function

    subroutine gaussj(n, m, a, b)
        integer, intent(in)                         :: n, m 
        real(real64), intent(inout), dimension(n,n) :: a 
        real(real64), intent(inout), dimension(n,m) :: b 

        integer, dimension(n)                       :: ipiv, indxr, indxc 
        logical, dimension(n)                       :: lpiv 
        real(real64), dimension(n)                  :: dumc
        integer, target, dimension(2)               :: irc
        integer, pointer                            :: irow, icol 
        real(real64) pivinv

        irow => irc(1)
        icol => irc(2)
        ipiv = 0 
        do i=1,n
            lpiv = (ipiv == 0)
            irc = maxloc(array=abs(a), mask=outerand(lpiv, lpiv))
            ipiv(icol) = ipiv(icol) + 1
            if(ipiv(icol) .gt. 1) then
                print '(a)', 'nrerror: gaussj: singular matrix (1)' 
            endif  
            if(irow .ne. icol) then
                call swap(a(irow,:), a(icol,:))
                call swap(b(irow,:), b(icol,:))
            endif 
            indxr(i) = irow 
            indxc(i) = icol 
            if(a(icol,icol) .eq. .0_real64) then 
                print '(a)', 'nrerror: gaussj: singular matrix (2)'
            endif 
            pivinv = 1.0_real64 / a(icol, icol)
            a(icol, icol) = 1.0_real64 
            a(icol, :) = a(icol, :) * pivinv 
            b(icol, :) = b(icol, :) * pivinv 
            dumc = a(:, icol)
            a(:, icol) = .0_real64 
            a(icol, icol) = pivinv 
            a(:icol-1, :) = a(:icol-1, :) - outerprod(dumc(:icol-1), a(icol,:))
            b(:icol-1, :) = b(:icol-1, :) - outerprod(dumc(:icol-1), b(icol,:))
            a(icol+1:, :) = a(icol+1:, :) - outerprod(dumc(icol+1:), a(icol,:))
            b(icol+1:, :) = b(icol+1:, :) - outerprod(dumc(icol+1:), b(icol,:))
        enddo   
        do l=n,1,-1
            call swap(a(:, indxr(l)), a(:, indxc(l))) 
        enddo 
    
    end subroutine

    subroutine tred2(n, a, d, e)
        integer, intent(in) :: n
        real(real64), dimension(:,:), intent(inout) :: a
        real(real64), dimension(:), intent(out) :: d, e
        
        real(real64), dimension(size(a,1)) :: gg
        real(real64) :: scale, h, f, g, hh
        logical :: yesvec = .true.

        do i=n,2,-1
            l=i-1
            h=0.0
            if (l > 1) then
                scale=sum(abs(a(i,1:l)))
                if (scale == 0.0) then
                    e(i)=a(i,l)
                else
                    a(i,1:l)=a(i,1:l)/scale
                    h=sum(a(i,1:l)**2)
                    f=a(i,l)
                    g=-sign(sqrt(h),f)
                    e(i)=scale*g
                    h=h-f*g
                    a(i,l)=f-g
                    if (yesvec) a(1:l,i)=a(i,1:l)/h
                    do j=1,l
                        e(j)=(                                  &   
                            dot_product(a(j,1:j),a(i,1:j)) +    & 
                            dot_product(a(j+1:l,j),a(i,j+1:l))  &
                        )/h
                    end do
                    f=dot_product(e(1:l),a(i,1:l))
                    hh=f/(h+h)
                    e(1:l)=e(1:l)-hh*a(i,1:l)
                    do j=1,l
                        a(j,1:j)=a(j,1:j)-a(i,j)*e(1:j)-e(j)*a(i,1:j)
                    end do
                end if
            else
                e(i)=a(i,l)
            end if
            d(i)=h
        end do
        if (yesvec) d(1)=0.0
        e(1)=0.0
        do i=1,n
            if (yesvec) then
                l=i-1
                if (d(i) /= 0.0) then
                    gg(1:l)=matmul(a(i,1:l),a(1:l,1:l))
                    a(1:l,1:l)=a(1:l,1:l)-outerprod(a(1:l,i),gg(1:l))
                end if
                d(i)=a(i,i)
                a(i,i)=1.0
                a(i,1:l)=0.0
                a(1:l,i)=0.0
            else
                d(i)=a(i,i)
            end if
        end do
    end subroutine

    subroutine tqli(n, d, e, z)
        integer, intent(in) :: n
        real(real64), dimension(:), intent(inout) :: d, e
        real(real64), dimension(:,:), optional, intent(inout) :: z
        
        real(real64), dimension(size(e)) :: ff
        real(real64) :: dd, r, g, s, c, f, b, p

        e(:)=eoshift(e(:),1)
        do l=1,n
            iter=0
            iterate: do
                do m=l,n-1
                    dd=abs(d(m))+abs(d(m+1))
                    if (abs(e(m))+dd == dd) exit
                end do
                if (m == l) exit iterate
                iter=iter+1
                g=(d(l+1)-d(l))/(2.0_real64*e(l))
                r=sqrt(1.0_real64 + g**2)
                g=d(m)-d(l)+e(l)/(g+sign(r,g))
                s=1.0
                c=1.0
                p=0.0
                do i=m-1,l,-1
                    f=s*e(i)
                    b=c*e(i)
                    r=sqrt(f**2 + g**2)
                    e(i+1)=r
                    if (r == 0.0) then
                        d(i+1)=d(i+1)-p
                        e(m)=0.0
                        cycle iterate
                    end if
                    s=f/r
                    c=g/r
                    g=d(i+1)-p
                    r=(d(i)-g)*s+2.0_real64*c*b
                    p=s*r
                    d(i+1)=g+p
                    g=c*r-b
                    if (present(z)) then
                        ff(1:n)=z(1:n,i+1)
                        z(1:n,i+1)=s*z(1:n,i)+c*ff(1:n)
                        z(1:n,i)=c*z(1:n,i)-s*ff(1:n)
                    end if
                end do
                d(l)=d(l)-p
                e(l)=g
                e(m)=0.0
            end do iterate
        end do
    end subroutine

    subroutine least_square_fit(n, x, y, m, q, r2, sem, seq)
        integer, intent(in) :: n
        real(real64), intent(in), dimension(n) :: x, y
        real(real64), intent(out) :: m, q, r2, sem, seq
        real(real64) :: xm, ym, s, ssxx, ssyy, ssxy
        xm = sum(x(:)) / dble(n)
        ym = sum(y(:)) / dble(n)
        ssxx = sum((x(:) - xm)**2)
        ssyy = sum((y(:) - ym)**2)
        ssxy = sum((x(:) - xm) * (y(:) - ym))
        m = ssxy / ssxx
        q = ym - m * xm
        r2 = ssxy**2 / (ssxx*ssyy)
        s = sqrt((ssyy - m * ssxy) / dble(n-2))
        sem = s / sqrt(ssxx)
        seq = s * sqrt(1.0_real64 / dble(n) + xm**2/ssxx)
    end subroutine

    subroutine calc_cofactor(mat, mi, mj, cf)
        real(real64), dimension(:,:), intent(in)  :: mat
        integer, intent(in) :: mi, mj
        real(real64), dimension(:,:), intent(out) :: cf
        
        integer, dimension(2) :: msize
        
        msize = shape(mat)
        n = msize(1)

        l = 0
        k = 1
        do i=1, n
            if (i /= mi) then
                l = 1
                do j=1, n
                    if (j /= mj) then
                        cf(k, l) = mat(i, j)
                        l = l + 1
                    end if
                end do
                k = k + 1
            end if
        end do
    end subroutine

    function outerprod(a, b)
        real(real64), dimension(:), intent(in)    :: a, b 
        real(real64), dimension(size(a), size(b)) :: outerprod
        outerprod = spread(a, dim=2, ncopies=size(b)) * spread(b, dim=1, ncopies=size(a)) 
    end function 

    function outerand(a, b)
        logical, dimension(:), intent(in)    :: a, b 
        logical, dimension(size(a), size(b)) :: outerand
        outerand = spread(a, dim=2, ncopies=size(b)) .and. spread(b, dim=1, ncopies=size(a)) 
    end function

    subroutine swap(a, b)
        real(real64), dimension(:), intent(inout) :: a, b 
        real(real64), dimension(size(a))          :: dum 
        dum = a 
        a = b 
        b = dum 
    end subroutine 

    subroutine dqelg(n, epstab, res, abserr, res3la, nres)
        integer, intent(inout) :: n, nres
        real(real64), intent(inout) :: epstab(52), res3la(3)
        real(real64), intent(out) :: res, abserr

        real(real64) :: delta1, delta2, delta3, epsinf
        real(real64) :: err1, err2, err3, e0, e1, e1abs
        real(real64) :: e2, e3, rnew, ss, tol1, tol2, tol3
        real(real64) :: error
        integer :: i, ib, ib2, ie, indx, k1, k2, k3, num
        integer :: newelm

        nres = nres + 1
        abserr = huge(1.0_real64)
        res = epstab(n)

        if (n >= 3) then
            epstab(n + 2) = epstab(n)
            newelm = (n - 1) / 2
            epstab(n) = huge(1.0_real64)
            num = n
            k1 = n

            do i = 1, newelm
                k2 = k1 - 1
                k3 = k1 - 2

                rnew = epstab(k1 + 2)
                e0 = epstab(k3)
                e1 = epstab(k2)
                e2 = rnew
                e1abs = abs(e1)

                delta2 = e2 - e1
                err2 = abs(delta2)
                tol2 = max(abs(e2), e1abs) * epsilon(1.0_real64)

                delta3 = e1 - e0
                err3 = abs(delta3)
                tol3 = max(e1abs, abs(e0)) * epsilon(1.0_real64)

                if (err2 > tol2 .or. err3 > tol3) then
                    e3 = epstab(k1)
                    epstab(k1) = e1

                    delta1 = e1 - e3
                    err1 = abs(delta1)
                    tol1 = max(e1abs, abs(e3)) * epsilon(1.0_real64)

                    if (err1 > tol1 .and. err2 > tol2 .and. err3 > tol3) then
                        ss = 1.0_real64 / delta1 + &
                             1.0_real64 / delta2 - 1.0_real64 / delta3

                        epsinf = abs(ss * e1)

                        if (epsinf > 0.1e-03_real64) then
                            rnew = e1 + 1.0_real64 / ss
                            epstab(k1) = rnew
                            k1 = k1 - 2

                            error = err2 + abs(rnew - e2) + err3

                            if (error <= abserr) then
                                abserr = error
                                res = rnew
                            end if

                            cycle
                        end if
                    end if

                    n = 2 * i - 1
                    exit
                else
                    res = rnew
                    abserr = err2 + err3
                    abserr = max(abserr, &
                                 5.0_real64 * epsilon(1.0_real64) * abs(res))
                    return
                end if
            end do

            if (n == 50) n = 49

            ib = 1
            if (mod(num, 2) == 0) ib = 2
            ie = newelm + 1

            do i = 1, ie
                ib2 = ib + 2
                epstab(ib) = epstab(ib2)
                ib = ib2
            end do

            if (num /= n) then
                indx = num - n + 1

                do i = 1, n
                    epstab(i) = epstab(indx)
                    indx = indx + 1
                end do
            end if

            if (nres >= 4) then
                abserr = abs(res - res3la(3)) + &
                         abs(res - res3la(2)) + &
                         abs(res - res3la(1))

                res3la(1) = res3la(2)
                res3la(2) = res3la(3)
                res3la(3) = res
            else
                res3la(nres) = res
                abserr = huge(1.0_real64)
            end if
        end if

        abserr = max(abserr, &
                     5.0_real64 * epsilon(1.0_real64) * abs(res))

    end subroutine dqelg

    subroutine gk21(fn, a, b, res, abserr)
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        procedure(integrand) :: fn
        real(real64), intent(in)  :: a, b
        real(real64), intent(out) :: res, abserr

        real(real64), parameter :: xgk(11) = [ &
            0.9956571630258081_real64, &
            0.9739065285171717_real64, &
            0.9301574913557082_real64, &
            0.8650633666889845_real64, &
            0.7808177265864169_real64, &
            0.6794095682990244_real64, &
            0.5627571346686047_real64, &
            0.4333953941292472_real64, &
            0.2943928627014602_real64, &
            0.1488743389816312_real64, &
            0.0_real64 ]

        real(real64), parameter :: wgk(11) = [ &
            0.01169463886737187_real64, &
            0.03255816230796473_real64, &
            0.05475589657435199_real64, &
            0.07503967481091995_real64, &
            0.09312545458369761_real64, &
            0.1093871588022976_real64, &
            0.1234919762620659_real64, &
            0.1347092173114733_real64, &
            0.1427759385770601_real64, &
            0.1477391049013385_real64, &
            0.1494455540029169_real64 ]

        real(real64), parameter :: wg(5) = [ &
            0.06667134430868814_real64, &
            0.1494513491505806_real64, &
            0.2190863625159820_real64, &
            0.2692667193099963_real64, &
            0.2955242247147529_real64 ]

        real(real64) :: center, half, fc, f1, f2, fsum
        real(real64) :: resg, resk, resabs, resasc, mean
        real(real64) :: fv1(10), fv2(10)
        integer :: j, ig

        center = 0.5_real64 * (a + b)
        half   = 0.5_real64 * (b - a)

        fc = fn(center)
        resg = 0.0_real64
        resk = wgk(11) * fc
        resabs = wgk(11) * abs(fc)

        ig = 0
        do j = 1, 10
            f1 = fn(center - half * xgk(j))
            f2 = fn(center + half * xgk(j))
            fv1(j) = f1
            fv2(j) = f2
            fsum = f1 + f2

            resk = resk + wgk(j) * fsum
            resabs = resabs + wgk(j) * (abs(f1) + abs(f2))

            if (mod(j, 2) == 0) then
                ig = ig + 1
                resg = resg + wg(ig) * fsum
            end if
        end do

        res = resk * half
        resabs = resabs * abs(half)

        mean = 0.5_real64 * resk
        resasc = wgk(11) * abs(fc - mean)

        do j = 1, 10
            resasc = resasc + wgk(j) * &
                (abs(fv1(j) - mean) + abs(fv2(j) - mean))
        end do
        resasc = resasc * abs(half)

        abserr = abs((resk - resg) * half)

        if (resasc /= 0.0_real64 .and. abserr /= 0.0_real64) then
            abserr = resasc * min(1.0_real64, &
                (200.0_real64 * abserr / resasc)**1.5_real64)
        end if

        if (resabs > tiny(1.0_real64) / &
            (50.0_real64 * epsilon(1.0_real64))) then
            abserr = max(abserr, &
                50.0_real64 * epsilon(1.0_real64) * resabs)
        end if

    end subroutine

    subroutine integr_method(fn, x0, x1, epsabs, epsrel, res, status)
        procedure(integrand) :: fn
        real(real64), intent(in) :: x0, x1, epsabs, epsrel
        real(real64), intent(out) :: res
        integer, intent(out) :: status

        integer, parameter :: max_intervals = 1000
        integer, parameter :: limexp = 50

        real(real64) :: a(max_intervals), b(max_intervals)
        real(real64) :: area(max_intervals), err(max_intervals)
        real(real64) :: epstab(limexp + 2), res3la(3)
        real(real64) :: rlist2(limexp + 2)
        real(real64) :: total_error, tol, global_area
        real(real64) :: mid, left_value, left_error
        real(real64) :: right_value, right_error
        real(real64) :: area12, error12, old_error
        real(real64) :: errmax, small, erlarg, erlast
        real(real64) :: extrap_res, extrap_err, best_res, best_err
        real(real64) :: correction
        real(real64) :: old_a, old_b
        integer :: count, worst, nres, nseq
        integer :: nrmax, ktmin
        integer :: iroff1, iroff2, iroff3
        logical :: extrap, noext

        res = 0.0_real64
        status = 0

        if (epsabs <= 0.0_real64 .and. &
            epsrel < max(50.0_real64 * epsilon(1.0_real64), &
                         0.5e-28_real64)) then
            status = 6
            return
        end if

        if (x0 == x1) return

        count = 1
        a(1) = min(x0, x1)
        b(1) = max(x0, x1)

        call gk21(fn, a(1), b(1), area(1), err(1))

        global_area = area(1)
        total_error = err(1)
        tol = max(epsabs, epsrel * abs(global_area))

        if (total_error <= tol) then
            res = global_area
            if (x1 < x0) res = -res
            return
        end if

        epstab = 0.0_real64
        res3la = 0.0_real64
        rlist2 = 0.0_real64

        rlist2(1) = global_area
        nseq = 1
        nres = 0
        nrmax = 1
        ktmin = 0

        extrap = .false.
        noext = .false.

        iroff1 = 0
        iroff2 = 0
        iroff3 = 0

        small = 0.375_real64 * abs(b(1) - a(1))
        erlarg = total_error
        erlast = total_error
        correction = 0.0_real64

        best_res = global_area
        best_err = huge(1.0_real64)

        do while (count < max_intervals)

            worst = maxloc(err(1:count), dim=1)
            errmax = err(worst)
            old_error = errmax
            erlast = old_error

            old_a = a(worst)
            old_b = b(worst)
            mid = .5_real64 * (old_a + old_b)

            call gk21(fn, old_a, mid, left_value, left_error)
            call gk21(fn, mid, old_b, right_value, right_error)

            area12 = left_value + right_value
            error12 = left_error + right_error

            total_error = max(0.0_real64, &
                              total_error + error12 - old_error)
            global_area = global_area + area12 - area(worst)

            tol = max(epsabs, epsrel * abs(global_area))

            if (abs(area(worst) - area12) <= &
                1.0e-5_real64 * abs(area12) .and. &
                error12 >= 0.99_real64 * old_error) then
                if (extrap) then
                    iroff2 = iroff2 + 1
                else
                    iroff1 = iroff1 + 1
                end if
            end if

            if (count >= 10 .and. error12 > old_error) then
                iroff3 = iroff3 + 1
            end if

            if (right_error > left_error) then
                a(worst) = mid
                area(worst) = right_value
                err(worst) = right_error

                count = count + 1
                a(count) = a(worst)
                b(count) = mid
                area(count) = left_value
                err(count) = left_error
            else
                b(worst) = mid
                area(worst) = left_value
                err(worst) = left_error

                count = count + 1
                a(count) = mid
                b(count) = b(worst)
                area(count) = right_value
                err(count) = right_error
            end if

            global_area = sum(area(1:count))
            total_error = sum(err(1:count))
            tol = max(epsabs, epsrel * abs(global_area))

            if (total_error <= tol) then
                res = global_area
                exit
            end if

            if (iroff1 + iroff2 >= 10 .or. iroff3 >= 20) then
                status = 2
                exit
            end if

            if (count >= max_intervals) then
                status = 1
                exit
            end if

            if (noext) cycle

            worst = maxloc(err(1:count), dim=1)

            if (.not. extrap) then
                if (abs(b(worst) - a(worst)) > small) cycle
                extrap = .true.
                nrmax = 2
            end if

            if (nseq < limexp) then
                nseq = nseq + 1
                rlist2(nseq) = global_area
                epstab(1:nseq) = rlist2(1:nseq)

                call dqelg(nseq, epstab, extrap_res, extrap_err, &
                           res3la, nres)

                ktmin = ktmin + 1

                if (ktmin > 5 .and. best_err < 1.0e-3_real64 * total_error) then
                    status = 5
                    exit
                end if

                if (extrap_err < best_err) then
                    ktmin = 0
                    best_err = extrap_err
                    best_res = extrap_res
                    correction = total_error
                    tol = max(epsabs, epsrel * abs(best_res))

                    if (best_err <= tol) then
                        res = best_res
                        exit
                    end if
                end if
            else
                noext = .true.
            end if

            extrap = .false.
            small = 0.5_real64 * small
            erlarg = total_error
        end do

        if (best_err < total_error) then
            res = best_res
        else
            res = global_area
        end if

        if (count >= max_intervals .and. total_error > tol) status = 1

        if (x1 < x0) res = -res

    end subroutine integr_method

end module linalg_mod
