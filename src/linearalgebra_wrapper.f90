module linalg_wrapper
    use iso_fortran_env, only: real64
    use iso_c_binding
    use linalg_mod

    implicit none

    abstract interface 
        function c_callback(x) result(funval) bind(c)
            import :: c_double
            real(c_double), value :: x
            real(c_double) :: funval
        end function
    end interface

    procedure(c_callback), pointer, save :: current_callback => null()

    contains

    function fn_adapter(x) result(y)
        use iso_fortran_env, only: real64
        real(real64), intent(in) :: x
        real(real64) :: y
        y = current_callback(x)
    end function

    subroutine c_integrate(fn, x0, x1, epsabs, epsrel, res, status) &
            bind(c, name='c_integrate')
        type(c_funptr), value, intent(in) :: fn
        real(c_double), intent(in) :: x0, x1
        real(c_double), intent(in) :: epsabs, epsrel
        real(c_double), intent(out) :: res
        integer(c_int), intent(out) :: status

        real(real64) :: r

        call c_f_procpointer(fn, current_callback)

        call integr_method( &
            fn_adapter, &
            real(x0, real64), real(x1, real64), &
            real(epsabs, real64), real(epsrel, real64), &
            r, status &
        )

        res = real(r, c_double)
        nullify(current_callback)

    end subroutine

    subroutine c_determinant(n, mat, det) bind(c, &
            name="calc_determinant")
        integer(c_int), value, intent(in) :: n
        real(c_double), intent(in), dimension(n, n) :: mat
        real(c_double), intent(out) :: det
        det = determinant(mat)
    end subroutine

    subroutine c_diagonalization(n, a, w) bind(c, &
            name="c_diagonalization")
        integer(c_int), value, intent(in) :: n
        real(c_double), intent(inout), dimension(n, n) :: a
        real(c_double), intent(out), dimension(n) :: w
        
        real(real64), dimension(n) :: d, e

        call tred2(n, a, d, e)
        call tqli(n, d, e, a)
        w = d
    end subroutine

    subroutine c_gaussjordan(n, m, a, b) bind(c, &
            name="c_gaussjordan")
        integer(c_int), value, intent(in)              :: n, m
        real(c_double), intent(inout), dimension(n, n) :: a
        real(c_double), intent(inout), dimension(n, m) :: b

        call gaussj(n, m, a, b)
    end subroutine

    subroutine c_least_square_fit(n, x, y, m, q, r2, sem, seq) bind(c, &
            name='c_least_square_fit')
        integer(c_int), value, intent(in)        :: n
        real(c_double), intent(in), dimension(n) :: x, y
        real(c_double), intent(out)              :: m, q, r2, sem, seq
        call least_square_fit(n, x, y, m, q, r2, sem, seq)
    end subroutine
end module linalg_wrapper
