from ctypes import CDLL, c_int, c_double, c_void_p, byref, CFUNCTYPE, POINTER, Array
from os import path
from typing import Protocol, Callable, Final

type _M[T] = list[list[T]]
type _V[T] = list[T]
type _FN[T1, T2] = Callable[[T1], T2]
type _LSFR = tuple[float, float, float, float, float]


class LinearAlgebraWrapper:

    _dll: CDLL

    def __init__(self):
        try:
            self._dll = CDLL(path.join("lib", "linearalgebra_wrapper.so"))
        except OSError as exc:
            print(f"unable to load shared object: {exc}")
            exit(1)
        self._dll.calc_determinant.argtypes = [
            c_int,
            c_void_p,
            POINTER(c_double),
        ]
        self._dll.calc_determinant.restype = None

        self._dll.c_diagonalization.argtypes = [
            c_int,
            c_void_p,
            c_void_p,
        ]
        self._dll.c_diagonalization.restype = None

        self._dll.c_gaussjordan.argtypes = [
            c_int,
            c_int,
            c_void_p,
            c_void_p,
        ]
        self._dll.c_gaussjordan.restype = None

        self._dll.c_least_square_fit.argtypes = [
            c_int,
            POINTER(c_double),
            POINTER(c_double),
            POINTER(c_double),
            POINTER(c_double),
            POINTER(c_double),
        ]
        self._dll.c_least_square_fit.restype = None

        self._dll.c_integrate.argtypes = [
            CFUNCTYPE(c_double, c_double),
            POINTER(c_double),
            POINTER(c_double),
            POINTER(c_double),
            POINTER(c_double),
            POINTER(c_double),
            POINTER(c_int),
        ]
        self._dll.c_integrate.restype = None

    @staticmethod
    def __c_rect_matrix(a: _M[float]) -> Array[Array[c_double]]:
        dim2 = len(a)
        dim1 = len(a[0])
        return (c_double * dim2 * dim1) (
                *((c_double * dim2) (*col) for col in zip(*a))
        )

    @staticmethod
    def __c_square_matrix(a: _M[float]) -> Array[Array[c_double]]:
        dim = len(a[0])
        return (c_double * dim * dim) (
            *((c_double * dim) (*col) for col in zip(*a))
        )

    @staticmethod
    def __c_vector(a: _V[float]) -> Array[c_double]:
        return (c_double * len(a)) (*a)

    def get_determinant(self, a: _M[float]) -> float:
        n = len(a[0])
        c = c_double()
        mat = self.__c_square_matrix(a)
        self._dll.calc_determinant(c_int(n), byref(mat), byref(c))
        return c.value

    def eigen_solve(self, a: _M[float]) -> tuple[_V[float], _M[float]]:
        n = len(a[0])
        w = [.0 for _ in range(n)]
        mat_a = self.__c_square_matrix(a)
        vec_w = self.__c_vector(w)
        self._dll.c_diagonalization(c_int(n), byref(mat_a), byref(vec_w))
        eigenvalues = [vec_w[i] for i in range(n)]
        eigenvectors = [[mat_a[j][i] for j in range(n)] for i in range(n)]
        return eigenvalues, eigenvectors

    def gauss_jordan_elimination(self, a: _M[float], b: _M[float]) -> tuple[_M[float], _M[float]]:
        assert len(a[0]) == len(a[:][0]), 'a must be a square matrix'
        n = len(a[:][0])
        m = len(b[0][:])
        mat_a = self.__c_square_matrix(a)
        mat_b = self.__c_rect_matrix(b)
        self._dll.c_gaussjordan(
            c_int(n),
            c_int(m), 
            byref(mat_a), 
            byref(mat_b)
        )
        inv_a = [[mat_a[j][i] for j in range(n)] for i in range(n)]
        sol = [[mat_b[j][i] for j in range(m)] for i in range(n)]
        return inv_a, sol

    def least_square_fit(self, x: _V[float], y: _V[float]) -> _LSFR:
        n1 = len(x)
        n2 = len(y)
        if n1 != n2:
            raise ValueError("x and y must have the same length")
        n = c_int(n1)
        m = c_double()
        q = c_double()
        erm = c_double()
        erq = c_double()
        r2 = c_double()
        self._dll.c_least_square_fit(
            n, 
            byref(m), 
            byref(q), 
            byref(erm), 
            byref(erq), 
            byref(r2)
        )
        return m.value, q.value, erm.value, erq.value, r2.value

    def integrate(self, fn: _FN[float, float], a: float, b: float, 
        epsabs: float = 1.49e-8, epsrel: float = 1.49e-8) -> float:
        xi = c_double(a)
        xf = c_double(b)
        abs_tol = c_double(epsabs)
        rel_tol = c_double(epsrel)
        status = c_int()
        fnp = CFUNCTYPE(c_double, c_double)(fn)
        res = c_double()
        self._dll.c_integrate(
            fnp, 
            byref(xi), 
            byref(xf), 
            byref(abs_tol), 
            byref(rel_tol), 
            byref(res), 
            byref(status)
        )

        if status.value != 0:
            raise RuntimeError(f"Integration falied: (status={status.value})")
        return res.value
    
