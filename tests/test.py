from time import perf_counter_ns
from operator import itemgetter
from typing import Callable

from numpy import isclose, transpose, allclose
from numpy.linalg import det, eigh, solve, inv
from scipy.integrate import quad

from app import linalg



class CountedFunction[T1, T2]:

    _fn: Callable
    _count: int

    def __init__(self, fn: Callable[[T1], T2]):
        self._fn = fn
        self._count = 0

    def __call__(self, x: T1) -> T2:
        self._count += 1
        return self._fn(x)
    
    @property
    def count(self) -> int:
        return self._count


def test_determinant(a: list[list[float]]):
    print()
    start = perf_counter_ns()
    d1 = linalg.get_determinant(a)
    end = perf_counter_ns()
    print(f"fortran wrapper determinant: {(end - start)/1000:.3f} μs")

    start = perf_counter_ns()
    d2 = det(a)
    end = perf_counter_ns()
    print(f"numpy linalg determinant: {(end - start)/1000:.3f} μs")
    
    assert isclose(d2, d1)


def test_eigen_solve(b: list[list[float]]):
    print()
    start = perf_counter_ns()
    eigenvalues, eigenvectors = linalg.eigen_solve(b)
    end = perf_counter_ns()
    print(f"fortran wrapper eigen_solve: {(end - start)/1000:.3f} μs")

    eigenvectors = list(map(list, zip(*eigenvectors)))
    idx, eigenvalues = zip(
        *sorted([
            (i, e) for i, e in enumerate(eigenvalues)
        ],
        key=itemgetter(1))
    )
    eigenvectors = [eigenvectors[i] for i in idx[::-1]]
    eigenvalues = eigenvalues[::-1]

    start = perf_counter_ns()
    numpy_eigenvalues, numpy_eigenvectors = eigh(b)
    end = perf_counter_ns()
    print(f"numpy linalg eigen_solve: {(end - start)/1000:.3f} μs")

    idx = numpy_eigenvalues.argsort()[::-1]
    numpy_eigenvalues, numpy_eigenvectors = numpy_eigenvalues[idx], numpy_eigenvectors[:, idx]
    numpy_eigenvectors = transpose(numpy_eigenvectors)

    for eigenvalue, numpy_eigenvalue in zip(eigenvalues, numpy_eigenvalues):
        assert isclose(eigenvalue, numpy_eigenvalue)
    for v, nv in zip(eigenvectors, numpy_eigenvectors):
        assert allclose(v, nv) or allclose(v, -nv)


def test_gauss_jordan_elimination(c: list[list[float]], d: list[list[float]]):
    print()
    start = perf_counter_ns()
    inv_c, sol = linalg.gauss_jordan_elimination(c, d)
    end = perf_counter_ns()
    print(f"fortran wrapper gauss_jordan: {(end - start)/1000:.3f} μs")

    start = perf_counter_ns()
    numpy_sol = solve(c, d)
    end = perf_counter_ns()
    numpy_inv_c = inv(c)
    print(f"numpy lingalg gauss:jordan: {(end - start)/1000:.3f} μs")

    assert allclose(sol, numpy_sol)
    assert allclose(inv_c, numpy_inv_c)


def test_integrate(fn: Callable[[float], float], extremes: tuple[float, float]):
    print()
    counted_fn = CountedFunction(fn)
    start = perf_counter_ns()
    ft_sol = linalg.integrate(counted_fn, *extremes)
    end = perf_counter_ns()
    print(f"fortran wrapper integrate: {(end - start)/1000:.3f} μs")
    print(f"callback calls: {counted_fn.count}")

    start = perf_counter_ns()
    np_sol, _, info = quad(fn, *extremes, full_output=1)
    end = perf_counter_ns()
    print(f"numpy linalg integrate: {(end - start)/1000:.3f} μs")
    print(f"callback calls: {info["neval"]}")

    assert isclose(ft_sol, np_sol)

