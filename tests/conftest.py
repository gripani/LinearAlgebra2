from typing import Callable

from numpy import sqrt, pi, cos

from pytest import fixture


@fixture
def a() -> list[list[float]]:
    return [
        [2., 9., 4.],
        [7., 5., 3.],
        [6., 1., 8.]
    ]


@fixture
def b() -> list[list[float]]:
    return [
        [1., sqrt(2.), 2.],
        [sqrt(2.), 3., sqrt(2.)],
        [2., sqrt(2.), 1.]
    ]


@fixture
def c() -> list[list[float]]:
    return [
        [2., 1., -1., 2.],
        [4., 5., -3., 6.],
        [-2., 5., -2., 6.],
        [4., 11., -4., 8.]
    ]


@fixture
def d() -> list[list[float]]:
    return [
        [ 5.,10.],
        [ 9.,18.],
        [ 4., 8.],
        [ 2., 4.]
    ]


@fixture
def fn() -> Callable[[float], float]:
    def _fn(x: float) -> float:
        return x * cos(1. * x**2) / (x**2 + 1.)
    return _fn


@fixture
def extremes() -> tuple[float, float]:
    return 0, pi
