import abc
from typing import Tuple, Type

from sympl._core.factory import AbstractFactory

from config import np


class StateFactory(AbstractFactory):
    @abc.abstractmethod
    def __call__(
        self, nx: int, ny: int, *, dtype: Type[np.float_] = np.float64
    ) -> Tuple[np.ndarray, np.ndarray]:
        pass


class Plateau(StateFactory):
    name = "plateau"

    def __call__(self, nx, ny, *, dtype=np.float64):
        nx = max(nx, 1)
        bx = nx // 4
        ny = max(ny, 1)
        by = ny // 4

        u = np.ones((nx, ny), dtype=dtype)
        u[bx : nx - bx, by : ny - by] = 4.0
        v = np.ones((nx, ny), dtype=dtype)
        v[bx : nx - bx, by : ny - by] = 4.0

        return u, v
