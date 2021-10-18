import abc

from sympl._core.factory import AbstractFactory

from config import np


class Boundary(AbstractFactory):
    def __init__(self, nx: int, ny: int, nb: int):
        self.nx = nx
        self.ny = ny
        self.nb = nb

    @property
    @abc.abstractmethod
    def ni(self) -> int:
        pass

    @property
    @abc.abstractmethod
    def nj(self) -> int:
        pass

    @abc.abstractmethod
    def get_numerical_field(self, phi: np.ndarray) -> np.ndarray:
        pass

    @abc.abstractmethod
    def get_physical_field(self, phi: np.ndarray) -> np.ndarray:
        pass

    @abc.abstractmethod
    def __call__(self, u: np.ndarray, v: np.ndarray) -> None:
        pass


class Plateau(Boundary):
    name = "plateau"

    @property
    def ni(self):
        return self.nx

    @property
    def nj(self):
        return self.ny

    def get_numerical_field(self, phi):
        return phi

    def get_physical_field(self, phi):
        return phi

    def __call__(self, u, v):
        u[: self.nb, :] = 1.0
        u[-self.nb :, :] = 1.0
        u[self.nb : -self.nb, : self.nb] = 1.0
        u[self.nb : -self.nb, -self.nb :] = 1.0
        v[: self.nb, :] = 1.0
        v[-self.nb :, :] = 1.0
        v[self.nb : -self.nb, : self.nb] = 1.0
        v[self.nb : -self.nb, -self.nb :] = 1.0
