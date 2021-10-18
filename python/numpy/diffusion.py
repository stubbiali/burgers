import abc
from typing import Optional, Tuple

from sympl._core.factory import AbstractFactory

from config import np


class Diffusion(AbstractFactory):
    nb: int = None

    def __init__(self, nu: float) -> None:
        self.nu = nu

    @abc.abstractmethod
    def __call__(
        self,
        dx: float,
        dy: float,
        u: np.ndarray,
        v: np.ndarray,
        *,
        nb: Optional[int] = None,
        out_diff_u: Optional[np.ndarray] = None,
        out_diff_v: Optional[np.ndarray] = None,
    ) -> Tuple[np.ndarray, np.ndarray]:
        pass


class FourthOrder(Diffusion):
    name = "fourth_order"
    nb = 2

    def __call__(
        self, dx, dy, u, v, *, nb=None, out_diff_u=None, out_diff_v=None
    ):
        nb = self.nb if nb is None else max(nb, self.nb)
        diff_u = out_diff_u if out_diff_u is not None else np.zeros_like(u)
        diff_v = out_diff_v if out_diff_v is not None else np.zeros_like(v)

        diff_u[nb:-nb, nb:-nb] = self.nu * (
            (
                -u[nb - 2 : -nb - 2, nb:-nb]
                + 16 * u[nb - 1 : -nb - 1, nb:-nb]
                - 30 * u[nb:-nb, nb:-nb]
                + 16 * u[nb + 1 : -nb + 1, nb:-nb]
                - u[nb + 2 : (-nb + 2 if nb > 2 else None), nb:-nb]
            )
            / (12 * (dx ** 2))
            + (
                -u[nb:-nb, nb - 2 : -nb - 2]
                + 16 * u[nb:-nb, nb - 1 : -nb - 1]
                - 30 * u[nb:-nb, nb:-nb]
                + 16 * u[nb:-nb, nb + 1 : -nb + 1]
                - u[nb:-nb, nb + 2 : (-nb + 2 if nb > 2 else None)]
            )
            / (12 * (dy ** 2))
        )
        diff_v[nb:-nb, nb:-nb] = self.nu * (
            (
                -v[nb - 2 : -nb - 2, nb:-nb]
                + 16 * v[nb - 1 : -nb - 1, nb:-nb]
                - 30 * v[nb:-nb, nb:-nb]
                + 16 * v[nb + 1 : -nb + 1, nb:-nb]
                - v[nb + 2 : (-nb + 2 if nb > 2 else None), nb:-nb]
            )
            / (12 * (dx ** 2))
            + (
                -v[nb:-nb, nb - 2 : -nb - 2]
                + 16 * v[nb:-nb, nb - 1 : -nb - 1]
                - 30 * v[nb:-nb, nb:-nb]
                + 16 * v[nb:-nb, nb + 1 : -nb + 1]
                - v[nb:-nb, nb + 2 : (-nb + 2 if nb > 2 else None)]
            )
            / (12 * (dy ** 2))
        )

        return diff_u, diff_v
