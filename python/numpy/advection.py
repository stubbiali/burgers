import abc
from typing import Optional, Tuple

from sympl._core.factory import AbstractFactory

from config import np


class Advection(AbstractFactory):
    nb: int = None

    @abc.abstractmethod
    def __call__(
        self,
        dx: float,
        dy: float,
        u: np.ndarray,
        v: np.ndarray,
        *,
        nb: Optional[int] = None,
        out_adv_u: Optional[np.ndarray] = None,
        out_adv_v: Optional[np.ndarray] = None,
    ) -> Tuple[np.ndarray, np.ndarray]:
        pass


class FifthOrderUpwind(Advection):
    name = "fifth_order_upwind"
    nb = 3

    def __call__(
        self, dx, dy, u, v, *, nb=None, out_adv_u=None, out_adv_v=None
    ):
        nb = self.nb if nb is None else max(nb, self.nb)
        adv_u = out_adv_u if out_adv_u is not None else np.zeros_like(u)
        adv_v = out_adv_v if out_adv_v is not None else np.zeros_like(v)

        def advection_x(vx, phi):
            return (
                vx[nb:-nb, nb:-nb]
                / (60 * dx)
                * (
                    45
                    * (
                        phi[nb + 1 : -nb + 1, nb:-nb]
                        - phi[nb - 1 : -nb - 1, nb:-nb]
                    )
                    - 9
                    * (
                        phi[nb + 2 : -nb + 2, nb:-nb]
                        - phi[nb - 2 : -nb - 2, nb:-nb]
                    )
                    + (
                        phi[nb + 3 : (-nb + 3 if nb > 3 else None), nb:-nb]
                        - phi[nb - 3 : -nb - 3, nb:-nb]
                    )
                )
            ) - np.abs(vx[nb:-nb, nb:-nb]) / (60 * dx) * (
                (
                    phi[nb + 3 : (-nb + 3 if nb > 3 else None), nb:-nb]
                    + phi[nb - 3 : -nb - 3, nb:-nb]
                )
                - 6
                * (
                    phi[nb + 2 : -nb + 2, nb:-nb]
                    + phi[nb - 2 : -nb - 2, nb:-nb]
                )
                + 15
                * (
                    phi[nb + 1 : -nb + 1, nb:-nb]
                    + phi[nb - 1 : -nb - 1, nb:-nb]
                )
                - 20 * phi[nb:-nb, nb:-nb]
            )

        def advection_y(vy, phi):
            return (
                vy[nb:-nb, nb:-nb]
                / (60 * dy)
                * (
                    45
                    * (
                        phi[nb:-nb, nb + 1 : -nb + 1]
                        - phi[nb:-nb, nb - 1 : -nb - 1]
                    )
                    - 9
                    * (
                        phi[nb:-nb, nb + 2 : -nb + 2]
                        - phi[nb:-nb, nb - 2 : -nb - 2]
                    )
                    + (
                        phi[nb:-nb, nb + 3 : (-nb + 3 if nb > 3 else None)]
                        - phi[nb:-nb, nb - 3 : -nb - 3]
                    )
                )
            ) - np.abs(vy[nb:-nb, nb:-nb]) / (60 * dy) * (
                (
                    phi[nb:-nb, nb + 3 : (-nb + 3 if nb > 3 else None)]
                    + phi[nb:-nb, nb - 3 : -nb - 3]
                )
                - 6
                * (
                    phi[nb:-nb, nb + 2 : -nb + 2]
                    + phi[nb:-nb, nb - 2 : -nb - 2]
                )
                + 15
                * (
                    phi[nb:-nb, nb + 1 : -nb + 1]
                    + phi[nb:-nb, nb - 1 : -nb - 1]
                )
                - 20 * phi[nb:-nb, nb:-nb]
            )

        adv_u[nb:-nb, nb:-nb] = advection_x(u, u) + advection_y(v, u)
        adv_v[nb:-nb, nb:-nb] = advection_x(u, v) + advection_y(v, v)

        return adv_u, adv_v
