import abc
from typing import Optional, TYPE_CHECKING, Type

from sympl._core.factory import AbstractFactory

from config import np

if TYPE_CHECKING:
    from .advection import Advection
    from .boundary import Boundary
    from .diffusion import Diffusion


class Stepper(AbstractFactory):
    def __init__(
        self,
        advection: "Advection",
        diffusion: "Diffusion",
        boundary: "Boundary",
        *,
        dtype: Type[np.float_]
    ):
        self.advection = advection
        self.diffusion = diffusion
        self.boundary = boundary
        self.dtype = dtype

    @abc.abstractmethod
    def __call__(
        self,
        dt: float,
        dx: float,
        dy: float,
        u_now: np.ndarray,
        v_now: np.ndarray,
        u_tmp: Optional[np.ndarray] = None,
        v_tmp: Optional[np.ndarray] = None,
        *,
        out_u_new: Optional[np.ndarray] = None,
        out_v_new: Optional[np.ndarray] = None
    ):
        pass


class ForwardEuler(Stepper):
    name = "forward_euler"

    def __init__(self, advection, diffusion, boundary, *, dtype=np.float64):
        super().__init__(advection, diffusion, boundary, dtype=dtype)
        ni, nj = self.boundary.ni, self.boundary.nj
        self.adv_u = np.empty((ni, nj), dtype=self.dtype)
        self.adv_v = np.empty((ni, nj), dtype=self.dtype)
        self.diff_u = np.empty((ni, nj), dtype=self.dtype)
        self.diff_v = np.empty((ni, nj), dtype=self.dtype)

    def __call__(
        self,
        dt,
        dx,
        dy,
        u_now,
        v_now,
        u_tmp=None,
        v_tmp=None,
        *,
        out_u_new=None,
        out_v_new=None
    ):
        u_tmp = u_tmp if u_tmp is not None else u_now
        v_tmp = v_tmp if v_tmp is not None else v_now
        u_new = out_u_new if out_u_new is not None else np.zeros_like(u_now)
        v_new = out_v_new if out_v_new is not None else np.zeros_like(v_now)

        ni, nj, nb = self.boundary.ni, self.boundary.nj, self.boundary.nb
        i, j = slice(nb, ni - nb), slice(nb, nj - nb)

        self.advection(
            dx,
            dy,
            u_tmp,
            v_tmp,
            nb=nb,
            out_adv_u=self.adv_u,
            out_adv_v=self.adv_v,
        )
        self.diffusion(
            dx,
            dy,
            u_tmp,
            v_tmp,
            nb=nb,
            out_diff_u=self.diff_u,
            out_diff_v=self.diff_v,
        )
        u_new[i, j] = u_now[i, j] - dt * (self.adv_u[i, j] - self.diff_u[i, j])
        v_new[i, j] = v_now[i, j] - dt * (self.adv_v[i, j] - self.diff_v[i, j])
        self.boundary(u_new, v_new)

        return u_new, v_new


class RK3WS(Stepper):
    name = "rk3ws"

    def __init__(self, advection, diffusion, boundary, *, dtype=np.float64):
        super().__init__(advection, diffusion, boundary, dtype=dtype)
        self.forward_euler = Stepper.factory(
            "forward_euler",
            self.advection,
            self.diffusion,
            self.boundary,
            dtype=self.dtype,
        )
        ni, nj = self.boundary.ni, self.boundary.nj
        self.u1 = np.empty((ni, nj), dtype=self.dtype)
        self.v1 = np.empty((ni, nj), dtype=self.dtype)
        self.u2 = np.empty((ni, nj), dtype=self.dtype)
        self.v2 = np.empty((ni, nj), dtype=self.dtype)

    def __call__(
        self,
        dt,
        dx,
        dy,
        u_now,
        v_now,
        u_tmp=None,
        v_tmp=None,
        *,
        out_u_new=None,
        out_v_new=None
    ):
        self.forward_euler(
            dt / 3, dx, dy, u_now, v_now, out_u_new=self.u1, out_v_new=self.v1
        )
        self.forward_euler(
            dt / 2,
            dx,
            dy,
            u_now,
            v_now,
            u_tmp=self.u1,
            v_tmp=self.v1,
            out_u_new=self.u2,
            out_v_new=self.v2,
        )
        return self.forward_euler(
            dt,
            dx,
            dy,
            u_now,
            v_now,
            u_tmp=self.u2,
            v_tmp=self.v2,
            out_u_new=out_u_new,
            out_v_new=out_v_new,
        )
