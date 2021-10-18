from typing import Optional, Tuple, Type

from sympl._core.time import Timer

from advection import Advection
from boundary import Boundary
from config import np
from diffusion import Diffusion
from state import StateFactory
from stepper import Stepper


class BurgersModel:
    def __init__(
        self,
        nu: float,
        nx: int,
        ny: int,
        advection_type: str,
        diffusion_type: str,
        stepper_type: str,
        state_factory_type: str,
        boundary_type: str,
        *,
        dtype: Type[np.float_] = np.float64
    ) -> None:
        self.nx = max(nx, 1)
        self.ny = max(ny, 1)
        self.dtype = dtype

        self.dx = 1 / (self.nx - 1)
        self.dy = 1 / (self.ny - 1)

        self.state_factory = StateFactory.factory(state_factory_type)

        self.advection = Advection.factory(advection_type)
        self.diffusion = Diffusion.factory(diffusion_type, nu)
        self.nb = max(self.advection.nb, self.diffusion.nb)
        self.boundary = Boundary.factory(boundary_type, nx, ny, self.nb)
        self.stepper = Stepper.factory(
            stepper_type,
            self.advection,
            self.diffusion,
            self.boundary,
            dtype=dtype,
        )

        self.u_new: Optional[np.ndarray] = None
        self.v_new: Optional[np.ndarray] = None

    def get_initial_state(self) -> Tuple[np.ndarray, np.ndarray]:
        u, v = self.state_factory(self.nx, self.ny, dtype=self.dtype)
        u = self.boundary.get_numerical_field(u)
        v = self.boundary.get_numerical_field(v)
        return u, v

    def run(
        self,
        dt: float,
        nt: int,
        u: Optional[np.ndarray] = None,
        v: Optional[np.ndarray] = None,
    ) -> Tuple[np.ndarray, np.ndarray]:
        if u is None or v is None:
            u, v = self.get_initial_state()

        u_now, v_now = u.copy(), v.copy()

        Timer.start("time_integration")
        for _ in range(nt):
            self.u_new, self.v_new = self.stepper(
                dt,
                self.dx,
                self.dy,
                u_now,
                v_now,
                out_u_new=self.u_new,
                out_v_new=self.v_new,
            )
            self.u_new, u_now = u_now, self.u_new
            self.v_new, v_now = v_now, self.v_new
        Timer.stop()

        return u_now, v_now
