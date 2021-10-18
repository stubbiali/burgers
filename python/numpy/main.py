import click

from sympl._core.time import Timer

from config import np
from model import BurgersModel


@click.command()
@click.option("--nu", type=float, default=0.1, help="Diffusivity.")
@click.option(
    "--nx",
    type=int,
    default=161,
    help="Number of grid points along x-direction.",
)
@click.option(
    "--ny",
    type=int,
    default=161,
    help="Number of grid points along y-direction.",
)
@click.option("--cfl", type=float, default=0.1, help="CFL number.")
@click.option("--nt", type=int, default=100, help="Number of time iterates.")
def main(nu, nx, ny, cfl, nt):
    # initialize model
    burgers = BurgersModel(
        nu,
        nx,
        ny,
        advection_type="fifth_order_upwind",
        diffusion_type="fourth_order",
        stepper_type="rk3ws",
        state_factory_type="plateau",
        boundary_type="plateau",
        dtype=np.float64,
    )

    # set time-step
    dt = cfl * (min(burgers.dx, burgers.dy) ** 2)

    # warm-up cache and get initial state
    u, v = burgers.run(dt=0.0, nt=1)

    # run
    uf, vf = burgers.run(dt=dt, nt=nt, u=u, v=v)

    # log
    print(f"Validation: max(u) = {uf.max():.8f}, min(u) = {uf.min():.8f}")
    print(f"Run time: {Timer.get_time('time_integration', 's'):.5f} seconds")


if __name__ == "__main__":
    main()
