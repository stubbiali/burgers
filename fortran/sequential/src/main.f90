program main
    use boundary, only: set_initial_conditions, apply_lateral_conditions
    use omp_lib, only: omp_get_wtime
    use space, only: advection, diffusion
    use time, only: forward_euler
    use utils, only: copy_fields, print_field
    implicit none

    ! control parameters
    real, parameter :: cfl = 0.1, nu = 0.1
    integer, parameter :: factor = 4, nt = 100
    integer, parameter :: nx = 10 * (2 ** factor) + 1, ny = 10 * (2 ** factor) + 1

    ! derived parameters
    real :: dt, dx, dy
    integer :: i, j, nb, t
    real :: u(nx, ny), u1(nx, ny), u2(nx, ny), u3(nx, ny)
    real :: v(nx, ny), v1(nx, ny), v2(nx, ny), v3(nx, ny)
    real :: start, finish
    dx = 1.0 / (nx - 1)
    dy = 1.0 / (ny - 1)
    dt = cfl * (dx ** 2)
    nb = 3

    ! set initial conditions
    call set_initial_conditions(u3, v3)
!    call print_field(u3)

    ! warm up cache
    call copy_fields(u3, v3, u, v)
    call forward_euler(nu, 0.0, dx, dy, nb, u, v, u, v, u1, v1)
    call apply_lateral_conditions(nb, u1, v1)
    call forward_euler(nu, 0.0, dx, dy, nb, u, v, u1, v1, u2, v2)
    call apply_lateral_conditions(nb, u2, v2)
    call forward_euler(nu, 0.0, dx, dy, nb, u, v, u2, v2, u3, v3)
    call apply_lateral_conditions(nb, u3, v3)

    call cpu_time(start)
    do t = 1, nt
        ! copy new fields into old fields
        call copy_fields(u3, v3, u, v)

        ! first step
        call forward_euler(nu, dt / 3.0, dx, dy, nb, u, v, u, v, u1, v1)
        call apply_lateral_conditions(nb, u1, v1)

        ! second step
        call forward_euler(nu, dt / 2.0, dx, dy, nb, u, v, u1, v1, u2, v2)
        call apply_lateral_conditions(nb, u2, v2)

        ! third step
        call forward_euler(nu, dt, dx, dy, nb, u, v, u2, v2, u3, v3)
        call apply_lateral_conditions(nb, u3, v3)

        ! print *, "Iteration ", t
    end do
    call cpu_time(finish)

    print *, "max(u) = ", maxval(u3), "min(u) = ", minval(u3)
    print *, "Run time: ",finish - start," seconds"
end program main