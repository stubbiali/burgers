program test_boundary
    use boundary, only: apply_lateral_conditions, set_initial_conditions
    use utils, only: print_field
    implicit none

    integer, parameter :: nx = 11, ny = 11, nb = 3
    real :: u(nx, ny), v(ny, ny)

    u(:, :) = 0.0
    v(:, :) = 0.0
    call apply_lateral_conditions(nb, u, v)

    call print_field(u)
    print *, ""
    call print_field(v)
end program test_boundary