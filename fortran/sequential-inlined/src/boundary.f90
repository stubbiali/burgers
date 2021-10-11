module boundary
    implicit none

    private
    public set_initial_conditions, apply_lateral_conditions

contains

    subroutine set_initial_conditions(u, v)
        implicit none

        real, intent(out) :: u(:, :), v(:, :)

        integer :: bx, by, nx, ny

        nx = size(u, 1)
        bx = nx / 4
        ny = size(u, 2)
        by = ny / 4

        u(:, :) = 1.0
        u(bx + 1:nx - bx, by + 1:ny - by) = 4.0
        v(:, :) = 1.0
        v(bx + 1:nx - bx, by + 1:ny - by) = 4.0
    end subroutine set_initial_conditions

    subroutine apply_lateral_conditions(nb, u, v)
        implicit none

        integer, intent(in) :: nb
        real, intent(out) :: u(:, :), v(:, :)

        integer :: nx, ny
        nx = size(u, 1)
        ny = size(u, 2)

        u(1:nb, :) = 1.0
        u(nx - nb + 1:nx, :) = 1.0
        u(nb + 1:nx - nb, 1:nb) = 1.0
        u(nb + 1:nx - nb, ny - nb + 1:ny) = 1.0
        v(1:nb, :) = 1.0
        v(nx - nb + 1:nx, :) = 1.0
        v(nb + 1:nx - nb, 1:nb) = 1.0
        v(nb + 1:nx - nb, ny - nb + 1:ny) = 1.0
    end subroutine apply_lateral_conditions
end module boundary