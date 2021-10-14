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

        integer :: i, j, nx, ny

        nx = size(u, 1)
        ny = size(u, 2)
        !$acc parallel loop present(u, v)
        do j = 1, ny
            if ((j <= nb) .or. (j > ny - nb)) then
                !$acc loop
                do i = 1, nx
                    u(i, j) = 1.0
                    v(i, j) = 1.0
                end do
            else
                !$acc loop
                do i = 1, nb
                    u(i, j) = 1.0
                    v(i, j) = 1.0
                end do
                !$acc loop
                do i = nx - nb + 1, nx
                    u(i, j) = 1.0
                    v(i, j) = 1.0
                end do
            end if
        end do
    end subroutine apply_lateral_conditions
end module boundary