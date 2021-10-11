module time
    use space, only: advection, diffusion
    implicit none

    private
    public forward_euler

contains

    subroutine forward_euler(nu, dt, dx, dy, nb, u_now, v_now, u_tmp, v_tmp, u_new, v_new)
        implicit none

        real, intent(in) :: nu, dt, dx, dy
        integer, intent(in) :: nb
        real, intent(in) :: u_now(:, :), v_now(:, :), u_tmp(:, :), v_tmp(:, :)
        real, intent(out) :: u_new(:, :), v_new(:, :)

        integer :: i, j
        real :: adv_u, adv_v, diff_u, diff_v

        do j = nb + 1, size(u_now, 2) - nb
            do i = nb + 1, size(u_now, 1) - nb
                call advection(i, j, dx, dy, u_tmp, v_tmp, adv_u, adv_v)
                call diffusion(i, j, nu, dx, dy, u_tmp, v_tmp, diff_u, diff_v)
                u_new(i, j) = u_now(i, j) - dt * (adv_u - diff_u)
                v_new(i, j) = v_now(i, j) - dt * (adv_v - diff_v)
            end do
        end do
    end subroutine forward_euler
end module time