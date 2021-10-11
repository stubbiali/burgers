module spacetime
    implicit none

    private
    public forward_euler

contains

    subroutine advection(i, j, dx, dy, u, v, adv_u, adv_v)
        implicit none

        integer, intent(in) :: i, j
        real, intent(in) :: dx, dy
        real, intent(in) :: u(:, :), v(:, :)
        real, intent(out) :: adv_u, adv_v

        real :: adv_u_x, adv_u_y, adv_v_x, adv_v_y

        adv_u_x = u(i, j) / (60. * dx) * (                  &
                + 45. * (u(i + 1, j) - u(i - 1, j))         &
                        -  9. * (u(i + 2, j) - u(i - 2, j))         &
                        +       (u(i + 3, j) - u(i - 3, j))         &
                ) - abs(u(i, j)) / (60. * dx) * (           &
                +       (u(i + 3, j) + u(i - 3, j))         &
                        -  6. * (u(i + 2, j) + u(i - 2, j))         &
                        + 15. * (u(i + 1, j) + u(i - 1, j))         &
                        - 20. *  u(i, j)                            &
                )
        adv_u_y = v(i, j) / (60. * dy) * (                  &
                + 45. * (u(i, j + 1) - u(i, j - 1))         &
                        -  9. * (u(i, j + 2) - u(i, j - 2))         &
                        +       (u(i, j + 3) - u(i, j - 3))         &
                ) - abs(v(i, j)) / (60. * dy) * (           &
                +       (u(i, j + 3) + u(i, j - 3))         &
                        -  6. * (u(i, j + 2) + u(i, j - 2))         &
                        + 15. * (u(i, j + 1) + u(i, j - 1))         &
                        - 20. *  u(i, j)                            &
                )
        adv_v_x = u(i, j) / (60. * dx) * (                  &
                + 45. * (v(i + 1, j) - v(i - 1, j))         &
                        -  9. * (v(i + 2, j) - v(i - 2, j))         &
                        +       (v(i + 3, j) - v(i - 3, j))         &
                ) - abs(u(i, j)) / (60. * dx) * (           &
                +       (v(i + 3, j) + v(i - 3, j))         &
                        -  6. * (v(i + 2, j) + v(i - 2, j))         &
                        + 15. * (v(i + 1, j) + v(i - 1, j))         &
                        - 20. *  v(i, j)                            &
                )
        adv_v_y = v(i, j) / (60. * dy) * (                  &
                + 45. * (v(i, j + 1) - v(i, j - 1))         &
                        -  9. * (v(i, j + 2) - v(i, j - 2))         &
                        +       (v(i, j + 3) - v(i, j - 3))         &
                ) - abs(v(i, j)) / (60. * dy) * (           &
                +       (v(i, j + 3) + v(i, j - 3))         &
                        -  6. * (v(i, j + 2) + v(i, j - 2))         &
                        + 15. * (v(i, j + 1) + v(i, j - 1))         &
                        - 20. *  v(i, j)                            &
                )

        adv_u = adv_u_x + adv_u_y
        adv_v = adv_v_x + adv_v_y
    end subroutine advection

    subroutine diffusion(i, j, nu, dx, dy, u, v, diff_u, diff_v)
        implicit none

        integer, intent(in) :: i, j
        real, intent(in) :: dx, dy, nu
        real, intent(in) :: u(:, :), v(:, :)
        real, intent(out) :: diff_u, diff_v

        diff_u = nu * (             &
                (                       &
                        -       u(i - 2, j) &
                                + 16. * u(i - 1, j) &
                                - 30. * u(i, j)     &
                                + 16. * u(i + 1, j) &
                                -       u(i + 2, j) &
                        ) / (12. * dx * dx) +   &
                        (                       &
                                -       u(i, j - 2) &
                                        + 16. * u(i, j - 1) &
                                        - 30. * u(i, j)     &
                                        + 16. * u(i, j + 1) &
                                        -       u(i, j + 2) &
                                ) / (12. * dy * dy)     &
                )
        diff_v = nu * (             &
                (                       &
                        -       v(i - 2, j) &
                                + 16. * v(i - 1, j) &
                                - 30. * v(i, j)     &
                                + 16. * v(i + 1, j) &
                                -       v(i + 2, j) &
                        ) / (12. * dx * dx) +   &
                        (                       &
                                -       v(i, j - 2) &
                                        + 16. * v(i, j - 1) &
                                        - 30. * v(i, j)     &
                                        + 16. * v(i, j + 1) &
                                        -       v(i, j + 2) &
                                ) / (12. * dy * dy)     &
                )
    end subroutine diffusion

    subroutine forward_euler(nu, dt, dx, dy, nb, u_now, v_now, u_tmp, v_tmp, u_new, v_new)
        implicit none

        real, intent(in) :: nu, dt, dx, dy
        integer, intent(in) :: nb
        real, intent(in) :: u_now(:, :), v_now(:, :), u_tmp(:, :), v_tmp(:, :)
        real, intent(out) :: u_new(:, :), v_new(:, :)

        integer :: i, j
        real :: adv_u, adv_v, diff_u, diff_v

        !$omp do schedule(static) private(i, j, adv_u, adv_v, diff_u, diff_v)
        do j = nb + 1, size(u_now, 2) - nb
            do i = nb + 1, size(u_now, 1) - nb
                call advection(i, j, dx, dy, u_tmp, v_tmp, adv_u, adv_v)
                call diffusion(i, j, nu, dx, dy, u_tmp, v_tmp, diff_u, diff_v)
                u_new(i, j) = u_now(i, j) - dt * (adv_u - diff_u)
                v_new(i, j) = v_now(i, j) - dt * (adv_v - diff_v)
            end do
        end do
        !$omp end do
    end subroutine forward_euler
end module spacetime