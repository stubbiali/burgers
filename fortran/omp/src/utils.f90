module utils
    implicit none

    private
    public copy_fields, print_field

contains

    subroutine copy_fields(u_src, v_src, u_trg, v_trg)
        implicit none

        real, intent(in) :: u_src(:, :), v_src(:, :)
        real, intent(out) :: u_trg(:, :), v_trg(:, :)

        integer :: i, j

        !$omp do schedule(static) private(i, j)
        do j = 1, size(u_src, 2)
            do i = 1, size(u_src, 1)
                u_trg(i, j) = u_src(i, j)
                v_trg(i, j) = v_src(i, j)
            end do
        end do
        !$omp end do
    end subroutine copy_fields

    subroutine print_field(field)
        implicit none

        real, intent(in) :: field(:, :)

        integer :: i, j

        do i = 1, size(field, 1)
            print "(40f6.2)", field(i, 1:size(field, 2))
        end do
    end subroutine print_field
end module utils