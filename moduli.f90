! Created by yungpeepo on 16/02/21
module sub
    implicit none
    public :: psi, update, E_dmc, read_input

contains
    subroutine psi(x,alpha,L,which, out)
        implicit none
        integer, intent(in) :: L, which
        integer, intent(in) :: x
        real(kind=8), dimension(:), intent(in) :: alpha
        real(kind=8), intent(out) :: out

        if (which == 0) then
            if (x<1 .or. x>L) then
                out = 0
            else
                out = 1
            end if
        else if (which == 1) then
            if (x<1 .or. x>L) then
                out = 0
            else
                out = exp( - alpha(1) * x)
            end if
        else if (which == 2) then
            if (x<1 .or. x>L) then
                out = 0
            else
                out = exp(- alpha(1)*x) * x**alpha(2)
            end if
        end if
    end subroutine psi

    subroutine update(x, alpha, t, t2, V, L, which, b_out, elx, P_left, P_right)

        implicit none
        integer, intent(in) :: x, L, which
        real(kind=8), intent(in) :: t, t2, V
        real(kind=8), intent(in), dimension(:) :: alpha
        real(kind=8), intent(out) :: b_out, P_left, P_right, elx
        real(kind=8) :: vsf, bx, Lambda, psix, psi_l,psi_r, psi_l2, psi_r2

        Lambda = V * L
        call psi(x, alpha, L, which, psix)
        call psi(x-1, alpha, L, which, psi_l)
        call psi(x+1, alpha, L, which, psi_r)
        call psi(x-2, alpha, L, which, psi_l2)
        call psi(x+2, alpha, L, which, psi_r2)
        vsf = t2 * (psi_l2 + psi_r2) / psix
        bx = -V * x - vsf + Lambda + t * (psi_l + psi_r) /psix

        elx = -t * (psi_l + psi_r) /psix + vsf + V * x
        P_left = t * psi_l / psix / bx
        P_right = t * psi_r /psix / bx
        b_out = bx / Lambda

    end subroutine update

    subroutine E_dmc(alpha, t, t2, V, L, which, n_it, x0, y, x, w, elx)
        implicit none
        real(kind=8), dimension(:), intent(in) :: alpha
        real(kind=8), intent(in) :: t, t2, V
        real(kind=8), dimension(:), intent(in) :: y
        integer, intent(in) :: L, which, n_it, x0
        real(kind=8), intent(out) :: w, elx
        integer, intent(out) :: x
        real(kind=8) :: bx, P_left, P_right, b_old
        integer :: i

        x = x0
        call update(x, alpha, t, t2, V, L, which, bx, elx, P_left, P_right)
        !weight are set to one
        w = 1
        do i = 1, n_it
            b_old = bx
            if (y(i) < P_left) then
                x = x - 1
                call update(x, alpha, t, t2, V, L, which, bx, elx, P_left, P_right)
            else if (y(i)< P_left + P_right) then
                x = x + 1
                call update(x, alpha, t, t2, V, L, which, bx, elx, P_left, P_right)
            end if
            w = w*b_old
        end do

    end subroutine E_dmc

    subroutine branching(x_w, w_w, nw)
        implicit none
        integer, dimension(:), intent(inout) :: x_w
        real(kind=8), dimension(:), intent(in) :: w_w
        real(kind=8), allocatable, dimension(:) :: prob
        integer, allocatable, dimension(:) :: x_out
        real(kind=8) :: y, prob_sum, z
        integer, intent(in) :: nw
        integer :: i, j

        allocate(prob(nw), x_out(nw))
        prob(:) = w_w(:) / sum(w_w)

        call random_number(y)

        do i = 1, nw
            z = ( i + y - 1)/nw
            j = 1
            prob_sum = prob(1)
            do while (z > prob_sum)
                j = j + 1
                prob_sum = prob_sum + prob(j)
            end do
        x_out(i) = x_w(j)
        end do
        x_w = x_out
    end subroutine branching

    subroutine read_input(alpha, t, t2, L, V, which, nit, nw, nbra)
        implicit none
        real(kind=8), allocatable, dimension(:), intent(out) :: alpha
        real(kind=8), intent(out) :: t
        real(kind=8), intent(out) :: t2
        integer,           intent(out) :: L
        real(kind=8),   intent(out) :: V
        integer,        intent(out) :: which
        integer,        intent(out) :: nit
        integer,        intent(out) :: nw
        integer,        intent(out) :: nbra

        character(len=100) :: buffer, label, equal
        integer :: pos
        integer, parameter :: fh = 15
        integer :: ios = 0
        integer :: line = 0

        open(fh, file='input.txt')

        ! ios is negative if an end of record condition is encountered or if
        ! an endfile condition was detected.  It is positive if an error was
        ! detected.  ios is zero otherwise.

        do while (ios == 0)
            read(fh, '(A)', iostat=ios) buffer
            if (ios == 0) then
                line = line + 1

                ! Find the first instance of whitespace.  Split label and data.
                pos = scan(buffer, ' 	')
                label = buffer(1:pos)
                buffer = buffer(pos+1:)

                select case (label)
                    case("t")
                        read(buffer,*) equal, t
                    case("t2")
                        read(buffer,*) equal, t2
                    case("L")
                        read(buffer,*) equal, L
                    case("V")
                        read(buffer,*) equal, V
                    case("which")
                        read(buffer,*) equal, which
                    case("alpha")
                        if (which == 2) then
                            allocate(alpha(2))
                            read(buffer,*) equal, alpha(1), alpha(2)
                        else
                            allocate(alpha(1))
                            read(buffer,*) equal, alpha(1)
                        end if
                    case("nsteps")
                        read(buffer,*) equal, nit
                    case("nbra")
                        read(buffer,*) equal, nbra
                    case("nw")
                        read(buffer,*) equal,  nw
                    case default
                        ! an unknown word will stop the execution
                        write(0,*) "Unknown keyword :",trim(label)
                        stop
                    end select
            end if
        end do

    end subroutine read_input
end module