! Created by yungpeepo on 16/02/21.

program dmc_mw
    use sub
    implicit none
    real(kind=8) :: t, t2, V, weight, energy, pos
    real(kind=8), allocatable, dimension(:) :: wwalker, el,alpha
    integer, allocatable, dimension(:) :: xxx
    real(kind=8), allocatable, dimension(:,:) :: y
    integer :: L, which, nit, nw, nbra, ij, w, x_new

    call read_input(alpha, t, t2, L, V, which, nit, nw, nbra)

    allocate(xxx(nw), wwalker(nw), el(nw), y(nbra,nw))
    xxx(:) = 1
    open(2,file='dmc_branch.data',status='unknown', action='write')
    write(2,*) "# mean_posi     mean_energy     mean_weight"
    do ij = 1, nit
        call random_number(y)
        !$omp parallel do default(shared) private(w,x_new)
        do w = 1, nw
            call E_dmc(alpha, t, t2, V, L, which, nbra, xxx(w), y(:,w), x_new, wwalker(w), el(w))
            xxx(w) = x_new
        end do
        !$omp end parallel do
        weight = sum(wwalker)
        energy = sum(el * wwalker) / weight
        pos = sum(xxx * wwalker) / weight
        weight = weight / nw
        call branching(xxx, wwalker, nw)
        write(2,*) pos, energy, weight
        if (mod(10*(ij+1),nit) == 0) then
            print*, 100*(ij+1)/nit, '% completed'
        end if
    end do
end program dmc_mw