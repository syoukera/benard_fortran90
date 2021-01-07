program main
    use params
    implicit none

    character*12 :: hedu, hedv, hedp
    data hedu /'U-VELOCITY  '/
    data hedv /'V-VELOCITY  '/
    data hedp /'PRESSURE    '/

    ! Chapter 1, paramenters and control indices
    integer :: niter = 0
    ! program control and monitor
    integer :: maxit = 1000
    integer :: imon = 5
    integer :: jmon = 5
    real(16) :: sormax = 1.0e-1
    real(16) :: sorce

    integer :: i, j
    real(16) :: arden, flowin, xmonin, tin

    ! Chapter 2, initial operations
    call grid
    call init
    flowin = 0.0
    xmonin = 0.0

    ! Chapter 3, iteration loop
    sorce = 1.0e10
    do while (sorce .GT. sormax)
        niter = niter + 1
        call calcu
        call calcv
        call calcp
        call calct

        flowin = 0.0
        xmonin = 0.0
        tin    = 0.0
        if (vmax .LT. vmax1) vmax = vmax1
        do i = 2, nim1
            arden = densit*sew(i)
            flowin = flowin + arden*vmax
            xmonin = xmonin + arden*vmax*vmax
            tin = tin + arden*vmax*sph*(tsouth-tnorth)
        end do

        resorm = resorm/flowin
        resoru = resoru/xmonin
        resorv = resorv/xmonin
        resort = resort/tin

        print '(i3, 8(E10.2))', niter, resoru, resorv, resorm, resort, &
                U(imon, jmon), V(imon, jmon), P(imon, jmon), T(imon, jmon)

        ! termination tests
        sorce = amax1(resorm, resoru, resorv, resort)
        if (niter .EQ. 20 .AND. sorce .GT. 1.0e7*sormax) then
            print *, 'Solution diverges or needs more iterations'
            stop
        end if

        ! if (niter .EQ. 2) exit
    end do

    ! Chapter 4, Final operations and output
    open(10, file='output.dat')

    do i = 1, ni
        do j = 1, nj
            write(10, *) X(i), Y(j), U(i, j), V(i, j), P(i, j), T(i, j) 
        end do
    end do

end program main