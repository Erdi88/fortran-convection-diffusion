program Obligatory2026
    implicit none
    integer::N,i
    double precision:: deltat,AA,BB,nu,alpha
    double precision, allocatable:: d(:),a(:),b(:),c(:),u(:)

!Added these
! print counter (every 50 steps)
    integer:: step
!grid spacing = 1/(N-1), time, max change, grid position
    double precision:: dx, time, eps, maxdev
!u from previous time step, initial condition for output, analytical solution
    double precision, allocatable:: uold(:), uinit(:), analytic(:), x(:)

    open(10,file="file01",status="old",action="read")
    read(10,*) deltat,AA,BB,nu,alpha,N
    close(10) 

!Allocate all arrays to size N (grid points)
allocate(u(N), uold(N), uinit(N), analytic(N), x(N))
allocate(a(N), b(N), c(N), d(N))

!Initialise grid spacing and time loop variables
dx = 1.d0/(N-1)
time = 0.d0
step = 0
eps = 2.d-6


! Call BOUNDARY to assign initial and boundary conditions to u
    call BOUNDARY(N,u, AA, BB)
    uinit = u 
! Test the assignment
    write(*,*) 'N=',N,'time step=',deltat,' A=',AA,'B=',BB, &
                   ' nu =',nu,' alpha=',alpha
    write(*,*) 'Initial conditions for u:'
    write(*,300) (u(i),i=1,N)
! Call subroutine THOMAS with N as arrays' dimension and u,a,b,c,d as arguments
! NB! u contains old values that will be updated by the subroutine

! Writing to file_step50 to get all outputs for every 50 time steps
open(10,file="file_step50", status="replace", action="write")
    do while (time < 250.d0)

        uold = u
        time = time + deltat
        step = step + 1

        ! Boundary rows
        a(1)=0.d0; b(1)=1.d0; c(1)=0.d0; d(1)=AA
        a(N)=0.d0; b(N)=1.d0; c(N)=0.d0; d(N)=BB

        !Time derivative uses backward Euler, while x-derivatives use central differences.
        do i=2,N-1
            a(i) = -alpha/dx**2 - nu/(2.d0*dx)
            b(i) =  1.d0/deltat + 2.d0*alpha/dx**2
            c(i) = -alpha/dx**2 + nu/(2.d0*dx)
            d(i) =  uold(i)/deltat
        end do

        call THOMAS(N,u,a,b,c,d)

        maxdev = 0.d0
        do i=1,N
            if (abs(u(i)-uold(i)) > maxdev) maxdev = abs(u(i)-uold(i))
        end do
        if (maxdev < eps) exit

        ! Print every 50 steps
        if (mod(step,50)==0) then
            write(10,200) step, time, maxdev
            write(10,300) (u(i),i=1,N)
            write(10,*)
        end if
    end do
close(10)

    ! Fill x array and call analytical subroutine
    do i=1,N
        x(i) = (i-1)*dx
    end do
    call ANALYTICAL(N, analytic, x, nu, alpha)

    ! Terminal output
    write(*,*)
    if (maxdev < eps) then
        write(*,*) 'Steady state reached.'
    else
        write(*,*) 'Maximum time reached.'
    end if

    write(*,*) 'Final time =', time, '  Steps =', step
    write(*,*) 'Max deviation =', maxdev
    write(*,*)
    write(*,*) 'Final numerical solution:'
    write(*,300) (u(i),i=1,N)
    write(*,*)
    write(*,*) 'Analytical solution:'
    write(*,300) (analytic(i),i=1,N)

    ! Write to file
    open(10,file="file02",status="replace",action="write")
    write(10,*) '# x   u_initial   u_final   u_analytical'
    do i=1,N
        write(10,400) x(i), uinit(i), u(i), analytic(i)
    end do
    close(10)
    write(*,*) 'Results written to file02'

200 FORMAT('Step=',i6,'  Time=',f10.4,'  MaxDiff=',e12.4)
300 FORMAT(11f10.5)
400 FORMAT(4f14.8)
end program Obligatory2026


! THOMAS FUNCTION

subroutine THOMAS(M,u,a,b,c,d)
! Explicit-shape (old-style) passing array dimensions as separate integer argument
    implicit none
    integer :: M, i
    double precision:: u(M),a(M),b(M),c(M),d(M)
    double precision:: bb(M),dd(M)
! Do some array manipulation ending with updating u in accordance with the Thomas algorithm
!COPY B AND D
    bb = b 
    dd = d
    do i=2, M 
        bb(i) = bb(i) - a(i)/bb(i-1)*c(i-1)   ! new main diagonal
        dd(i) = dd(i) - a(i)/bb(i-1)*dd(i-1)  ! new right hand side
    end do
    u(M) = dd(M)/bb(M)

    ! Back substitution
    do i=M-1,1,-1
        u(i) = (dd(i) - c(i)*u(i+1))/bb(i)
    end do
end subroutine THOMAS


! BOUNDARY CONDITIONS

subroutine BOUNDARY(M,u,A, B)
! Explicit-shape (old-style) passing array dimensions as separate integer argument
    implicit none
    integer :: M, i
    double precision:: u(M),x, A, B
! Assign the initial conditions as the specified exponential function but set the rhs boundsary exactly to zero
    do i=1,M-1
        x = 1./(M-1)*(i-1)
        u(i) = dexp(-15*x*x)
    enddo
    u(M) = B ! Endpoint
    u(1) = A ! Startpoint
end 


! ANALYTICAL SOLUTION

subroutine ANALYTICAL(N, u, x, nu, alpha)
    implicit none
    integer :: N, i
    double precision :: u(N), x(N), nu, alpha
    double precision :: ratio, eR

    ratio = nu/alpha
    eR = dexp(ratio)

    do i = 1, N
        u(i) = (eR - dexp(ratio*x(i))) / (eR - 1.0d0)
    end do
end