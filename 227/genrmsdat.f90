

program genrmsdat

implicit none

integer, parameter :: steps=10000
real(8), parameter :: pi = atan(1.)*4., dt=0.01
integer :: i, j
real(8) :: lx = 4*pi, lz = 2*pi, aor, aoi
real(8) :: w, lambda, B=100, kx = 0.5, kz=1.0, k=sqrt(1.25), Pe=100, Re=1000, uo, to
real(8) :: u(steps), temp(steps), t(steps)

w = -0.5 * sqrt(4.*((B*kx**2)/(k**2) - (k**4)*(1./(Pe**2) + 1./(Re**2) + 2./(Pe*Re))))
lambda = 0.5*k**2.*(1./Pe + 1./Re)
aor = 10.
aoi = (10./w)*(-lambda+(k**2)/Pe)

to = aor/sqrt(2.)
uo = (10.*w + aoi*((k**2)/Pe - lambda))/(kx*sqrt(2.))*k

do i = 1, steps
    t(i) = 0 + (i-1)*dt
end do

temp = abs(to*exp(-lambda*t)*cos(w*t))
u = abs(uo*exp(-lambda*t)*sin(w*t))

print *, (-aoi*w - 10.*lambda + (10.*k**2)/Pe)

open(40, file="rmsdat", status="new")

do i = 1, steps
    ! need to make this write formatted, FT is NOT defined yet
    write (40, '(3E20.7)') t(i), u(i), temp(i)

end do 
close(40)

end program genrmsdat
