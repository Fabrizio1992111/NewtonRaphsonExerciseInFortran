program Newton

!uso de modulos
use math

!declaracion de variables

integer :: ier
real(8),dimension(3) :: x
real(8),dimension(3) :: fx

write(*,'("resuelvo un sistema no lineal 2x2 mediante el metodo de Newton")')
write(*,'(" f(x,y,z) = 2x + y - 3z - 7 = 0 ")')
write(*,'(" g(x,y,z) = 5x - y + z + 19 = 0 ")')
write(*,'(" h(x,y,z) = x - y - 4z -4 = 0 ")')
 

x=(/-1.d0,3.d0,-2.d0/) !estimacion inicial

ier=NewtonRaphson_d(SNL2x2,x,fx,1.d-06,100)

write(*,'(" La solucion es: ",2e12.4)') x(1)
write(*,'(" La solucion es: ",2e12.4)') x(2)
write(*,'(" La solucion es: ",2e12.4)') x(3)
write(*, '(" codigo de error:",i)') ier

pause

CONTAINS

function SNL2x2(x,Jac) result(fx)

!variables externas
real(8),dimension(:) :: x
real(8), dimension(:,:), optional :: Jac
real(8), dimension(size(x)) :: fx

fx(1)= 2.d0 * x(1) + x(2) - 3.d0*x(3) - 7.d0
fx(2)= 5.d0 * x(1) - x(2) + x(3) + 19.d0 
fx(3)= x(1) - x(2) - 4.d0*x(3) - 4.d0 

if (present(Jac)) then
Jac(1,1) = 2.d0
Jac(1,2) = 1.d0
Jac(1,3) = -3.d0
Jac(2,1) = 5.d0
Jac(2,2) = -1.d0
Jac(2,3) = 1.d0
Jac(3,1) = 1.d0
Jac(3,2) = -1.d0
Jac(3,3) = -4.d0
end if

end function

end program