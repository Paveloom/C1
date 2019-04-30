program main
implicit none

real(8), allocatable, dimension(:,:) :: A ! Матрица исходных данных
real(8), allocatable, dimension(:) :: r ! Вектор коэффициентов корреляции
real(8) s1, s2, s3, s4, s5 ! Суммы для вычисления выборочного коэффициента корреляции
real(8) koef               ! Коэффициент в вычислении ^
integer N                  ! Размер выборки
integer i, l               ! Вспомогательные переменные
real(8) N_d, l_d           ! Перевод integer в double real

! Указать размер выборки
N = 5860

allocate(A(0:N-1,2), r(0:N-1))

do i=0, N-1
read(*,*) A(i,1), A(i,2)
enddo

! Коррелограмма: l - лаг
open(10, file="output1")
do l=0, N-1, 1
   
   s1 = 0; s2 = 0; s3 = 0; s4 = 0; s5 = 0
   do i = 0, N-l-1
   s1 = s1 + A(i,2)*A(i+l,2)
   s2 = s2 + A(i,2)
   s3 = s3 + A(i+l,2)
   s4 = s4 + A(i,2)*A(i,2)
   s5 = s5 + A(i+l,2)*A(i+l,2)
   enddo

   N_d = N; l_d = l

   koef = N_d-l_d
   r(l) = (koef*s1-s2*s3)/(dsqrt(koef*s4-s2*s2)*dsqrt(koef*s5-s3*s3))

 ! write(*,'(e26.20, 1x, e26.20)')  A(m,1), c(m)
   write(10,'(e28.20, 1x, e28.20)') l_d, abs(r(l))

enddo
close(10)

deallocate(A, r)

end
