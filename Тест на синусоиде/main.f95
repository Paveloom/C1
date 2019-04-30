program main
implicit none

real(8), allocatable, dimension(:,:) :: A ! Матрица исходных данных
real(8), allocatable, dimension(:) :: r ! Вектор коэффициентов корреляции
real(8) s1, s2, s3, s4, s5 ! Суммы для вычисления выборочного коэффициента корреляции
real(8) koef               ! Коэффициент в вычислении ^
integer N                  ! Размер выборки
integer i, l               ! Вспомогательные переменные
real(8) N_d, l_d           ! Перевод integer в double real

real(8), allocatable, dimension(:) :: r_hat ! Смещённая оценка коррелограммы
real(8), allocatable, dimension(:) :: nu ! Система частот
real(8), allocatable, dimension(:) :: D ! Вектор коэффициентов периодограммы
real(8) s, delta_t, nu_w, pi
integer j, k
real(8) j_d, k_d

real(8) i_r

! Определение pi
pi = 4d0*datan(1d0)

! Указать размер выборки
N = 200

allocate(A(0:N-1,2), r(0:N-1), r_hat(0:N-1), nu(0:N-1), D(0:N-1))

open(12, file='sinus')
do i=0, N-1
   i_r = i
   i_r = (-1.d0 + 0.02*i_r)*pi
   A(i,1) = i+1
   A(i,2) = dsin(i_r)
   write(12, '(e28.20, 1x, e28.20)') A(i,1), A(i,2) 
enddo
close(12)

! Коррелограмма: l - лаг
open(10, file="sinus1")
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

   r(N-1) = 0

 ! write(*,'(e26.20, 1x, e26.20)')  A(m,1), c(m)
   write(10,'(e28.20, 1x, e28.20)') l_d, abs(r(l))

enddo
close(10)

! Смещённая оценка корреляционной функции
do l=0, N-1, 1
N_d=N; l_d=l
r_hat(l) = (N_d-l_d)/N_d*r(l)
enddo

! Шаг выборки
delta_t = A(2,1) - A(1,1)

! Рабочая частота
nu_w = 1d0/delta_t

! Система частот
do j=0, N-1, 1
N_d=N; j_d=j
nu(j) = nu_w/N_d*j_d
enddo

! Периодограмма
open(11, file="sinus2")
do l=0, N-1, 1
   N_d=N
   s=0
      do k=0, N-1, 1
      k_d=k
      s=s+r_hat(k)*dcos(2d0*pi*nu(l)*k_d*delta_t)
      enddo
   D(l) = 1d0/N_d*(2d0*s-r_hat(0))
 ! write(*,*)  A(l,1), D(l)
   write(11,*) delta_t*l, D(l)
enddo
close(11)

deallocate(A, r, r_hat, nu, D)

end
