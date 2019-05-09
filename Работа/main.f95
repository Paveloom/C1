program main
implicit none

     real(8), allocatable, dimension(:,:) :: A ! Матрица исходных данных
     integer N                                 ! Размер выборки
     integer N_if                              ! Число исключений
     
     ! Массив индексов-исключений
     integer(4), allocatable, dimension(:) :: N_if_array

     ! Коррелограмма:

     real(8), allocatable, dimension(:) :: r ! Вектор коэффициентов корреляции
     real(8) x_mean                          ! Среднее значение x
     real(8) koef_mean                       ! Коэффициент 1/N

     ! Периодограмма:

     real(8), allocatable, dimension(:) :: I_p ! Вектор значений периодограммы
     real(8) p_step                            ! Шаг дискретизации частот
     integer p_num                             ! Число желаемых данных
     real(8) p_cur                             ! Текущее значение p

     real(8) koef                 ! Коэффициент для тригонометрических функций
     real(8) cos_value, sin_value ! Вычисленные значения косинуса и синуса
     real(8) diff                 ! Разность (A(t,2) - x_mean)

     ! Общие переменные

     real(8) s1, s2              ! Временные суммы для элементов выражений
     real(8) pi                  ! Число pi
     real(8) k_d, N_d, p_d, t_d  ! Перевод integer в double real
     integer k, t, i, p, p_num_d ! Вспомогательные переменные 

     ! Указать размер выборки
     N = 5860
     
     ! Указать число исключений
     N_if = 301
     
     allocate(N_if_array(6))
     
     N_if_array = (/ 2, (i, i = 4,8) /)
     
     write(*,*) N_if_array

     allocate(A(1:N,2)) ! Исходные данные
     allocate(r(1:N-1)) ! Для вычисления коррелограммы

     ! Считывание исходных данных

        do i = 1, N
                read(*,*) A(i,1), A(i,2)
        enddo

     ! Вычисление среднего значения выборки

     x_mean = 0d0

     N_d = N
     koef_mean = 1d0/(N_d - N_if)

     do t = 1, N
     
          if ((t .eq. 2) .or. ((t .le. 8) .and. (t .ge. 4))) cycle
          if ((t .eq. 10) .or. (t .eq. 54) .or. (t .eq. 133)) cycle
          if ((t .ge. 1541) .and. (t .le. 1546)) cycle
          if ((t .ge. 2141) .and. (t .le. 2149)) cycle
          if ((t .ge. 2425) .and. (t .le. 2426)) cycle
          if ((t .ge. 2772) .and. (t .le. 2777)) cycle
          if ((t .ge. 2807) .and. (t .le. 2809)) cycle
          if ((t .ge. 2863) .and. (t .le. 2867)) cycle
          if ((t .ge. 2896) .and. (t .le. 2897)) cycle
          if (t .eq. 3004) cycle
          if ((t .ge. 3117) .and. (t .le. 3123)) cycle
          if ((t .ge. 3537) .and. (t .le. 3545)) cycle
          if ((t .ge. 3551) .and. (t .le. 3556)) cycle
          if ((t .ge. 3537) .and. (t .le. 3545)) cycle
          if ((t .ge. 3586) .and. (t .le. 3594)) cycle
          if ((t .ge. 3602) .and. (t .le. 3607)) cycle
          if ((t .ge. 3795) .and. (t .le. 3801)) cycle
          if ((t .ge. 3810) .and. (t .le. 3953)) cycle
          if ((t .ge. 3961) .and. (t .le. 4026)) cycle
          if ((t .ge. 3537) .and. (t .le. 3545)) cycle
          if ((t .eq. 5284) .or. (t .eq. 5287)) cycle
          if ((t .eq. 5613) .or. (t .eq. 5669)) cycle 
          
          x_mean = x_mean + A(t,2)
     
     enddo

     x_mean = koef_mean * x_mean

     ! Вычисление коррелограммы

!     open(10, file="output1")
!     do k = 1, N - 1
! 
!          k_d = k
!                
!          s1 = 0d0
!  
!          do t = 1, N - k
!          
!               s1 = s1 + (A(t,2) - x_mean) * (A(t+k,2) - x_mean)
!          
!          enddo

!          s2 = 0d0

!          do t = 1, N
!               
!               s2 = s2 + (A(t,2) - x_mean) * (A(t,2) - x_mean)
!          
!          enddo

!          r(k) = s1/s2

!          write(10,'(e28.20, 1x, e28.20)') k_d, r(k)

!     enddo
!     close(10)

     ! Определение pi
     pi = 4d0*datan(1d0)

     ! Вычисление периодограммы

     p_num = 58600
     p_num_d = p_num
     p_step = N_d / p_num_d

     allocate(I_p(1:p_num))   ! Для вычисления периодограммы вне зависимости от c

!     open(11, file="output2")
!     do p = 1, p_num, 1

!          p_d = p
!          p_cur = 0d0 + p_d * p_step

!          s1 = 0d0
!          s2 = 0d0

!          do t = 1, N

!                    if ((t .eq. 2) .or. ((t .le. 8) .and. (t .ge. 4))) cycle
!                    if ((t .eq. 10) .or. (t .eq. 54) .or. (t .eq. 133)) cycle
!                    if ((t .ge. 1541) .and. (t .le. 1546)) cycle
!                    if ((t .ge. 2141) .and. (t .le. 2149)) cycle
!                    if ((t .ge. 2425) .and. (t .le. 2426)) cycle
!                    if ((t .ge. 2772) .and. (t .le. 2777)) cycle
!                    if ((t .ge. 2807) .and. (t .le. 2809)) cycle
!                    if ((t .ge. 2863) .and. (t .le. 2867)) cycle
!                    if ((t .ge. 2896) .and. (t .le. 2897)) cycle
!                    if (t .eq. 3004) cycle
!                    if ((t .ge. 3117) .and. (t .le. 3123)) cycle
!                    if ((t .ge. 3537) .and. (t .le. 3545)) cycle
!                    if ((t .ge. 3551) .and. (t .le. 3556)) cycle
!                    if ((t .ge. 3537) .and. (t .le. 3545)) cycle
!                    if ((t .ge. 3586) .and. (t .le. 3594)) cycle
!                    if ((t .ge. 3602) .and. (t .le. 3607)) cycle
!                    if ((t .ge. 3795) .and. (t .le. 3801)) cycle
!                    if ((t .ge. 3810) .and. (t .le. 3953)) cycle
!                    if ((t .ge. 3961) .and. (t .le. 4026)) cycle
!                    if ((t .ge. 3537) .and. (t .le. 3545)) cycle
!                    if ((t .eq. 5284) .or. (t .eq. 5287)) cycle
!                    if ((t .eq. 5613) .or. (t .eq. 5669)) cycle
!     
!                    t_d = t
!                             
!                    koef = 2d0 * pi * p_cur * t_d / (N_d - N_if)
!!                    koef = 2d0 * pi * p_cur * t_d / N_d
!                             
!                    cos_value = dcos(koef)
!                    sin_value = dsin(koef)
!                             
!                    ! Проверка на ошибку округления sin_value и cos_value
!                    if (abs(cos_value) .le. 1e-3) cos_value = 0d0
!                    if (abs(sin_value) .le. 1e-3) sin_value = 0d0       
!                                                             
!                    diff = (A(t,2) - x_mean)
!                             
!                    s1 = s1 + diff * cos_value
!                    s2 = s2 + diff * sin_value
!                             
!          enddo

!          I_p(p) = (s1 * s1 + s2 * s2)/( (N_d - N_if) * pi)
!!          I_p(p) = (s1 * s1 + s2 * s2)/(N_d * pi)

!          write(11,'(e16.7, 1x, e16.7, 1x, e16.7)') N_d/p_cur, p_cur/N_d, I_p(p)

!     enddo
!     close(11)

     deallocate(A, r, I_p, N_if_array)

end
