program main
implicit none

real(8), allocatable, dimension(:,:) :: A ! Матрица исходных данных
integer N                                 ! Размер выборки

! Коррелограмма:

real(8), allocatable, dimension(:) :: r ! Вектор коэффициентов корреляции
real(8) x_mean                          ! Среднее значение x
real(8) koef_mean                       ! Коэффициент 1/N

! Периодограмма:

real(8), allocatable, dimension(:) :: I_p, I_p_2 ! Вектор значений периодограммы
real(8) p_step                                   ! Шаг дискретизации частот
integer p_num                                    ! Число желаемых данных
real(8) p_cur                                    ! Текущее значение p

real(8), allocatable, dimension(:) :: c ! Коэффициенты автокорреляции

! Общие переменные

real(8) s1, s2, s3          ! Временные суммы для элементов выражений
real(8) pi                  ! Число pi
real(8) k_d, N_d, p_d, t_d  ! Перевод integer в double real
integer k, t, i, p, p_num_d ! Вспомогательные переменные 

! Указать размер выборки
N = 5860

allocate(A(1:N,2)) ! Исходные данные
allocate(r(1:N-1)) ! Для вычисления коррелограммы
allocate(c(0:N-1))   ! Коэффициенты автокорреляции

! Считывание исходных данных

        do i=1, N
                read(*,*) A(i,1), A(i,2)
        enddo

! Вычисление среднего значения выборки

        x_mean = 0d0

        N_d = N
        koef_mean = 1d0/N_d

        do t=1, N
           x_mean = x_mean + A(t,2)
        enddo
        x_mean = koef_mean * x_mean

! Вычисление коррелограммы

        open(10, file="output1")
        do k = 1, N - 1
 
                k_d = k
                
                s1 = 0d0
  
                do t = 1, N - k
                        s1 = s1 + (A(t,2) - x_mean) * (A(t+k,2) - x_mean)
                enddo

                c(k) = s1/N_d

                s2 = 0d0

                do t = 1, N
                        s2 = s2 + (A(t,2) - x_mean) * (A(t,2) - x_mean)
                enddo

                r(k) = s1/s2

                write(10,'(e28.20, 1x, e28.20)') k_d, r(k)

        enddo
        close(10)

! Определение pi
pi = 4d0*datan(1d0)

! Вычисление c(0)
        c(0) = 0d0
        do t = 1, N
                c(0) = c(0) + (A(t,2) - x_mean) * (A(t,2) - x_mean)
        enddo
        c(0) = c(0) / N_d

! Вычисление периодограммы

p_num = 58600
p_num_d = p_num
p_step = N_d/p_num_d

allocate(I_p(1:p_num-1))   ! Для вычисления периодограммы вне зависимости от c
allocate(I_p_2(1:p_num-1)) ! Для вычисления периодограммы через c

        open(11, file="output2")
        do p=1, p_num-1, 1

                p_d = p
                p_cur = 0d0 + p_d * p_step

                s1 = 0d0
                s2 = 0d0

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

                        t_d = t
                        s1 = s1 + (A(t,2) - x_mean) * dcos(2d0 * pi * p_cur * t_d / N_d)
                        s2 = s2 + (A(t,2) - x_mean) * dsin(2d0 * pi * p_cur * t_d / N_d)

                enddo

                s3 = 0d0
             
                do k = 1, N - 1
                        k_d = k
                        s3 = s3 + c(k) * dcos(k_d * 2d0 * pi * p_cur / N_d)
                enddo

                I_p(p) = (s1 * s1 + s2 * s2)/(N_d * pi)
                I_p_2(p) = (c(0) + 2d0 * s3)/pi

                write(11,'(e28.20, 1x, e28.20, 1x, e28.20, 1x, e28.20)') N_d/p_d, p_d/N_d, I_p(p), I_p_2(p)

        enddo
        close(11)

deallocate(c)
deallocate(A)
deallocate(r)
deallocate(I_p, I_p_2)

end
