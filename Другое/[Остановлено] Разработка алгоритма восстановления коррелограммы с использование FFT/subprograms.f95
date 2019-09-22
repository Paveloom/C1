module subprograms
implicit none

     contains
     
     ! [Функция для вычисления логарифма с основанием 2]
     real(8) function dlog2(x)
     implicit none
      
     real(8), intent(in) :: x

     dlog2 = dlog(x) / dlog(2d0)
     
     end function
          
     ! [Вычисление коррелограммы через применение обратного преобразования Фурье к периодограмме]
     subroutine F4_correlogram_fourier_transform(I_p, leftbound, leftbound_d, rightbound, p_step, t_koef, N, N_d, &
     &pi, bias_fix)
     implicit none
     
     integer(4), intent(in) :: N                         ! Размер выборки
     integer(4), intent(in) :: t_koef                    ! Множитель дискретизации множителей t (для периодов)
     integer(4), intent(in) :: leftbound, rightbound     ! Границы рабочего диапазона частот
     integer(4), intent(in) :: bias_fix                  ! Приводить коэффициенты корреляции к несмещённой оценке?
     real(8), intent(in)    :: I_p(leftbound:rightbound) ! Массив значений периодограммы
     real(8), intent(in)    :: p_step                    ! Шаг дискретизации множителей p (для частот)
     real(8), intent(in)    :: N_d                       ! Овеществление N
     real(8), intent(in)    :: leftbound_d               ! Овеществление leftbound
     real(8), intent(in)    :: pi                        ! Число pi
     
     real(8) p_d, t_d, t_koef_d, rightbound_d ! Овеществления
     
     integer(4) p, t, ier ! Вспомогательные переменные
     real(8) s1, s2       ! Временные держатели сумм для элементов выражений
     
     real(8) p_cur  ! Текущее значение p
     real(8) t_cur  ! Текущее значение t
     real(8) t_step ! Шаг дискретизации множителей t (для периодов)
     
     real(8) arg       ! Аргумент тригонометрических функций
     real(8) cos_value ! Значение косинуса от аргумента arg
     
     real(8) C(0:t_koef*(N - 1)) ! Вектор значений коррелограммы
     
     ! Массив значений периодограммы, дополненный нулями
     complex(8), allocatable, dimension(:) :: I_ext
     
     real(8) ext          ! Держатель показателя log_2(x)
     integer(4) ext_new   ! Правая граница для I_ext
     logical(1) ext_l     ! Равен ли ext величине aint(ext)
     
     ! Где x = rightbound_d - leftbound_d + 1
     
     ! Вычисление коэффициента автокорреляции c(0)
     
     s2 = 0d0
     s2 = sum(I_p)
     
     ! Проверка на представимость размера массива I_p в виде 2**m
     
     rightbound_d = rightbound
     
     ext = dlog2(rightbound_d - leftbound_d + 1)
     ext_l = ext .ne. aint(ext)
     
     if (ext_l) then
     
          ext_new = 2**(floor(ext) + 1) + leftbound - 1
     
          ! Массив значений периодограммы, дополненный нулями
          allocate(I_ext(leftbound:ext_new), stat = ier)
          if (ier .ne. 0) stop 'Не удалось выделить память для массива I_ext'
          
          I_ext(leftbound:rightbound) = I_p(leftbound:rightbound)
          I_ext(rightbound + 1:ext_new) = 0d0
          
     endif
     
     ! Вычисление коэффициентов автокорреляции c(k)
     ! и деление их на значение коэффициента c(0)
     
     t_koef_d = t_koef
     t_step = 1d0 / t_koef_d
     
     open(12, file = 'reverse_correlogram.dat')
     write(12, '(a, /, a, 5x, a, /, a)') '#', '# Time Lag, Days, k', 'Autocorrelation Coefficient, r(k)', '#'
     
     if (bias_fix .eq. 0) then ! Приводить коэффициенты корреляции к несмещённой оценке?
     
          do t = 0, t_koef * (N - 1), 1
     
               t_d = t
               t_cur = 0d0 + t_d * t_step
          
               s1 = 0d0
          
               do p = leftbound, rightbound, 1

                    p_d = p
                    p_cur = leftbound_d - 1d0 + p_d * p_step
          
                    arg = 2d0 * pi * p_cur * t_cur / N_d

                    cos_value = dcos(arg)

                    ! Проверка на ошибку округления sin_value и cos_value
                    if (abs(cos_value) .le. 1e-3) cos_value = 0d0
               
                    s1 = s1 + I_p(p) * cos_value

               enddo
          
               ! Вычисление коррелограммы через обратное преобразование Фурье
               ! к периодограмме по формуле для коэффициентов автокорреляции
               ! 73 (вычисляя величину c(k)/c(0)) из Витязева - 
               ! Спектрально-корреляционный анализ равномерных рядов, стр. 25
          
               C(t) = s1 * N_d / (N_d - t_cur) / s2

               if (abs(C(t)) .le. 1e-15) C(t) = 0d0

               write(12, '(f11.1, 17x, e23.15)') t_cur, C(t)
     
          enddo     
          close(12)
     
     else
          
          open(12, file = 'reverse_correlogram.dat')
          do t = 0, t_koef * (N - 1), 1
          
               t_d = t
               t_cur = 0d0 + t_d / t_koef
               
               s1 = 0d0
               
               do p = leftbound, rightbound, 1
     
                    p_d = p
                    p_cur = leftbound_d - 1d0 + p_d * p_step
               
                    arg = 2d0 * pi * p_cur * t_cur / N_d
     
                    cos_value = dcos(arg)

                    ! Проверка на ошибку округления cos_value
                    if (abs(cos_value) .le. 1e-3) cos_value = 0d0
               
                    s1 = s1 + I_p(p) * cos_value

               enddo
          
               ! Вычисление коррелограммы через обратное преобразование Фурье
               ! к периодограмме по формуле для коэффициентов автокорреляции
               ! 73 (вычисляя величину c(k)/c(0)) из Витязева - 
               ! Спектрально-корреляционный анализ равномерных рядов, стр. 25
          
               C(t) = s1 / s2

               if (abs(C(t)) .le. 1e-15) C(t) = 0d0

               write(12, '(f11.1, 17x, e23.15)') t_cur, C(t)
     
          enddo     
          close(12)
     
     endif
     
     if (ext_l) then
     
          call FFT(I_ext, leftbound, leftbound_d, t_step)
          
          open(13, file = 'FFT_reverse_correlogram.dat')
          write(13, '(e23.15)') I_ext
          
          close(13)
          deallocate(I_ext)
     
     endif
     
     end subroutine
     
     
     ! [Процедура быстрого преобразования Фурье]
     recursive subroutine FFT(C, leftbound, leftbound_d, t_step)
     implicit none
     
     complex(8), intent(inout), dimension(leftbound:) :: C ! Входной / выходной массив данных
     integer(4), intent(in) :: leftbound ! Левая граница рабочего диапазона частот
     real(8), intent(in) :: leftbound_d ! Овеществление leftbound
     real(8), intent(in) :: t_step ! Шаг дискретизации множителей t (для периодов)
     
     ! Массивы данных по нечётным и чётным индексам
     complex(8), allocatable, dimension(:) :: odd, even
     
     integer(4) M      ! Размер входного массива
     integer(4) M_half ! Половина от M
     
     real(8) arg   ! Аргумент тригонометрических функций
     real(8) arg_0 ! Начальный аргумент для тригонометрических функций
     
     real(8) cos_value   ! Значение косинуса от аргумента arg
     real(8) sin_value   ! Значение синуса от аргумента arg
     real(8) cos_value_0 ! Значение косинуса от аргумента arg_0
     real(8) sin_value_0 ! Значение синуса от аргумента arg_0
     
     complex(8) w  ! Держатель значений корня из единицы
     complex(8) wn ! Комплексное число cos(koef) + i * sin(koef)
     
     integer(4) ier, i ! Вспомогательные переменные
     real(8) M_d       ! Овеществление M
     
     real(8) pi ! Число pi
     
     ! Определение pi
     pi = 4d0 * datan(1d0)
     
     M = size(C)
     
     if (M .eq. 1) return
     
     M_half = M / 2
     
     allocate(odd(leftbound:leftbound + M_half - 1), stat = ier)
     if (ier .ne. 0) stop 'Не удалось выделить память для массива odd'
     
     allocate(even(leftbound:leftbound + M_half - 1), stat = ier)
     if (ier .ne. 0) stop 'Не удалось выделить память для массива even'
     
     even(leftbound:leftbound + M_half - 1)  = C(leftbound:leftbound + M - 1:2)
     odd(leftbound:leftbound + M_half - 1)   = C(1:leftbound + M - 1:2)
     
     call FFT(even, leftbound, leftbound_d, t_step)
     call FFT(odd, leftbound, leftbound_d, t_step)
     
     M_d = M
     
     arg = 2d0 * pi / M_d
     arg_0 = arg * (leftbound_d * (1d0 + p_step) - 1d0)
     
     cos_value = dcos(arg * t_step)
     sin_value = dsin(arg * t_step)
     
     cos_value_0 = dcos(arg_0)
     sin_value_0 = dsin(arg_0)
     
     ! Проверка на ошибку округления cos_value и sin_value
     if (abs(cos_value) .le. 1e-3) cos_value = 0d0
     if (abs(sin_value) .le. 1e-3) sin_value = 0d0
     
     ! Проверка на ошибку округления cos_value_0 и sin_value_0
     if (abs(cos_value_0) .le. 1e-3) cos_value_0 = 0d0
     if (abs(sin_value_0) .le. 1e-3) sin_value_0 = 0d0
     
     w  = complex(1d0, 0)
     wn = complex(cos_value, sin_value)
     
     do i = 0, leftbound + M_half - 1, 1
     
          ! "Преобразование бабочки"
          C(i) = even(i) + w * odd(i)
          C(i + M_half) = even(i) - w * odd(i)
          
          w = w * wn
     
     enddo
     
     deallocate(odd, even)
     
     end subroutine
     
end module
