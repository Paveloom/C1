module subprograms
implicit none

     contains
     
     ! [Формирование массива N_index_array и считывание исходных данных]
     subroutine F0_get_index_array(N_index_array, N_wif, N, N_if_array, A)
     implicit none
     
     integer(4), intent(in) :: N_wif ! Размер выборки с исключениями (N - N_if)
     integer(4), intent(in) :: N     ! Размер выборки

     integer(4), intent(out) :: N_Index_array(0:N_wif-1) ! Массив индексов с учётом исключений
     
     ! Вспомогательные переменные
     integer(4) i, k              
     
     ! Опциональные переменные
     integer(4), optional, intent(inout) :: N_if_array(N_wif) ! Массив индексов-исключений
     real(8), optional, intent(in)       :: A(0:N-1,2) ! Матрица исходных данных
     
     if (present(N_if_array)) then
     
          ! Заполнение массива индексов-исключений
          N_if_array = (/ 2, (i, i = 4,8), 10, 54, 133, (i, i = 1541,1546), (i, i = 2141, 2149),&
          & 2425, 2426, (i, i = 2772, 2777), (i, i = 2807, 2809), (i, i = 2863, 2867), 2896, 2897,&
          & 3004, (i, i = 3117, 3123), (i, i = 3537, 3545), (i, i = 3551, 3556), (i, i = 3586, 3594),&
          & (i, i = 3602, 3607), (i, i = 3795, 3801), (i, i = 3810, 3953), (i, i = 3961, 4026),&
          & 5284, 5287, 5613, 5669 /)
     
          k = 1 ! Сдвиг при обнаружении элемента из N_if_array,
                ! уменьшаем таким образом массив индексов 1:N до размера N - N_if

          ! Заполнение массива N_index_array

          do i = 0, N - 1

               if (i .ne. N_if_array(k) - 1) then

                    N_index_array(i - k + 1) = i

               else

                    k = k + 1; cycle

               endif

          enddo
     
     else
     
          k = 1 ! Сдвиг при обнаружении нулевого элемента из столбца A(:,2),
                ! уменьшаем таким образом массив индексов 1:N до размера N - N_if
     
          do i = 0, N - 1

               if (A(i,2) .ne. 0d0) then

                    N_index_array(i - k + 1) = i

               else

                    k = k + 1; cycle

               endif

          enddo
     
     endif
     
     end subroutine
     

     ! [Вычисление среднего значения выборки]
     subroutine F1_mean(A, x_mean, N_index_array, N_wif, N, N_d, use_if)
     implicit none
     
     integer(4), intent(in) :: N_wif                    ! Размер выборки с исключениями (N - N_if)
     integer(4), intent(in) :: N_Index_array(0:N_wif-1) ! Массив индексов с учётом исключений
     integer(4), intent(in) :: N                        ! Размер выборки
     integer(4), intent(in) :: use_if                   ! Использовать массив исключений?
     real(8), intent(in)    :: A(0:N-1,2)               ! Матрица исходных данных
     real(8), intent(in)    :: N_d                      ! Овеществление N
     real(8), intent(out)   :: x_mean                   ! Среднее значение выборки с учётом исключений
     
     integer(4) t, j ! Вспомогательные переменные
     
     x_mean = 0d0
     
     if (use_if .eq. 0) then

          do t = 0, N_wif - 1

               j = N_index_array(t)

               x_mean = x_mean + A(j,2)

          enddo
     
          x_mean = x_mean / N_wif
     
     else
     
          x_mean = sum(A(:,2))
          x_mean = x_mean / N_d
     
     endif

     end subroutine
     
     
     ! [Вычисление коррелограммы прямым способом]
     subroutine F2_correlogram_direct(A, x_mean, N, N_d, bias_fix)
     implicit none
     
     integer(4), intent(in) :: N          ! Размер выборки
     integer(4), intent(in) :: bias_fix   ! Приводить коэффициенты корреляции к несмещённой оценке?
     real(8), intent(in)    :: A(0:N-1,2) ! Матрица исходных данных
     real(8), intent(in)    :: x_mean     ! Среднее значение выборки с учётом исключений
     real(8), intent(in)    :: N_d        ! Овеществление N
     
     integer(4) t, k ! Вспомогательные переменные
     real(8) k_d     ! Овеществление k
     real(8) s1, s2  ! Временные держатели сумм для элементов выражений
     real(8) diff    ! Разность (A(t,2) - x_mean)
     
     real(8) r(0:N-1) ! Вектор значений коррелограммы
     
     ! Вычисление коэффициента автокорреляции c(0)
     
     s2 = 0d0
     
     do t = 0, N - 1

          diff = (A(t,2) - x_mean)

          s2 = s2 + diff * diff

     enddo
     
     ! Вычисление коэффициентов автокорреляции c(k) по 
     ! формуле 4.1 из Chatfield - The Analysis of 
     ! Time Series, стр. 49, и деление их на значение
     ! коэффициента c(0) согласно формуле 4.4
     
     ! Нормировка для несмещённой оценки взята из
     ! Витязева - Спектрально-корреляционный анализ 
     ! равномерных рядов
     
     open(10, file="direct_correlogram.dat")
     write(10, '(a, /, a, 5x, a, /, a)') '#', '# Time Lag, Days, k', 'Autocorrelation Coefficient, r(k)', '#'
     
     if (bias_fix .eq. 0) then ! Приводить коэффициенты корреляции к несмещённой оценке?
     
          do k = 0, N - 1

               k_d = k

               s1 = 0d0

               do t = 0, N - k - 1

                    s1 = s1 + (A(t,2) - x_mean) * (A(t+k,2) - x_mean)

               enddo

               r(k) = s1 * N_d / (N_d - k_d) / s2

               write(10, '(f11.1, 17x, e23.15)') k_d, r(k)

          enddo
          close(10)
     
     else
     
          do k = 0, N - 1

               k_d = k

               s1 = 0d0

               do t = 0, N - k - 1

                    s1 = s1 + (A(t,2) - x_mean) * (A(t+k,2) - x_mean)

               enddo

               r(k) = s1 / s2

               write(10, '(f11.1, 17x, e23.15)') k_d, r(k)

          enddo
          close(10)
     
     endif
     
     end subroutine


     ! [Вычисление периодограммы]
     subroutine F3_periodogram(A, leftbound, leftbound_d, rightbound, p_step, N_index_array, N, N_d, N_wif, x_mean, pi, I_p, use_if)
     implicit none
     
     integer(4), intent(in) :: leftbound, rightbound    ! Границы рабочего диапазона частот
     integer(4), intent(in) :: N                        ! Размер выборки
     integer(4), intent(in) :: N_wif                    ! Размер выборки с исключениями (N - N_if)
     integer(4), intent(in) :: use_if                   ! Использовать массив исключений?
     integer(4), intent(in) :: N_index_array(0:N_wif-1) ! Массив индексов с учётом исключений
     
     real(8), intent(in)    :: A(0:N-1,2)  ! Матрица исходных данных
     real(8), intent(in)    :: p_step      ! Шаг дискретизации множителей p (для частот)
     real(8), intent(in)    :: x_mean      ! Среднее значение выборки с учётом исключений
     real(8), intent(in)    :: pi          ! Число pi
     real(8), intent(in)    :: N_d         ! Овеществление N
     real(8), intent(in)    :: leftbound_d ! Овеществление leftbound
     
     real(8), intent(out) :: I_p(leftbound:rightbound) ! Массив значений периодограммы
     
     real(8) p_d, j_d ! Овеществления
     
     integer(4) p, t, j ! Вспомогательные переменные
     
     real(8) s1, s2    ! Временные держатели сумм для элементов выражений
     real(8) diff      ! Разность (A(j,2) - x_mean)
     real(8) p_cur     ! Текущее значение p
     real(8) arg       ! Аргумент тригонометрических функций
     real(8) cos_value ! Значение косинуса от аргумента arg
     real(8) sin_value ! Значение синуса от аргумента arg
     
     open(11, file="periodogram.dat")
     write(11, '(a, /, a, 6x, a, 6x, a, /, a)') '#', '#    Period, Days, T', 'Frequency, 1/Days, v', 'Periodogram, I(v)', '#'
     
     if (use_if .eq. 0) then ! Использовать массив исключений?
     
          do p = leftbound, rightbound, 1

               p_d = p
               p_cur = leftbound_d - 1d0 + p_d * p_step
          
               s1 = 0d0
               s2 = 0d0

               do t = 0, N_wif - 1

                    ! Применение массива исключений:
                    ! суммирование происходит только 
                    ! по невырожденным элементам
               
                    j = N_index_array(t)

                    j_d = j
                    arg = 2d0 * pi * p_cur * j_d / N_wif

                    cos_value = dcos(arg)
                    sin_value = dsin(arg)

                    ! Проверка на ошибку округления sin_value и cos_value
                    if (abs(cos_value) .le. 1e-3) cos_value = 0d0
                    if (abs(sin_value) .le. 1e-3) sin_value = 0d0

                    diff = (A(j,2) - x_mean)

                    s1 = s1 + diff * cos_value
                    s2 = s2 + diff * sin_value

               enddo
          
               ! Вычисление периодограммы по формуле 7.17 из
               ! Chatfield - The Analysis of Time Series, стр. 111,
               ! с предварительным центрированием ряда (то есть с
               ! вычитанием среднего из выборки)

               I_p(p) = (s1 * s1 + s2 * s2) / (N_wif * pi)
          
               if (abs(I_p(p)) .le. 1e-15) I_p(p) = 0d0

               write(11, '(e23.15, 1x, e23.15, 1x, e23.15)') N_d / p_cur, p_cur / N_d, I_p(p)

          enddo
          close(11)
     
     else
     
          do p = leftbound, rightbound, 1

               p_d = p
               p_cur = leftbound_d - 1d0 + p_d * p_step
          
               s1 = 0d0
               s2 = 0d0

               do t = 0, N - 1
     
                    ! Применение массива исключений:
                    ! суммирование происходит только 
                    ! по невырожденным элементам
               
                    j_d = t
                    arg = 2d0 * pi * p_cur * j_d / N_d

                    cos_value = dcos(arg)
                    sin_value = dsin(arg)

                    ! Проверка на ошибку округления sin_value и cos_value
                    if (abs(cos_value) .le. 1e-3) cos_value = 0d0
                    if (abs(sin_value) .le. 1e-3) sin_value = 0d0

                    diff = (A(t,2) - x_mean)

                    s1 = s1 + diff * cos_value
                    s2 = s2 + diff * sin_value

               enddo
          
               ! Вычисление периодограммы по формуле 7.17 из
               ! Chatfield - The Analysis of Time Series, стр. 111,
               ! с предварительным центрированием ряда (то есть с
               ! вычитанием среднего из выборки)

               I_p(p) = (s1 * s1 + s2 * s2) / (N_d * pi)
          
               if (abs(I_p(p)) .le. 1e-15) I_p(p) = 0d0

               write(11, '(e23.15, 1x, e23.15, 1x, e23.15)') N_d / p_cur, p_cur / N_d, I_p(p)

          enddo
          close(11)
     
     endif
     
     end subroutine
     
     
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
     
     real(8) p_d, t_d, t_koef_d ! Овеществления
     
     integer(4) p, t   ! Вспомогательные переменные
     real(8) s1, s2    ! Временные держатели сумм для элементов выражений
     real(8) p_cur     ! Текущее значение p
     real(8) t_cur     ! Текущее значение t
     real(8) arg       ! Аргумент тригонометрических функций
     real(8) cos_value ! Значение косинуса от аргумента arg
     
     real(8) C(0:t_koef*(N - 1)) ! Вектор значений коррелограммы
     
     ! Вычисление коэффициента автокорреляции c(0)
     
     s2 = 0d0
     s2 = sum(I_p)
     
     ! Вычисление коэффициентов автокорреляции c(k)
     ! и деление их на значение коэффициента c(0)
     
     t_koef_d = t_koef
     
     open(12, file = 'reverse_correlogram.dat')
     write(12, '(a, /, a, 5x, a, /, a)') '#', '# Time Lag, Days, k', 'Autocorrelation Coefficient, r(k)', '#'
     
     if (bias_fix .eq. 0) then ! Приводить коэффициенты корреляции к несмещённой оценке?
     
          do t = 0, t_koef * (N - 1), 1
     
               t_d = t
               t_cur = 0d0 + t_d / t_koef_d
          
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
               t_cur = 0d0 + t_d / t_koef_d
               
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
     
     end subroutine
     
end module
