module subprograms
implicit none

     contains

     ! [Вычисление среднего значения выборки]
     subroutine F1_mean(A, x_mean, N_index_array, N_wif, N, N_d, use_if)
     implicit none
     
     integer(4), intent(in) :: N_Index_array(0:N_wif-1) ! Массив индексов с учётом исключений
     integer(4), intent(in) :: N_wif                    ! Размер выборки с исключениями (N - N_if)
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
     
          do t = 0, N - 1

               x_mean = x_mean + A(t,2)

          enddo
          
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
     
     ! Вычисление коэффициентов автокорреляции c(k)
     ! и деление их на значение коэффициента c(0)
     
     if (bias_fix .eq. 0) then
     
          open(10, file="output1")
          do k = 0, N - 1

               k_d = k

               s1 = 0d0

               do t = 0, N - k - 1

                    s1 = s1 + (A(t,2) - x_mean) * (A(t+k,2) - x_mean)

               enddo

               r(k) = s1 * N_d / (N_d - k_d) / s2

               write(10,'(e16.7, 1x, e16.7)') k_d, r(k)

          enddo
          close(10)
     
     else
     
          open(10, file="output1")
          do k = 0, N - 1

               k_d = k

               s1 = 0d0

               do t = 0, N - k - 1

                    s1 = s1 + (A(t,2) - x_mean) * (A(t+k,2) - x_mean)

               enddo

               r(k) = s1 / s2

               write(10,'(e16.7, 1x, e16.7)') k_d, r(k)

          enddo
          close(10)
     
     endif
     
     end subroutine


     ! [Вычисление периодограммы]
     subroutine F3_periodogram(A, leftbound, leftbound_d, rightbound, p_step, N_index_array, N, N_d, N_wif, x_mean, pi, I_p, use_if)
     implicit none
     
     real(8), intent(out) :: I_p(leftbound:rightbound) ! Массив значений периодограммы
     
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
     
     real(8) p_d, j_d ! Овеществления
     
     integer(4) p, t, j ! Вспомогательные переменные
     
     real(8) s1, s2    ! Временные держатели сумм для элементов выражений
     real(8) diff      ! Разность (A(j,2) - x_mean)
     real(8) p_cur     ! Текущее значение p
     real(8) koef      ! Аргумент тригонометрических функций
     real(8) cos_value ! Значение косинуса от аргумента koef
     real(8) sin_value ! Значение синуса от аргумента koef
     
     if (use_if .eq. 0) then
     
          open(11, file="output2")
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

                    koef = 2d0 * pi * p_cur * j_d / N_wif

                    cos_value = dcos(koef)
                    sin_value = dsin(koef)

                    ! Проверка на ошибку округления sin_value и cos_value
                    if (abs(cos_value) .le. 1e-3) cos_value = 0d0
                    if (abs(sin_value) .le. 1e-3) sin_value = 0d0

                    diff = (A(j,2) - x_mean)

                    s1 = s1 + diff * cos_value
                    s2 = s2 + diff * sin_value

               enddo
          
               ! Вычисление периодограммы по формуле 7.17 из
               ! Chatfield - The Analysis of Time Series, стр. 111
               ! с предварительным центрированием ряда (то есть с
               ! вычитанием среднего из выборки)

               I_p(p) = (s1 * s1 + s2 * s2) / ( N_wif * pi )
          
               if (abs(I_p(p)) .le. 1e-15) I_p(p) = 0d0

               write(11,'(e16.7, 1x, e16.7, 1x, e16.7)') N_d / p_cur, p_cur / N_d, I_p(p)

          enddo
          close(11)
     
     else
     
          open(11, file="output2")
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

                    koef = 2d0 * pi * p_cur * j_d / N_d

                    cos_value = dcos(koef)
                    sin_value = dsin(koef)

                    ! Проверка на ошибку округления sin_value и cos_value
                    if (abs(cos_value) .le. 1e-3) cos_value = 0d0
                    if (abs(sin_value) .le. 1e-3) sin_value = 0d0

                    diff = (A(t,2) - x_mean)

                    s1 = s1 + diff * cos_value
                    s2 = s2 + diff * sin_value

               enddo
          
               ! Вычисление периодограммы по формуле 7.17 из
               ! Chatfield - The Analysis of Time Series, стр. 111
               ! с предварительным центрированием ряда (то есть с
               ! вычитанием среднего из выборки)

               I_p(p) = (s1 * s1 + s2 * s2) / ( N_d * pi )
          
               if (abs(I_p(p)) .le. 1e-15) I_p(p) = 0d0

               write(11,'(e16.7, 1x, e16.7, 1x, e16.7)') N_d / p_cur, p_cur / N_d, I_p(p)

          enddo
          close(11)
     
     endif
     
     end subroutine
     
     
     ! [Вычисление коррелограммы через обратное преобразование Фурье к периодограмме]
     subroutine F4_correlogram_fourier_transform(I_p, leftbound, leftbound_d, rightbound, p_num, p_step, t_koef, N, N_d, &
     &pi, bias_fix)
     implicit none
     
     integer(4), intent(in) :: p_num                     ! Общее число множителей p
     integer(4), intent(in) :: N                         ! Размер выборки
     integer(4), intent(in) :: t_koef                    ! Множитель дискретизации множителей t (для периодов)
     integer(4), intent(in) :: leftbound, rightbound     ! Границы рабочего диапазона частот
     integer(4), intent(in) :: bias_fix                  ! Приводить коэффициенты корреляции к несмещённой оценке?
     real(8), intent(in)    :: I_p(leftbound:rightbound) ! Массив значений периодограммы
     real(8), intent(in)    :: p_step                    ! Шаг дискретизации множителей p (для частот)
     real(8), intent(in)    :: N_d                       ! Овеществление N
     real(8), intent(in)    :: leftbound_d               ! Овеществление leftbound
     real(8), intent(in)    :: pi                        ! Число pi
     
     real(8) p_d, t_d ! Овеществления
     
     integer(4) p, t   ! Вспомогательные переменные
     real(8) s1, s2    ! Временные держатели сумм для элементов выражений
     real(8) p_cur     ! Текущее значение p
     real(8) t_cur     ! Текущее значение t
     real(8) koef      ! Аргумент тригонометрических функций
     real(8) cos_value ! Значение косинуса от аргумента koef
     
     real(8) C(0:t_koef*(p_num - 10)) ! Вектор значений коррелограммы
     
     ! Вычисление коэффициента автокорреляции c(0)
     
     s2 = 0d0
     
     do p = leftbound, rightbound - 1, 1
               
               s2 = s2 + I_p(p)

     enddo
     
     ! Вычисление коэффициентов автокорреляции c(k)
     ! и деление их на значение коэффициента c(0)
     
     if (bias_fix .eq. 0) then
     
          open(12, file = 'output3')
          do t = 0, t_koef * (N - 1), 1
     
               t_d = t
               t_cur = 0d0 + t_d / t_koef
          
               s1 = 0d0
          
               do p = leftbound, rightbound - 1, 1

                    p_d = p
                
                    p_cur = leftbound_d - 1d0 + p_d * p_step
          
                    koef = 2d0 * pi * p_cur * t_cur / N_d

                    cos_value = dcos(koef)

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

               write(12,'(e16.7, 1x, e16.7)') t_cur, C(t)
     
          enddo     
          close(12)
     
     else
          
          open(12, file = 'output3')
          do t = 0, t_koef * (N - 1), 1
          
               t_d = t
               t_cur = 0d0 + t_d / t_koef
!               t_cur = t
               
               s1 = 0d0
               
               do p = leftbound, rightbound - 1, 1
     
                    p_d = p
                     
                    p_cur = leftbound_d - 1d0 + p_d * p_step
               
                    koef = 2d0 * pi * p_cur * t_cur / N_d
     
                    cos_value = dcos(koef)

                    ! Проверка на ошибку округления sin_value и cos_value
                    if (abs(cos_value) .le. 1e-3) cos_value = 0d0
               
                    s1 = s1 + I_p(p) * cos_value

               enddo
          
               ! Вычисление коррелограммы через обратное преобразование Фурье
               ! к периодограмме по формуле для коэффициентов автокорреляции
               ! 73 (вычисляя величину c(k)/c(0)) из Витязева - 
               ! Спектрально-корреляционный анализ равномерных рядов, стр. 25
          
               C(t) = s1 / s2

               if (abs(C(t)) .le. 1e-15) C(t) = 0d0

               write(12,'(e16.7, 1x, e16.7)') t_cur, C(t)
     
          enddo     
          close(12)
     
     endif
     
     end subroutine
     
end module
