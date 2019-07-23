module subprograms
use mpi_subprograms
implicit none

     character(*), parameter :: result_path = "./Результат/"

     contains

     ! [Формирование массива N_index_array и считывание исходных данных]
     subroutine F0_get_index_array(N_index_array, N_wif, N_if, N, A)
     implicit none
     
     integer(4), intent(in) :: N_wif ! Размер выборки с исключениями (N - N_if)
     integer(4), intent(in) :: N_if  ! Число исключений
     integer(4), intent(in) :: N     ! Размер выборки

     integer(4), intent(out) :: N_Index_array(0:N_wif-1) ! Массив индексов с учётом исключений
     integer(4), allocatable, dimension(:) :: N_if_array ! Массив индексов-исключений
     
     ! Вспомогательные переменные
     integer(4) i, k, ier
     
     ! Опциональная переменная
     real(8), optional, intent(in) :: A(0:N-1,2) ! Матрица исходных данных
     
     if (.not.(present(A))) then
     
          allocate(N_if_array(1:N_if + 1), stat = ier)
          if (ier .ne. 0) stop 'Не удалось выделить память для массива N_if_array'
     
          ! Заполнение массива индексов-исключений: первые N_if элементов заполняются
          ! индексами-исключениями, а последний элемент намеренно приравнивается значению 0
          
          N_if_array = (/ 2, (i, i = 4,8), 10, 54, 133, (i, i = 1541,1546), (i, i = 2141, 2149),&
          & 2425, 2426, (i, i = 2772, 2777), (i, i = 2807, 2809), (i, i = 2863, 2867), 2896, 2897,&
          & 3004, (i, i = 3117, 3123), (i, i = 3537, 3545), (i, i = 3551, 3556), (i, i = 3586, 3594),&
          & (i, i = 3602, 3607), (i, i = 3795, 3801), (i, i = 3810, 3953), (i, i = 3961, 4026),&
          & 5284, 5287, 5613, 5669, 0 /)
     
          k = 1 ! Уменьшаемое в сдвиге (k - 1) при обнаружении элемента из N_if_array;
                ! приводим таким образом массив индексов 1:N к размеру 1:N - N_if

          ! Заполнение массива N_index_array

          do i = 0, N - 1

               if (i .ne. N_if_array(k) - 1) then

                    N_index_array(i - k + 1) = i

               else

                    k = k + 1

               endif

          enddo
          
          deallocate(N_if_array)
     
     else
     
          k = 1 ! Уменьшаемое в сдвиге (k - 1) при обнаружении элемента из N_if_array;
                ! приводим таким образом массив индексов 1:N к размеру 1:N - N_if
     
          do i = 0, N - 1

               if (A(i,2) .ne. 0d0) then

                    N_index_array(i - k + 1) = i

               else

                    k = k + 1

               endif

          enddo
     
     endif
     
     end subroutine
     

     ! [Вычисление среднего значения выборки]
     subroutine F1_mean(A, x_mean, N_wif, N, N_d, N_index_array, mpiSize, mpiRank)
     implicit none
     
     integer(4), intent(in) :: N_wif                    ! Размер выборки с исключениями (N - N_if)
     integer(4), intent(in) :: N                        ! Размер выборки
     real(8), intent(in)    :: A(0:N-1,2)               ! Матрица исходных данных
     real(8), intent(in)    :: N_d                      ! Овеществление N
     real(8), intent(out)   :: x_mean                   ! Среднее значение выборки с учётом исключений
     
     ! Переменные для деления первого подпространства итераций на порции
     integer(4) t_size                    ! Общий размер пространства интераций
     integer(4) partion_size              ! Размер порции
     integer(4) partion_size_mod          ! Остаток от деления t_size на mpiSize
     integer(4) t_leftbound, t_rightbound ! Границы индексов для данного ранга
     integer(4) partion_shift             ! Величина mpiRank * partion_size
     
     ! Стандартные и вспомогательные переменные MPI
     integer(4), intent(in) :: mpiSize ! Размер коммуникатора
     integer(4), intent(in) :: mpiRank ! Ранг процесса
     integer(4) ierr   ! Переменная ошибки
     
     real(8) partion_x_mean ! Сумма A(2,j) в порции по массиву исключений
     
     ! Опциональный массив: использование зависит от ответа на вопрос о массиве исключений
     integer(4), optional, intent(in) :: N_Index_array(0:N_wif-1) ! Массив индексов с учётом исключений
     
     integer(4) t, j ! Вспомогательные переменные
     
     x_mean = 0d0
     partion_x_mean = 0d0
     
     if (present(N_index_array)) then
     
          ! Вычисление размеров порций и их границ [1]
     
          t_size = N_wif
     
          call partion_sizes(t_size, 0, mpiSize, mpiRank, partion_size, partion_shift, partion_size_mod, t_leftbound, t_rightbound)

          ! Выполнение процесcом своей порции [1]

          do t = t_leftbound, t_rightbound, 1

               j = N_index_array(t)

               partion_x_mean = partion_x_mean + A(j,2)

          enddo
     
          call mpi_allreduce(partion_x_mean, x_mean, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
          x_mean = x_mean / N_wif
     
     else
     
          ! Вычисление размеров порций и их границ [2]
     
          t_size = N
     
          call partion_sizes(t_size, 0, mpiSize, mpiRank, partion_size, partion_shift, partion_size_mod, t_leftbound, t_rightbound)
     
          do t = t_leftbound, t_rightbound, 1

               partion_x_mean = partion_x_mean + A(t,2)

          enddo
     
          call mpi_allreduce(partion_x_mean, x_mean, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
          x_mean = x_mean / N_d
     
     endif

     end subroutine
     
     
     ! [Вычисление коррелограммы прямым способом]
     subroutine F2_correlogram_direct(A, x_mean, N, N_d, bias_fix, mpiSize, mpiRank)
     implicit none
     
     integer(4), intent(in) :: N          ! Размер выборки
     integer(4), intent(in) :: bias_fix   ! Приводить коэффициенты корреляции к несмещённой оценке?
     real(8), intent(in)    :: A(0:N-1,2) ! Матрица исходных данных
     real(8), intent(in)    :: x_mean     ! Среднее значение выборки с учётом исключений
     real(8), intent(in)    :: N_d        ! Овеществление N
     
     ! Переменные для деления первого подпространства итераций на порции
     integer(4) k_size                    ! Общий размер пространства интераций
     integer(4) partion_size              ! Размер порции
     integer(4) partion_size_mod          ! Остаток от деления k_size на mpiSize
     integer(4) k_leftbound, k_rightbound ! Границы индексов для данного ранга
     integer(4) partion_shift             ! Величина mpiRank * partion_size
     
     ! Вспомогательные переменные при обмене сообщениями
     integer(4) send_leftbound, send_rightbound ! Границы индексов для данного ранга при ранге i
     
     ! Стандартные и вспомогательные переменные MPI
     integer(4), intent(in) :: mpiSize  ! Размер коммуникатора
     integer(4), intent(in) :: mpiRank  ! Ранг процесса
     integer(4) status(MPI_STATUS_SIZE) ! Переменная статуса передачи
     integer(4) ierr   ! Переменная ошибки
     
     integer(4) t, k, i         ! Вспомогательные переменные
     real(8) s1, s2, partion_s2 ! Временные держатели сумм для элементов выражений
     real(8) diff               ! Разность (A(t,2) - x_mean)
     
     real(8) r(0:N - 1)    ! Вектор значений коррелограммы
     real(8) k_d(0:N - 1)  ! Вектор значений лага (k)
     
     ! Вычисление размеров порций и их границ
     
     k_size = N
     
     call partion_sizes(k_size, 0, mpiSize, mpiRank, partion_size, partion_shift, partion_size_mod, k_leftbound, k_rightbound)
     
     ! Вычисление коэффициента автокорреляции c(0)
     
     ! Выполнение процесcом своей порции [1]
     
     s2 = 0d0
     partion_s2 = 0d0
     
     do t = k_leftbound, k_rightbound

          diff = (A(t,2) - x_mean)
          partion_s2 = partion_s2 + diff * diff

     enddo
     
     ! Суммирование элементов массива diff и передача результата всем процессам
     call mpi_allreduce(partion_s2, s2, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
     
     ! Вычисление коэффициентов автокорреляции c(k) по 
     ! формуле 4.1 из Chatfield - The Analysis of 
     ! Time Series, стр. 49, и деление их на значение
     ! коэффициента c(0) согласно формуле 4.4
     
     ! Нормировка для несмещённой оценки взята из
     ! Витязева - Спектрально-корреляционный анализ 
     ! равномерных рядов
     
     ! Выполнение процесcом своей порции [2]
     
     if (bias_fix .eq. 0) then ! Приводить коэффициенты корреляции к несмещённой оценке?
     
          do k = k_leftbound, k_rightbound, 1

               k_d(k) = k

               s1 = 0d0

               do t = 0, N - k - 1, 1

                    s1 = s1 + (A(t,2) - x_mean) * (A(t+k,2) - x_mean)

               enddo

               r(k) = s1 * N_d / (N_d - k_d(k)) / s2

          enddo
     
     else
     
          do k = k_leftbound, k_rightbound, 1

               k_d(k) = k

               s1 = 0d0

               do t = 0, N - k - 1, 1

                    s1 = s1 + (A(t,2) - x_mean) * (A(t+k,2) - x_mean)

               enddo

               r(k) = s1 / s2

          enddo
     
     endif
     
     ! Передача всех порций процессу 0
     
     if (mpiRank .gt. 0) then
     
          call mpi_send(k_d(k_leftbound:k_rightbound), partion_size, MPI_REAL8, 0, mpiRank, MPI_COMM_WORLD, ierr)
          call mpi_send(r(k_leftbound:k_rightbound), partion_size, MPI_REAL8, 0, mpiRank, MPI_COMM_WORLD, ierr)
          
     else
     
          do i = 1, mpiSize - 1
               
               send_leftbound = i * partion_size
               
               if (i .eq. mpiSize - 1 .and. partion_size_mod .ne. 0) then
               
                    send_rightbound = (partion_size + partion_size_mod) + i * partion_size - 1
                    
                    call mpi_recv(k_d(send_leftbound:send_rightbound), (partion_size + partion_size_mod), MPI_REAL8, i, i, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(r(send_leftbound:send_rightbound), (partion_size + partion_size_mod), MPI_REAL8, i, i, MPI_COMM_WORLD, status, ierr)
               
               else
               
                    send_rightbound = partion_size + i * partion_size - 1
                    
                    call mpi_recv(k_d(send_leftbound:send_rightbound), partion_size, MPI_REAL8, i, i, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(r(send_leftbound:send_rightbound), partion_size, MPI_REAL8, i, i, MPI_COMM_WORLD, status, ierr)
               
               endif
               
          enddo
          
     endif
     
     ! Сообщение результатов другим процессам
     
     call mpi_bcast(k_d, k_size, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
     call mpi_bcast(r, k_size, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
     
     ! Вывод результата первым процессом
     
     if (mpiRank .eq. 0) then
     open(17, file = result_path//"direct_correlogram.dat")
     
     write(17, '(a, /, a, 5x, a, /, a)') '#', '# Time Lag, Days, k', 'Autocorrelation Coefficient, r(k)', '#'     
     
     do k = 0, N - 1, 1
      
          write(17, '(f11.1, 17x, e23.15)') k_d(k), r(k)
     
     enddo
     
     close(17)
     endif
     
     end subroutine


     ! [Вычисление периодограммы]
     subroutine F3_periodogram(A, leftbound, leftbound_d, rightbound, p_step, N, N_d, N_wif, x_mean, pi, I_p, N_index_array, mpiSize, mpiRank)
     implicit none
     
     integer(4), intent(in) :: leftbound, rightbound    ! Границы рабочего диапазона частот
     integer(4), intent(in) :: N                        ! Размер выборки
     integer(4), intent(in) :: N_wif                    ! Размер выборки с исключениями (N - N_if)
     
     real(8), intent(in)    :: A(0:N-1,2)  ! Матрица исходных данных
     real(8), intent(in)    :: p_step      ! Шаг дискретизации множителей p (для частот)
     real(8), intent(in)    :: x_mean      ! Среднее значение выборки с учётом исключений
     real(8), intent(in)    :: pi          ! Число pi
     real(8), intent(in)    :: N_d         ! Овеществление N
     real(8), intent(in)    :: leftbound_d ! Овеществление leftbound
     
     ! Переменные для деления первого подпространства итераций на порции
     integer(4) p_size                    ! Общий размер пространства интераций
     integer(4) partion_size              ! Размер порции
     integer(4) partion_size_mod          ! Остаток от деления p_size на mpiSize
     integer(4) p_leftbound, p_rightbound ! Границы индексов для данного ранга
     integer(4) partion_shift             ! Величина mpiRank * partion_size
     
     ! Вспомогательные переменные при обмене сообщениями
     integer(4) send_leftbound, send_rightbound ! Границы индексов для данного ранга при ранге i
     
     ! Стандартные и вспомогательные переменные MPI
     integer(4), intent(in) :: mpiSize  ! Размер коммуникатора
     integer(4), intent(in) :: mpiRank  ! Ранг процесса
     integer(4) status(MPI_STATUS_SIZE) ! Переменная статуса передачи
     integer(4) ierr   ! Переменная ошибки
     
     real(8), intent(out) :: I_p(leftbound:rightbound) ! Массив значений периодограммы
     real(8) period(leftbound:rightbound)              ! Вектор значений периодов
     real(8) frequency(leftbound:rightbound)           ! Вектор значений частот
     
     ! Опциональный массив: использование зависит от ответа на вопрос о массиве исключений
     integer(4), optional, intent(in) :: N_Index_array(0:N_wif-1) ! Массив индексов с учётом исключений     

     real(8) p_d, j_d ! Овеществления
     
     integer(4) p, t, j, i ! Вспомогательные переменные
     
     real(8) s1, s2    ! Временные держатели сумм для элементов выражений
     real(8) diff      ! Разность (A(j,2) - x_mean)
     real(8) p_cur     ! Текущее значение p
     real(8) arg       ! Аргумент тригонометрических функций
     real(8) cos_value ! Значение косинуса от аргумента arg
     real(8) sin_value ! Значение синуса от аргумента arg
     
     ! Вычисление размеров порций и их границ
     
     p_size = rightbound - leftbound + 1
     
     call partion_sizes(p_size, leftbound, mpiSize, mpiRank, partion_size, partion_shift, partion_size_mod, p_leftbound, p_rightbound)
     
     ! Выполнение процесcом своей порции     
     
     if (present(N_index_array)) then ! Использовать массив исключений?
     
          do p = p_leftbound, p_rightbound, 1

               p_d = p
               p_cur = leftbound_d - 1d0 + p_d * p_step
          
               s1 = 0d0
               s2 = 0d0

               do t = 0, N_wif - 1, 1

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
               
               period(p) = N_d / p_cur
               frequency(p) = p_cur / N_d
          
          enddo
     
     else
     
          do p = p_leftbound, p_rightbound, 1

               p_d = p
               p_cur = leftbound_d - 1d0 + p_d * p_step
          
               s1 = 0d0
               s2 = 0d0

               do t = 0, N - 1, 1
     
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
               
               period(p) = N_d / p_cur
               frequency(p) = p_cur / N_d
          
          enddo
     
     endif
     
     ! Передача всех порций процессу 0
     
     if (mpiRank .gt. 0) then
     
          call mpi_send(period(p_leftbound:p_rightbound), partion_size, MPI_REAL8, 0, mpiRank, MPI_COMM_WORLD, ierr)
          call mpi_send(frequency(p_leftbound:p_rightbound), partion_size, MPI_REAL8, 0, mpiRank, MPI_COMM_WORLD, ierr)
          call mpi_send(I_p(p_leftbound:p_rightbound), partion_size, MPI_REAL8, 0, mpiRank, MPI_COMM_WORLD, ierr)
          
     else
     
          do i = 1, mpiSize - 1
               
               send_leftbound = leftbound + i * partion_size
               
               if (i .eq. mpiSize - 1 .and. partion_size_mod .ne. 0) then
               
                    send_rightbound = leftbound + (partion_size + partion_size_mod) + i * partion_size - 1
                    
                    call mpi_recv(period(send_leftbound:send_rightbound), (partion_size + partion_size_mod), MPI_REAL8, i, i, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(frequency(send_leftbound:send_rightbound), (partion_size + partion_size_mod), MPI_REAL8, i, i, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(I_p(send_leftbound:send_rightbound), (partion_size + partion_size_mod), MPI_REAL8, i, i, MPI_COMM_WORLD, status, ierr)
               
               else
               
                    send_rightbound = leftbound + partion_size + i * partion_size - 1
                    
                    call mpi_recv(period(send_leftbound:send_rightbound), partion_size, MPI_REAL8, i, i, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(frequency(send_leftbound:send_rightbound), partion_size, MPI_REAL8, i, i, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(I_p(send_leftbound:send_rightbound), partion_size, MPI_REAL8, i, i, MPI_COMM_WORLD, status, ierr)
               
               endif
               
          enddo
          
     endif
     
     ! Сообщение результатов другим процессам
     
     call mpi_bcast(period, p_size, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
     call mpi_bcast(frequency, p_size, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
     call mpi_bcast(I_p, p_size, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
     
     ! Вывод результата первым процессом
     
     if (mpiRank .eq. 0) then
     open(18, file = result_path//"periodogram.dat")
     
     write(18, '(a, /, a, 6x, a, 6x, a, /, a)') '#', '#    Period, Days, T', 'Frequency, 1/Days, v', 'Periodogram, I(v)', '#'
     
     do p = leftbound, rightbound, 1
      
          write(18, '(e23.15, 1x, e23.15, 1x, e23.15)') period(p), frequency(p), I_p(p)
     
     enddo
     
     close(18)
     endif
     
     end subroutine
     
     
     ! [Вычисление коррелограммы через применение обратного преобразования Фурье к периодограмме]
     subroutine F4_correlogram_fourier_transform(I_p, leftbound, leftbound_d, rightbound, p_step, t_koef, N, N_d, pi, bias_fix, mpiSize, mpiRank)
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
     
     ! Переменные для деления первого подпространства итераций на порции
     integer(4) t_size                    ! Общий размер пространства интераций
     integer(4) partion_size              ! Размер порции
     integer(4) partion_size_mod          ! Остаток от деления t_size на mpiSize
     integer(4) t_leftbound, t_rightbound ! Границы индексов для данного ранга
     integer(4) partion_shift             ! Величина mpiRank * partion_size
     
     ! Вспомогательные переменные при обмене сообщениями
     integer(4) send_leftbound, send_rightbound ! Границы индексов для данного ранга при ранге i
     
     ! Стандартные и вспомогательные переменные MPI
     integer(4), intent(in) :: mpiSize  ! Размер коммуникатора
     integer(4), intent(in) :: mpiRank  ! Ранг процесса
     integer(4) status(MPI_STATUS_SIZE) ! Переменная статуса передачи
     integer(4) ierr   ! Переменная ошибки
     
     real(8) p_d, t_d, t_koef_d ! Овеществления
     
     integer(4) p, t, i ! Вспомогательные переменные
     real(8) s1, s2     ! Временные держатели сумм для элементов выражений
     real(8) p_cur      ! Текущее значение p
     real(8) arg        ! Аргумент тригонометрических функций
     real(8) cos_value  ! Значение косинуса от аргумента arg
     
     real(8) t_cur(0:t_koef * (N - 1)) ! Вектор текущих значений t
     real(8) C(0:t_koef * (N - 1))     ! Вектор значений коррелограммы
     
     ! Вычисление размеров порций и их границ
     
     t_size = t_koef * (N - 1) + 1
     
     call partion_sizes(t_size, 0, mpiSize, mpiRank, partion_size, partion_shift, partion_size_mod, t_leftbound, t_rightbound)
     
     ! Вычисление коэффициента автокорреляции c(0)
     
     s2 = 0d0
     s2 = sum(I_p)
     
     ! Вычисление коэффициентов автокорреляции c(k)
     ! и деление их на значение коэффициента c(0)
     
     ! Выполнение процесcом своей порции
     
     t_koef_d = t_koef
     
     if (bias_fix .eq. 0) then ! Приводить коэффициенты корреляции к несмещённой оценке?
     
          do t = t_leftbound, t_rightbound, 1
     
               t_d = t
               t_cur(t) = 0d0 + t_d / t_koef_d
               
               s1 = 0d0
          
               do p = leftbound, rightbound, 1

                    p_d = p
                    p_cur = leftbound_d - 1d0 + p_d * p_step
          
                    arg = 2d0 * pi * p_cur * t_cur(t) / N_d

                    cos_value = dcos(arg)

                    ! Проверка на ошибку округления sin_value и cos_value
                    if (abs(cos_value) .le. 1e-3) cos_value = 0d0
               
                    s1 = s1 + I_p(p) * cos_value

               enddo
          
               ! Вычисление коррелограммы через обратное преобразование Фурье
               ! к периодограмме по формуле для коэффициентов автокорреляции
               ! 73 (вычисляя величину c(k)/c(0)) из Витязева - 
               ! Спектрально-корреляционный анализ равномерных рядов, стр. 25
          
               C(t) = s1 * N_d / (N_d - t_cur(t)) / s2

               if (abs(C(t)) .le. 1e-15) C(t) = 0d0

          enddo     
     
     else
          
          do t = t_leftbound, t_rightbound, 1
          
               t_d = t
               t_cur(t) = 0d0 + t_d / t_koef_d
               
               s1 = 0d0
               
               do p = leftbound, rightbound, 1
     
                    p_d = p
                    p_cur = leftbound_d - 1d0 + p_d * p_step
               
                    arg = 2d0 * pi * p_cur * t_cur(t) / N_d
     
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

          enddo  
          
     endif   
     
     ! Передача всех порций процессу 0
     
     if (mpiRank .gt. 0) then
     
          call mpi_send(t_cur(t_leftbound:t_rightbound), partion_size, MPI_REAL8, 0, mpiRank, MPI_COMM_WORLD, ierr)
          call mpi_send(C(t_leftbound:t_rightbound), partion_size, MPI_REAL8, 0, mpiRank, MPI_COMM_WORLD, ierr)
          
     else
     
          do i = 1, mpiSize - 1
               
               send_leftbound = i * partion_size
               
               if (i .eq. mpiSize - 1 .and. partion_size_mod .ne. 0) then
               
                    send_rightbound = (partion_size + partion_size_mod) + i * partion_size - 1
                    
                    call mpi_recv(t_cur(send_leftbound:send_rightbound), (partion_size + partion_size_mod), MPI_REAL8, i, i, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(C(send_leftbound:send_rightbound), (partion_size + partion_size_mod), MPI_REAL8, i, i, MPI_COMM_WORLD, status, ierr)
               
               else
               
                    send_rightbound = partion_size + i * partion_size - 1
                    
                    call mpi_recv(t_cur(send_leftbound:send_rightbound), partion_size, MPI_REAL8, i, i, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(C(send_leftbound:send_rightbound), partion_size, MPI_REAL8, i, i, MPI_COMM_WORLD, status, ierr)
               
               endif
               
          enddo
          
     endif
     
     ! Вывод результата первым процессом
     
     if (mpiRank .eq. 0) then
     open(19, file = result_path//'reverse_correlogram.dat')
     
     write(19, '(a, /, a, 5x, a, /, a)') '#', '# Time Lag, Days, k', 'Autocorrelation Coefficient, r(k)', '#'
     
     do t = 0, t_koef * (N - 1), 1
      
          write(19, '(f11.1, 17x, e23.15)') t_cur(t), C(t)
     
     enddo
     
     close(19)
     endif
     
     end subroutine
     
end module
