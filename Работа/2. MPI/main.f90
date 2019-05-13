program main
use mpi
implicit none

     real(8), allocatable, dimension(:,:) :: A ! Матрица исходных данных
     
     integer N     ! Размер выборки
     integer N_if  ! Число исключений
     integer N_wif ! Размер выборки с исключениями (N - N_if)                      
     
     ! Массив индексов-исключений
     integer(4), allocatable, dimension(:) :: N_if_array
     
     ! Массив индексов без индексов исключений
     integer(4), allocatable, dimension(:) :: N_index_array
     
     real(8) x_mean    ! Среднее значение x
     real(8) koef_mean ! Коэффициент 1/N
     
     ! Коррелограмма:

     real(8), allocatable, dimension(:) :: r     ! Вектор коэффициентов корреляции
     real(8), allocatable, dimension(:) :: k_arr ! Вектор значений множителей k

     ! Периодограмма:

     real(8), allocatable, dimension(:) :: I_p   ! Вектор значений периодограммы
     real(8), allocatable, dimension(:) :: p_arr ! Вектор значений множителей p
     
     real(8) p_step ! Шаг дискретизации частот
     integer p_num  ! Число желаемых данных
     real(8) p_cur  ! Текущее значение p

     real(8) koef                 ! Коэффициент для тригонометрических функций
     real(8) cos_value, sin_value ! Вычисленные значения косинуса и синуса
     real(8) diff                 ! Разность (A(t,2) - x_mean)

     ! Общие переменные

     real(8) s1, s2 ! Временные суммы для элементов выражений
     real(8) pi ! Число pi
     
     ! Вспомогательные переменные
     integer k, t, i, j, p, ier
     real(8) k_d, N_d, p_d, j_d, p_num_d ! Овеществления
     
     ! Выбор рабочего диапазона частот
     logical full_range ! Использовать полный или указанный рабочие диапазоны? (.true., если полный)
     integer leftbound ! Левая граница рабочего диапазона частот
     integer rightbound ! Праввая граница рабочего диапазона частот
     
     ! Переменные для деления первого подпространства итераций на порции
     integer(4) p_size                    ! Длина рабочего диапазона частот
     integer(4) p_partion_size            ! Размер порции
     integer(4) p_partion_size_mod        ! Остаток от деления p_size на mpiSize
     integer(4) p_leftbound, p_rightbound ! Границы индексов для данного ранга
     integer(4) p_partion_shift           ! Величина mpiRank * p_partion_size
     
     integer(4) k_partion_size            ! Размер порции
     integer(4) k_partion_size_mod        ! Остаток от деления k_size на mpiSize
     integer(4) k_leftbound, k_rightbound ! Границы индексов для данного ранга
     integer(4) k_partion_shift           ! Величина mpiRank * k_partion_size
     
     ! Вспомогательные переменные при обмене сообщениями
     integer(4) send_leftbound, send_rightbound ! Границы индексов для данного ранга при ранге i  

     ! Вспомогательные переменные MPI
     integer(4) mpiErr, mpiSize, mpiRank
     integer(4) ierr, status

     ! Указать размер выборки
     N = 5860
     
     ! Указать число исключений
     N_if = 301
     
     ! Использовать указанный диапазон частот для вычисления периодограммы 
     ! или считать по полному диапазону?
     full_range = .true.
     
     ! Настройки рабочего диапазона частот в разделе периодограмма
     
     ! Вычисление размера выборки с исключениями
     N_wif = N - N_if
     
     allocate(N_if_array(N_if), stat = ier)
     if (ier .ne. 0) stop 'Не удалось выделить память для массива N_if_array'
     
     allocate(N_index_array(N_wif), stat = ier)
     if (ier .ne. 0) stop 'Не удалось выделить память для массива N_index_array'
     
     N_if_array = (/ 2, (i, i = 4,8), 10, 54, 133, (i, i = 1541,1546), (i, i = 2141, 2149),&
     & 2425, 2426, (i, i = 2772, 2777), (i, i = 2807, 2809), (i, i = 2863, 2867), 2896, 2897,&
     & 3004, (i, i = 3117, 3123), (i, i = 3537, 3545), (i, i = 3551, 3556), (i, i = 3586, 3594),&
     & (i, i = 3602, 3607), (i, i = 3795, 3801), (i, i = 3810, 3953), (i, i = 3961, 4026),&
     & 5284, 5287, 5613, 5669 /)
     
     k = 1 ! Сдвиг при обнаружении элемента из N_if_array,
           ! уменьшаем таким образом массив индексов 1:N до размера N - N_if
     
     ! Заполнение массива N_index_array
     
     do i = 1, N
     
          if (i .ne. N_if_array(k)) then
               
               N_index_array(i - k + 1) = i
               
               else
               
               k = k + 1; cycle
               
          endif
          
     enddo
     
     ! Выделение памяти под рабочие массивы
     
     ! Массив исходных данных
     allocate(A(1:N,2), stat = ier)
     if (ier .ne. 0) stop 'Не удалось выделить память для массива A'

     ! Считывание исходных данных
     
     call mpi_init(mpiErr)

     call mpi_comm_size(MPI_COMM_WORLD, mpiSize, mpiErr)
     call mpi_comm_rank(MPI_COMM_WORLD, mpiRank, mpiErr)

     if (mpiRank .eq. 0) then
     
     do i = 1, N
        
                 read(*,*) A(i,1), A(i,2)

     enddo
     
     endif
     
     call mpi_bcast(A, N*2, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

     ! Вычисление среднего значения выборки

     x_mean = 0d0

     N_d = N
     koef_mean = 1d0/(N_d - N_if)

     do t = 1, N_wif
     
          j = N_index_array(t)
     
          x_mean = x_mean + A(j,2)
     
     enddo

     x_mean = koef_mean * x_mean

     ! [Вычисление коррелограммы прямым методом]
     
     ! Массив коэффициентов автокорреляции
     allocate(r(1:N-1), stat = ier)
     if (ier .ne. 0) stop 'Не удалось выделить память для массива r'
     
     allocate(k_arr(1:N-1), stat = ier)
     if (ier .ne. 0) stop 'Не удалось выделить память для массива k_arr'
     
     ! Коррелограмма: вычисление размеров порции
     
     k_partion_size_mod = mod(N-1,mpiSize)
     
     if (k_partion_size_mod .eq. 0) then
               
          k_partion_size = (N - 1) / mpiSize
          k_partion_shift = mpiRank * k_partion_size
               
     elseif (mpiRank .eq. mpiSize - 1) then
               
          k_partion_size = ((N - 1) + (mpiSize - k_partion_size_mod)) / mpiSize - (mpiSize - k_partion_size_mod)
          k_partion_shift = mpiRank * (k_partion_size + (mpiSize - k_partion_size_mod))
          
     else
     
          k_partion_size = ((N - 1) + (mpiSize - k_partion_size_mod)) / mpiSize
          k_partion_shift = mpiRank * k_partion_size
          
     endif
     
     ! Коррелограмма: вычисление порции
     
     k_leftbound = 1 + k_partion_shift
     k_rightbound = k_partion_size + k_partion_shift
     
     do k = k_leftbound, k_rightbound
 
          k_d = k
          
          k_arr(k) = k_d
                
          s1 = 0d0
  
          do t = 1, N_wif - k
          
               j = N_index_array(t)
          
               s1 = s1 + (A(j,2) - x_mean) * (A(j+k,2) - x_mean)
          
          enddo

          s2 = 0d0

          do t = 1, N_wif
          
               j = N_index_array(t)
               
               s2 = s2 + (A(j,2) - x_mean) * (A(j,2) - x_mean)
          
          enddo

          r(k) = s1 / s2

     enddo
     
     ! Коррелограмма: передача всех порций процессу 0
     
     if (mpiRank .gt. 0) then
     
          call mpi_send(r(k_leftbound:k_rightbound), k_partion_size, MPI_REAL8, 0, mpiRank, MPI_COMM_WORLD, ierr)
          call mpi_send(k_arr(k_leftbound:k_rightbound), k_partion_size, MPI_REAL8, 0, mpiRank, MPI_COMM_WORLD, ierr)
          
     else
     
          do i = 1, mpiSize - 1
               
               send_leftbound = 1 + i * k_partion_size
               
               if (i .eq. mpiSize - 1 .and. k_partion_size_mod .ne. 0) then
               
                    send_rightbound = k_partion_size + i * k_partion_size - (mpiSize - k_partion_size_mod)
               
               else
               
                    send_rightbound = k_partion_size + i * k_partion_size
               
               endif
               
               call mpi_recv(r(send_leftbound:send_rightbound), k_partion_size, MPI_REAL8, i, i, MPI_COMM_WORLD, status, ierr)
               call mpi_recv(k_arr(send_leftbound:send_rightbound), k_partion_size, MPI_REAL8, i, i, MPI_COMM_WORLD, status, ierr)

          enddo
          
     endif
     
     ! Коррелограмма: вывод результата в файл
     
     if (mpiRank .eq. 0) then
     
          open(10, file="output1")
          
               write(10,'(e28.20, 1x, e28.20)') (k_arr(i), r(i), i = 1, N - 1)
          
          close(10)
     
     endif

     ! Определение pi
     pi = 4d0*datan(1d0)

     ! [Вычисление периодограммы вне зависимости от коэффициентов автокорреляции]

     p_num = 586000
     p_num_d = p_num
     p_step = N_d / p_num_d
     
     ! Указание рабочего диапазона частот (в первом условии)
     
     if (.not. full_range) then
     
          leftbound = 1
          rightbound = 301
     
     else
     
          leftbound = 1
          rightbound = p_num
     
     endif
     
     ! Массив значений периодограммы
     allocate(I_p(leftbound:rightbound), stat = ier)
     if (ier .ne. 0) stop 'Не удалось выделить память для массива I_p'

     ! Массив значений множителей p
     allocate(p_arr(leftbound:rightbound), stat = ier)
     if (ier .ne. 0) stop 'Не удалось выделить память для массива p_arr'

     ! Периодограмма: вычисление размера порции
     
     p_size = rightbound - leftbound + 1
     
     p_partion_size_mod = mod(p_size,mpiSize)
     
     if (p_partion_size_mod .eq. 0) then
               
          p_partion_size = p_size / mpiSize
          p_partion_shift = mpiRank * p_partion_size
               
     elseif (mpiRank .eq. mpiSize - 1) then
               
          p_partion_size = (p_size + (mpiSize - p_partion_size_mod)) / mpiSize - (mpiSize - p_partion_size_mod)
          p_partion_shift = mpiRank * (p_partion_size + (mpiSize - p_partion_size_mod))
          
     else
     
          p_partion_size = (p_size + (mpiSize - p_partion_size_mod)) / mpiSize
          p_partion_shift = mpiRank * p_partion_size
          
     endif
     
     ! Периодограмма: заполнение процесcом своей порции
     
     p_leftbound = leftbound + p_partion_shift
     p_rightbound = leftbound + p_partion_size + p_partion_shift - 1
     
     do p = p_leftbound, p_rightbound, 1
     
          p_d = p
          p_cur = 0d0 + p_d * p_step
     
          p_arr(p) = p_cur

          s1 = 0d0
          s2 = 0d0

          do t = 1, N_wif

               j = N_index_array(t)
     
               j_d = j
                             
               koef = 2d0 * pi * p_cur * j_d / (N_d - N_if)
!               koef = 2d0 * pi * p_cur * j_d / N_d
                             
               cos_value = dcos(koef)
               sin_value = dsin(koef)
                             
               ! Проверка на ошибку округления sin_value и cos_value
               if (abs(cos_value) .le. 1e-3) cos_value = 0d0
               if (abs(sin_value) .le. 1e-3) sin_value = 0d0       
                                                             
               diff = (A(j,2) - x_mean)
                             
               s1 = s1 + diff * cos_value
               s2 = s2 + diff * sin_value
                             
          enddo

          I_p(p) = (s1 * s1 + s2 * s2)/( (N_d - N_if) * pi)
!          I_p(p) = (s1 * s1 + s2 * s2)/(N_d * pi)

     enddo

     ! Периодограмма: передача всех порций процессу 0
     
     if (mpiRank .gt. 0) then
     
          call mpi_send(I_p(p_leftbound:p_rightbound), p_partion_size, MPI_REAL8, 0, mpiRank, MPI_COMM_WORLD, ierr)
          call mpi_send(p_arr(p_leftbound:p_rightbound), p_partion_size, MPI_REAL8, 0, mpiRank, MPI_COMM_WORLD, ierr)
          
     else
     
          do i = 1, mpiSize - 1
               
               send_leftbound = leftbound + i * p_partion_size
               
               if (i .eq. mpiSize - 1 .and. p_partion_size_mod .ne. 0) then
               
                    send_rightbound = leftbound + p_partion_size + i * p_partion_size - 1 - (mpiSize - p_partion_size_mod)
               
               else
               
                    send_rightbound = leftbound + p_partion_size + i * p_partion_size - 1
               
               endif
               
               call mpi_recv(I_p(send_leftbound:send_rightbound), p_partion_size, MPI_REAL8, i, i, MPI_COMM_WORLD, status, ierr)
               call mpi_recv(p_arr(send_leftbound:send_rightbound), p_partion_size, MPI_REAL8, i, i, MPI_COMM_WORLD, status, ierr)     

          enddo
          
     endif
     
     ! Периодограмма: вывод результата в файл
     
     if (mpiRank .eq. 0) then
     
          open(11, file="output2")
          
               write(11,'(e16.7, 1x, e16.7, 1x, e16.7)') (N_d/p_arr(i), p_arr(i)/N_d, I_p(i), i = leftbound, rightbound)
          
          close(11)
     
     endif
     
     call mpi_finalize(mpiErr)

     deallocate(A, r, I_p, N_if_array, N_index_array, p_arr, k_arr)

end
