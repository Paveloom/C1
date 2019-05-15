program main
implicit none

     real(8), allocatable, dimension(:,:) :: A ! Матрица исходных данных

     integer N     ! Размер выборки
     integer N_if  ! Число исключений
     integer N_wif ! Размер выборки с исключениями (N - N_if)

     ! Массив индексов-исключений
     integer(4), allocatable, dimension(:) :: N_if_array

     ! Массив индексов без индексов исключений
     integer(4), allocatable, dimension(:) :: N_index_array

     ! Коррелограмма [прямой метод]:

     real(8), allocatable, dimension(:) :: r ! Вектор коэффициентов корреляции
     real(8) x_mean                          ! Среднее значение x
     real(8) koef_mean                       ! Коэффициент 1/N

     ! Периодограмма:

     real(8), allocatable, dimension(:) :: I_p ! Вектор значений периодограммы
     real(8), allocatable, dimension(:) :: nu  ! Вектор частот
     real(8) p_step                            ! Шаг дискретизации частот
     integer p_num                             ! Общее число множителей p
     real(8) p_cur                             ! Текущее значение p

     real(8) koef                 ! Коэффициент для тригонометрических функций
     real(8) cos_value, sin_value ! Вычисленные значения косинуса и синуса
     real(8) diff                 ! Разность (A(t,2) - x_mean)
     
     ! Коррелограмма [обратное преобразование Фурье к периодограмме]:
     
     real(8), allocatable, dimension(:) :: C ! Вектор значений корреляционной функции
     
     integer(4) cut_p_num ! Число множителей p с учетом вырезания высоких частот (указывается)
     real(8)    cut_N     ! Число множителей p с учетом вырезания высоких частот и смещением
     real(8)    t_cur     ! Текущее значение t

     ! Общие переменные

     real(8) s1, s2 ! Временные суммы для элементов выражений
     real(8) pi ! Число pi

     ! Вспомогательные переменные
     integer k, t, i, j, p, ier
     
     ! Овеществления
     real(8) k_d, N_d, p_d, j_d, t_d, p_num_d, cut_p_num_d, leftbound_d, rightbound_d

     ! Выбор рабочего диапазона частот
     logical full_range ! Использовать полный или указанный рабочие диапазоны? (.true., если полный)
     integer leftbound  ! Левая граница рабочего диапазона частот
     integer rightbound ! Праввая граница рабочего диапазона частот

     ! Указать размер выборки
     N = 5860

     ! Указать число исключений
     N_if = 301

     ! Использовать указанный диапазон частот для вычисления периодограммы
     ! или считать по полному диапазону?
     full_range = .false.

     ! Настройки рабочего диапазона частот в разделе периодограмма

     ! Вычисление размера выборки с исключениями
     N_wif = N - N_if

     allocate(N_if_array(N_if), stat = ier)
     if (ier .ne. 0) stop 'Не удалось выделить память для массива N_if_array'

     allocate(N_index_array(0:N_wif-1), stat = ier)
     if (ier .ne. 0) stop 'Не удалось выделить память для массива N_index_array'

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
     
     ! Выделение памяти под рабочие массивы

     ! Исходные данные
     allocate(A(0:N-1,2), stat = ier)
     if (ier .ne. 0) stop 'Не удалось выделить память для массива A'

     ! Для вычисления коррелограммы
     allocate(r(0:N-1), stat = ier)
     if (ier .ne. 0) stop 'Не удалось выделить память для массива r'

     ! Считывание исходных данных

        do i = 0, N - 1

                read(*,*) A(i,1), A(i,2)

        enddo
        
     ! Вычисление среднего значения выборки

     x_mean = 0d0

     N_d = N
     koef_mean = 1d0 / N_wif

     do t = 0, N_wif - 1

          j = N_index_array(t)

          x_mean = x_mean + A(j,2)

     enddo

     x_mean = koef_mean * x_mean

     ! Вычисление коррелограммы
     
     s2 = 0d0
     
     do t = 0, N - 1

!          j = N_index_array(t)
          
          diff = (A(t,2) - x_mean)

          s2 = s2 + diff * diff

     enddo
     
     open(10, file="output1")
     do k = 0, N - 1

          k_d = k

          s1 = 0d0

          do t = 0, N - k - 1

!               j = N_index_array(t)

               s1 = s1 + (A(t,2) - x_mean) * (A(t+k,2) - x_mean)

          enddo

          r(k) = s1 / s2

          write(10,'(e16.7, 1x, e16.7)') k_d, r(k)

     enddo
     close(10)

     ! Определение pi
     pi = 4d0 * datan(1d0)

     ! Вычисление периодограммы вне зависимости от коэффициентов автокорреляции

     p_num = 58600
     p_num_d = p_num
     p_step = N_d / p_num_d

     ! Указание рабочего диапазона частот (в первом условии)

     if (.not. full_range) then

          leftbound = 1
          rightbound = 5860

     else

          leftbound = 1
          rightbound = p_num

     endif

     ! Массив значений периодограммы
     allocate(I_p(leftbound:rightbound), stat = ier)
     if (ier .ne. 0) stop 'Не удалось выделить память для массива I_p'
     
     allocate(nu(leftbound:rightbound), stat = ier)
     if (ier .ne. 0) stop 'Не удалось выделить память для массива I_p'

     open(11, file="output2")
     do p = leftbound, rightbound, 1

          p_d = p
          p_cur = leftbound - 1d0 + p_d * p_step
          
          s1 = 0d0
          s2 = 0d0

          do t = 0, N_wif - 1

               j = N_index_array(t)

               j_d = j

               koef = 2d0 * pi * p_cur * j_d / N_wif
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

          I_p(p) = (s1 * s1 + s2 * s2) / ( N_wif * pi)
          nu(p) = p_cur / N_d
          
          if (abs(I_p(p)) .le. 1e-15) I_p(p) = 0d0

!          I_p(p) = (s1 * s1 + s2 * s2)/(N_d * pi)

          write(11,'(e16.7, 1x, e16.7, 1x, e16.7)') N_d / p_cur, nu(p), I_p(p)

     enddo
     close(11)

     ! [Вычисление корреляционной функции через периодограмму]
     
     cut_p_num = 293
     cut_p_num_d = cut_p_num
     leftbound_d = leftbound
     rightbound_d = rightbound
     
     cut_N = cut_p_num_d - leftbound_d + 1d0
     
     s2 = 0d0
     
     do p = leftbound, rightbound - 1, 1

!               p_d = p
                
!               p_cur = leftbound_d - 1d0 + p_d * p_step

!               ! Проверка на ошибку округления sin_value и cos_value
!               if (abs(cos_value) .le. 1e-3) cos_value = 0d0
               
               s2 = s2 + I_p(p)

     enddo
     
     allocate(C(0:p_num - 10), stat = ier)
     if (ier .ne. 0) stop 'Не удалось выделить память для массива I_p'
     
     open(12, file = 'output3')
     do t = 0, p_num - 10, 1
     
          t_d = t
          t_cur = 0d0 + t_d * p_step
!          t_cur = t
          
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
          
!          C(t) = s1 * N_d / (rightbound_d - leftbound_d + 1d0) / c_0
!          C(t) = C(t) * N_d / (N_d - t_cur + 1d0) / c_0
          C(t) = s1 * N_d / (N_d - t_cur) / s2

          if (abs(C(t)) .le. 1e-15) C(t) = 0d0

          write(12,'(e16.7, 1x, e16.7)') t_cur, C(t)
     
     enddo     
     close(12)

     deallocate(A, r, I_p, N_if_array, N_index_array, nu, C)

end
