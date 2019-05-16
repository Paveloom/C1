program main
use subprograms
implicit none

     ! Общие переменные

     real(8), allocatable, dimension(:,:) :: A ! Матрица исходных данных

     integer(4) N      ! Размер выборки
     integer(4) N_if   ! Число исключений
     integer(4) N_wif  ! Размер выборки с исключениями (N - N_if)
     real(8)    x_mean ! Среднее значение выборки с учётом исключений
     real(8)    pi     ! Число pi

     ! Массив индексов-исключений
     integer(4), allocatable, dimension(:) :: N_if_array

     ! Массив индексов с учётом исключений
     integer(4), allocatable, dimension(:) :: N_index_array

     ! Периодограмма:

     real(8), allocatable, dimension(:) :: I_p ! Массив значений периодограммы
     
     real(8) p_step ! Шаг дискретизации частот
     integer p_num  ! Общее число множителей p
     
     ! Коррелограмма (через обратное преобразование)
     integer(4) t_koef ! Множитель дискретизации периодов

     ! Вспомогательные переменные
     integer k, i, ier
     
     ! Овеществления
     real(8) N_d, p_num_d, leftbound_d

     ! Переменные для выбора рабочего диапазона частот
     logical full_range    ! Использовать полный или указанный рабочие диапазоны? (.true., если полный)
     integer(4) leftbound  ! Левая граница рабочего диапазона частот
     integer(4) rightbound ! Праввая граница рабочего диапазона частот

     ! Указать размер выборки
     N = 5860
     N_d = N

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

     ! Считывание исходных данных

        do i = 0, N - 1

                read(*,*) A(i,1), A(i,2)

        enddo
        
     ! [Вычисление среднего значения выборки]
     call F1_mean(A, x_mean, N_index_array, N_wif, N)

     ! [Вычисление коррелограммы прямым способом]
     call F2_correlogram_direct(A, x_mean, N)

     ! Определение pi
     pi = 4d0 * datan(1d0)

     ! Вычисление шага дисретизации

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
     
     leftbound_d = leftbound

     ! Массив значений периодограммы
     allocate(I_p(leftbound:rightbound), stat = ier)
     if (ier .ne. 0) stop 'Не удалось выделить память для массива I_p'

     ! [Вычисление периодограммы]
     
     t_koef = 10
     call F3_periodogram(A, leftbound, leftbound_d, rightbound, p_step, N_index_array, N, N_d, N_wif, x_mean, pi, I_p)

     ! [Вычисление коррелограммы через обратное преобразование Фурье к периодограмме]
     
     call F4_correlogram_fourier_transform(I_p, leftbound, leftbound_d, rightbound, p_num, p_step, t_koef, N, N_d, pi)
     
     deallocate(A, I_p, N_if_array, N_index_array)

end
