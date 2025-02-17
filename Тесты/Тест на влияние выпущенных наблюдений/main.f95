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
     
     integer(4) p_num  ! Общее число множителей p
     integer(4) p_koef ! Множитель дискретизации частот
     real(8)    p_step ! Шаг дискретизации частот
     
     ! Коррелограмма (через обратное преобразование)
     integer(4) t_koef ! Множитель дискретизации периодов

     ! Вспомогательные переменные
     integer(4) k, i, ier
     
     ! Овеществления
     real(8) N_d, p_num_d, leftbound_d

     ! Переменные для выбора рабочего диапазона частот
     integer(4) leftbound  ! Левая граница рабочего диапазона частот
     integer(4) rightbound ! Праввая граница рабочего диапазона частот
     
     ! Ответы на вопросы
     integer(4) full_range ! Использовать полный или указанный рабочие диапазоны?
     integer(4) use_if     ! Использовать массив исключений?
     integer(4) bias_fix   ! Приводить коэффициенты корреляции к несмещённой оценке?

     open(9, file = 'settings')
     
     ! Считывание размера выборки
     read(9,'()'); read(9,*) N
     
     ! Использовать массив исключений?
     read(9,'(//)'); read(9,*) use_if
     
     ! Считывание числа исключений
     read(9,'(/)'); read(9,*) N_if
     
     ! Использовать указанный диапазон частот для вычисления периодограммы
     ! или считать по полному диапазону?
     read(9,'(///)'); read(9,*) full_range
     
     ! Приводить коэффициенты корреляции к несмещённой оценке?
     read(9,'(//)'); read(9,*) bias_fix
     
     ! Считывание множителя дискретизации множителей p (для частот)
     read(9,'(/)'); read(9,*) p_koef
     
     ! Считывание множителя дискретизации множителей t (для периодов)
     read(9,'(/)'); read(9,*) t_koef
     
     ! Считывание границ рабочего диапазона частот
     read(9,'(///)'); read(9,*) leftbound
                      read(9,*) rightbound
     
     close(9)

     ! Овеществление N
     N_d = N
     
     ! Настройки рабочего диапазона частот в разделе периодограмма

     ! Вычисление размера выборки с исключениями
     N_wif = N - N_if

     if (use_if .eq. 0) then

     ! Массив индексов-исключений
     allocate(N_if_array(N_if), stat = ier)
     if (ier .ne. 0) stop 'Не удалось выделить память для массива N_if_array'

     ! Массив индексов с учётом исключений
     allocate(N_index_array(0:N_wif-1), stat = ier)
     if (ier .ne. 0) stop 'Не удалось выделить память для массива N_index_array'

     ! Заполнение массива индексов-исключений
     N_if_array = (/ 4, 5, 6, 11, 12, 13, 14, 15, 27, 28, 29, 48, 49, 50, 51, 52, 53, &
     &62, 79, 80, 81, 82, 88, 92, 99 /)
!     N_if_array = (/ 11, 48, 81, 88, 99 /)

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
     
     deallocate(N_if_array)
     
     endif
     
     ! Выделение памяти под рабочие массивы

     ! Исходные данные
     allocate(A(0:N-1,2), stat = ier)
     if (ier .ne. 0) stop 'Не удалось выделить память для массива A'

     ! Считывание исходных данных

        do i = 0, N - 1

                read(*,*) A(i,1), A(i,2)

        enddo
        
     ! [Вычисление среднего значения выборки]
     call F1_mean(A, x_mean, N_index_array, N_wif, N, N_d, use_if)

     ! [Вычисление коррелограммы прямым способом]
     call F2_correlogram_direct(A, x_mean, N, N_d, bias_fix)

     ! Определение pi
     pi = 4d0 * datan(1d0)

     ! Вычисление шага дисретизации для множителей частот

     p_num = N * p_koef
     p_num_d = p_num
     p_step = N_d / p_num_d

     ! Использование (если указано) полного рабочего диапазона частот

     if (full_range .eq. 0) then

          leftbound = 1
          rightbound = p_num

     endif
     
     leftbound_d = leftbound

     ! Массив значений периодограммы
     allocate(I_p(leftbound:rightbound), stat = ier)
     if (ier .ne. 0) stop 'Не удалось выделить память для массива I_p'

     ! [Вычисление периодограммы]
     call F3_periodogram(A, leftbound, leftbound_d, rightbound, p_step, N_index_array, N, N_d, N_wif, x_mean, pi, I_p, use_if)
     deallocate(A)
     if (use_if .eq. 0) deallocate(N_index_array)

     ! [Вычисление коррелограммы через обратное преобразование Фурье к периодограмме]
     call F4_correlogram_fourier_transform(I_p, leftbound, leftbound_d, rightbound, p_num, p_step, t_koef, N, N_d, pi, &
     &bias_fix)
     deallocate(I_p)

end
