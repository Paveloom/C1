program main
use subprograms
implicit none

     real(8), allocatable, dimension(:) :: A    ! Исходный массив данных
     real(8), allocatable, dimension(:) :: B    ! Выходной массив данных
     complex(8), allocatable, dimension(:) :: C ! Массив комплексных чисел
     integer(4) N / 4096 / ! Исходный размер массива данных
     integer(4) ier, i     ! Вспомогательные переменные
     
     allocate(A(0:N - 1), stat = ier)
     if (ier .ne. 0) stop 'Не удалось выделить память для массива A'
     
     allocate(B(0:N - 1), stat = ier)
     if (ier .ne. 0) stop 'Не удалось выделить память для массива B'
     
     allocate(C(0:N - 1), stat = ier)
     if (ier .ne. 0) stop 'Не удалось выделить память для массива С'
     
     read(*,*) A ! Считывание исходного массива данных
     
     do i = 0, N - 1, 1
     
          C(i) = complex(A(i), 0d0) ! Приведение массива вещественных чисел к комплексному
     
     enddo
     
     B = 0d0
     
     call DFT(A, B, N) ! Вызов процедуры дискретного преобразования Фурье
     
     open(11, file = 'DFT_result.dat') ! Вывод результата
     write(11, '(e23.15)') B
     
     deallocate(A, B)
     
     call FFT(C) ! Вызов процедуры быстрого преобразования Фурье
     
     open(12, file = 'FFT_result.dat') ! Вывод результата
     write(12, '(e23.15)') dble(C)
     
     deallocate(C)

end
