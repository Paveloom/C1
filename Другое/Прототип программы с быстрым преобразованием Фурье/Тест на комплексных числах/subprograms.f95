module subprograms
implicit none

     contains
     
     subroutine DFT(A, B, N)
     implicit none
     
     real(8), intent(in), dimension(0:) :: A    ! Входной массив данных
     real(8), intent(inout), dimension(0:) :: B ! Выходной массив данных
     integer(4), intent(in) :: N                ! Исходный размер массива данных
     
     integer(4) i, k ! Вспомогательные переменные
     real(8) koef          ! Аргумент тригонометрических функций
     real(8) cos_value     ! Значение косинуса от аргумента koef
     real(8) pi            ! Число pi
     real(8) i_d, k_d, N_d ! Овеществления i, k и N
     
     ! Определение pi
     pi = 4d0 * datan(1d0)
     
     N_d = N
     
     do k = 0, N - 1, 1
     
          k_d = k
          
          do i = 0, N - 1, 1
          
               i_d = i
               
               koef = 2d0 * pi * k_d * i_d / N_d
               cos_value = dcos(koef)
                         
               ! Проверка на ошибку округления cos_value
               if (abs(cos_value) .le. 1e-3) cos_value = 0d0
          
               B(k) = B(k) + A(i) * cos_value
          
          enddo
     
     enddo
     
     end subroutine



     recursive subroutine FFT(C) ! Процедура быстрого преобразования Фурье
     implicit none
     
     complex(8), intent(inout), dimension(0:) :: C ! Входной / выходной массив данных
     
     complex(8), allocatable, dimension(:) :: odd, even ! Массивы данных по нечётным и чётным индексам
     
     complex(8) wn     ! Комплексное число cos(koef) + i * sin(koef)
     complex(8) w      ! Держатель значений корня из единицы
     real(8) arg       ! Аргумент тригонометрических функций
     real(8) pi        ! Число pi
     integer(4) N      ! Размер входного массива
     integer(4) N_half ! Половина от N
     
     integer(4) ier, i ! Вспомогательные переменные
     real(8) N_d  ! Овеществления i и N
     
     ! Определение pi
     pi = 4d0 * datan(1d0)
     
     N = size(C)
     
     if (N .eq. 1) return
     
     N_half = N / 2
     
     allocate(odd(0:N_half - 1), stat = ier)
     if (ier .ne. 0) stop 'Не удалось выделить память для массива odd'
     
     allocate(even(0:N_half - 1), stat = ier)
     if (ier .ne. 0) stop 'Не удалось выделить память для массива even'
     
     even(0:N_half - 1)  = C(0:N - 1:2)
     odd(0:N_half - 1)   = C(1:N - 1:2)
     
     call FFT(even)
     call FFT(odd)
     
     N_d = N
     
     arg = 2d0 * pi / N_d
     
     w  = complex(1d0, 0d0)
     wn = complex(cos(arg), sin(arg))
     
     do i = 0, N_half - 1, 1
     
          ! "Преобразование бабочки"
          C(i) = even(i) + w * odd(i)
          C(i + N_half) = even(i) - w * odd(i)
          
          w = w * wn
     
     enddo
     
     deallocate(odd, even)
     
     end subroutine

end
