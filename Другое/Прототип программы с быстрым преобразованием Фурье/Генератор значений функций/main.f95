program main ! Генератор значений функции на равномерной сетке на промежутке
implicit none

     real(8) x    ! Аргументы функции
     real(8) y    ! Значения функции
     real(8) a, b ! Границы промежутка
     real(8) h    ! Шаг
     
     integer(4) i     ! Вспомогательная переменная
     integer(4) N     ! Число промежутков деления
     real(8) i_d, N_d ! Овеществления
     real(8) pi       ! Число pi

     pi = 4d0*datan(1d0)

     a = pi / 2d0
     b = 1001d0 * pi / 2d0
     
     N = 8191
     N_d = N

     h = (b - a) / N_d

     open(10, file = 'result')
     do i = 0, N, 1

          i_d = i
          x = a + h * i_d

          ! Указание желаемой функции
          y = dsin(x / 20d0) + 1000d0

          !write(10,'(i5, e17.7, 1x, e17.7)') i, x, y 
          !write(10,'(e17.7, 1x, e17.7)') x, y
          write(10,'(e23.15)') y 
     
     enddo
     close(10)

     !write(*,'(/,a, e16.7)') 'Шаг:', h
     !write(*,'(a,i4,/)') 'Число значений:', N + 1

end
