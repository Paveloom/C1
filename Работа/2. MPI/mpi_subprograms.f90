module mpi_subprograms
use mpi
implicit none

     contains

     ! [Процедура вычисления размеров порций и их границ]
     subroutine partion_sizes(full_size, full_leftbound, mpiSize, mpiRank, partion_size, partion_shift, partion_size_mod, partion_leftbound, partion_rightbound)
     implicit none
     
     integer(4), intent(in) :: full_size      ! Общий размер пространства итераций
     integer(4), intent(in) :: full_leftbound ! Левая граница пространства итераций (первый индекс)
     integer(4), intent(in) :: mpiSize        ! Размер коммуникатора 
     integer(4), intent(in) :: mpiRank        ! Ранг процесса
     
     integer(4), intent(inout) :: partion_size       ! Размер порции
     integer(4), intent(inout) :: partion_shift      ! Величина mpiRank * partion_size
     integer(4), intent(inout) :: partion_size_mod   ! Остаток от деления full_size на mpiSize
     integer(4), intent(inout) :: partion_leftbound  ! Левая граница порции для ранга mpiRank
     integer(4), intent(inout) :: partion_rightbound ! Правая граница порции для ранга mpiRank
     
     partion_size_mod = mod(full_size, mpiSize)
     
     if (mpiRank .eq. mpiSize - 1) then
               
          partion_size = full_size / mpiSize + partion_size_mod
          partion_shift = mpiRank * (partion_size - partion_size_mod)
               
     else
     
          partion_size = full_size / mpiSize
          partion_shift = mpiRank * partion_size
          
     endif
     
     partion_leftbound = full_leftbound + partion_shift
     partion_rightbound = full_leftbound + partion_size + partion_shift - 1
     
     end subroutine

end
