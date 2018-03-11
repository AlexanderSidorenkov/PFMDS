!> \brief Модуль настройки OpenMP параллелизма.
!> \warning omp_chunk_size пока нигде не используется.
module perfomance_settings
implicit none

public

integer:: omp_chunk_size !< Массивы данных будут разбиваться на части длины omp_chunk_size при использовании omp parallel for schedule(dynamic ,omp_chunk_size)

contains

!> Настраивает количество OpenMP потоков и omp_chunk_size.
subroutine set_openmp_perfomance(num_of_omp_treads,N)
	integer:: chunks_per_thread !< Количество частей массива на поток при использовании omp parallel for schedule(dynamic ,omp_chunk_size)
	integer:: num_of_omp_treads !< Количество OpenMP потоков
	integer:: N !< Количество частиц в моделируемой системе
	
	call omp_set_num_threads(num_of_omp_treads)
	!call omp_thread_limit(num_of_omp_treads)
	
	chunks_per_thread = 8
	omp_chunk_size = int(real(N)/chunks_per_thread/num_of_omp_treads)
	if(omp_chunk_size==0) omp_chunk_size = 1
	
end subroutine

end module perfomance_settings