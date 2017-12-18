module perfomance_settings
implicit none

public

integer:: chunk_size

contains

subroutine set_openmp_perfomance(num_of_omp_treads,N)
	integer:: chunks_per_thread,num_of_omp_treads,N
	
	chunks_per_thread = 8
	call omp_set_num_threads(num_of_omp_treads)
	!call omp_thread_limit(num_of_omp_treads)
	chunk_size = int(real(N)/chunks_per_thread/num_of_omp_treads)
	
end subroutine

end module perfomance_settings