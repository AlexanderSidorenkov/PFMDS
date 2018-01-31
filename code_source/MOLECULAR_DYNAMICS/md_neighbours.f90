module md_neighbours
use md_general
implicit none

contains

subroutine create_neighbour_list(nl)
	type(neighbour_list):: nl

	if (nl%N>0) then
		allocate(nl%particle_index(nl%N))
		allocate(nl%nnum(nl%N))
		allocate(nl%lessnnum(nl%N))
		allocate(nl%nlist(nl%neighb_num_max,nl%N))
		allocate(nl%dr(3,nl%neighb_num_max,nl%N))
		allocate(nl%moddr(nl%neighb_num_max,nl%N))
		nl%particle_index = 0
		nl%nnum = 0
		nl%lessnnum = 0
		nl%nlist = 0
		nl%dr = 0.
		nl%moddr = 0.
	endif
	
	return
end subroutine create_neighbour_list

subroutine update_neighbour_list(md_step,nl,atoms,group1,group2,box,exe_time_nlsearch,exe_time_nldistance)
	type(particles)::	atoms
	type(particle_group):: group1,group2
	type(simulation_cell):: box
	type(neighbour_list):: nl
	integer:: md_step
	real:: exe_t,exe_time_nlsearch,exe_time_nldistance,omp_get_wtime

	if (group1%N>0 .and. group2%N>0 .and. nl%N>0) then
		if (mod(md_step,nl%update_period)==0) then
			exe_t = omp_get_wtime()
			call find_neighbours(nl,atoms,group1,group2,box)
			exe_time_nlsearch = exe_time_nlsearch+omp_get_wtime()-exe_t
		else
			exe_t = omp_get_wtime()
			call find_neighbour_distances(nl,atoms,group1,group2,box)
			exe_time_nldistance = exe_time_nldistance+omp_get_wtime()-exe_t
		endif
	endif
	
end subroutine update_neighbour_list

subroutine find_neighbours(nl,atoms,group1,group2,box)
	type(particles)::	atoms
	type(particle_group):: group1,group2
	type(simulation_cell):: box
	type(neighbour_list):: nl
	real:: dr(3),dr2
	integer:: i,j,k,ind,jnd,nnumind,lessnnumind,group1N,group2N

	group1N = size(group1%indexes)
	group2N = size(group2%indexes)
	
	if (group1N/=nl%N) then; write(*,*) 'error: group1N/=nl%N',group1N,nl%N; stop; endif
	!$OMP PARALLEL firstprivate(i,j,k,ind,jnd,dr,dr2,nnumind,lessnnumind)
	!$OMP DO
	do ind=1,group1N
		i = group1%indexes(ind)
		nl%particle_index(ind) = group1%indexes(ind)
		!nl%nlist(:nl%nnum(ind),ind) = 0
		nnumind = 0
		lessnnumind = -1
		do jnd=1,group2N
			j = group2%indexes(jnd)
			if (i/=j) then
				call find_distance(dr,dr2,atoms%positions(:,i),atoms%positions(:,j),box)
				if (dr2<nl%r_cut*nl%r_cut) then
					if (lessnnumind==-1 .and. i<j) lessnnumind = nnumind
					nnumind = nnumind+1
					if (nnumind>nl%neighb_num_max) then; write(*,*) 'error: too many neighbours',ind,nnumind; stop; endif
					nl%nlist(nnumind,ind) = jnd
					nl%dr(1,nnumind,ind) = dr(1)
					nl%dr(2,nnumind,ind) = dr(2)
					nl%dr(3,nnumind,ind) = dr(3)
					nl%moddr(nnumind,ind) = sqrt(dr2)
				endif
			endif
		enddo
		if(lessnnumind==-1) then
			nl%lessnnum(ind) = nnumind
		else
			nl%lessnnum(ind) = lessnnumind
		endif
		nl%nnum(ind) = nnumind
	enddo
	!$OMP END DO
	!$OMP END PARALLEL

	return
end subroutine find_neighbours

subroutine find_neighbour_distances(nl,atoms,group1,group2,box)
	type(particles)::	atoms
	type(particle_group):: group1,group2
	type(simulation_cell):: box
	type(neighbour_list):: nl
	real:: dr2
	integer:: i,p

	!$OMP PARALLEL firstprivate(i,p,dr2)
	!$OMP DO
	do i=1,nl%N
		do p=1,nl%nnum(i)
			call find_distance(nl%dr(:,p,i),dr2,atoms%positions(:,group1%indexes(i)),atoms%positions(:,group2%indexes(nl%nlist(p,i))),box)
			nl%moddr(p,i) = sqrt(dr2)
		enddo
	enddo
	!$OMP END DO
	!$OMP END PARALLEL

	return
end subroutine find_neighbour_distances

subroutine converce_neighbour_list(cnl,group,nl)
	type(neighbour_list):: nl,cnl
	type(particle_group):: group
	integer::		i,j,p,groupN

	groupN = size(group%indexes)
	
	if (groupN>0 .and. nl%N>0 .and. cnl%N>0) then
		if (groupN/=cnl%N) then; write(*,*) 'error: groupN/=cnl%N',groupN,cnl%N; stop; endif
		!$OMP PARALLEL firstprivate(i)
		!$OMP DO
		do i=1,group%N
			cnl%particle_index(i) = group%indexes(i)
		enddo
		!$OMP END DO
		!$OMP DO
		do i=1,cnl%N
			cnl%nnum(i) = 0
		enddo
		!$OMP END DO
		!$OMP END PARALLEL
		do i=1,nl%N	!can not be parallel!
			do p=1,nl%nnum(i)
				j = nl%nlist(p,i)
				cnl%nnum(j) = cnl%nnum(j)+1
				if (cnl%nnum(j)>cnl%neighb_num_max) then; write(*,*) 'error: too many neighbours',j,cnl%nnum(j); stop; endif
				cnl%nlist(cnl%nnum(j),j) = i
				cnl%dr(:,cnl%nnum(j),j) = -nl%dr(:,p,i)
				cnl%moddr(cnl%nnum(j),j) = nl%moddr(p,i)
			enddo
		enddo
	endif
	
	return
end subroutine converce_neighbour_list

end module md_neighbours