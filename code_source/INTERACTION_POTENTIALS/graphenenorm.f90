module graphenenorm
use md_general
use md_neighbours
implicit none

contains

subroutine find_gr_nearest_neighbors(nl_nn,nl)
	type(neighbour_list):: nl,nl_nn
	integer:: nnum_nn,i,p,k	
	
	nnum_nn = 3
	if(nl_nn%neighb_num_max/=nnum_nn) then; write(*,*) 'error: nl_nn%neighb_num_max/=nnum_nn'; stop; endif
	!$OMP PARALLEL firstprivate(i,p,k)
	!$OMP DO
	do i=1,nl%N
		k = 0
		do p=1,nl%nnum(i)
			if (nl%moddr(p,i)<nl_nn%r_cut) then
				k = k+1
				if (k<=nnum_nn) then
					nl_nn%nlist(k,i) = nl%nlist(p,i)
					nl_nn%moddr(k,i) = nl%moddr(p,i)
					nl_nn%dr(:,k,i) = nl%dr(:,p,i)
				else
					write(*,*) 'error: too many gr nearest neibs',i,p,nl%moddr(p,i); stop;
				endif
			endif
		enddo
		if (k/=nnum_nn) then; write(*,*) 'error: not enough gr nearest neibs',i,nl_nn%moddr(:,i),nl_nn%nlist(:,i); stop; endif;
		nl_nn%nnum(i) = nnum_nn
	enddo
	!$OMP END DO
	!$OMP END PARALLEL
	
end subroutine find_gr_nearest_neighbors

subroutine find_norm_in_graphene(gr_norm,dr_nn)
	integer:: nnum_nn,i,p,k
	real:: dr_nn(:,:,:),gr_norm(:,:),drj12(3),drj31(3)
	
	nnum_nn = 3
	!$OMP PARALLEL firstprivate(i,p,k,drj12,drj31)
	!$OMP DO
	do i=1,size(gr_norm(1,:))
		drj12 = dr_nn(:,2,i)-dr_nn(:,1,i)
		drj31 = dr_nn(:,1,i)-dr_nn(:,3,i)
		do k=1,3
			gr_norm(k,i) = (drj12(mod(k,3)+1))*(drj31(mod(k+1,3)+1))-(drj12(mod(k+1,3)+1))*(drj31(mod(k,3)+1))
		enddo
		if (gr_norm(3,i)<0.) gr_norm(:,i)=-gr_norm(:,i)
		gr_norm(:,i) = gr_norm(:,i)/sqrt(sum(gr_norm(:,i)**2))
	enddo
	!$OMP END DO
	!$OMP END PARALLEL	
end subroutine find_norm_in_graphene

subroutine update_nearest_neighbours_in_graphene(md_step,nl_nn,nl,atoms,group,box)
	type(particles)::	atoms
	type(neighbour_list):: nl,nl_nn
	type(particle_group):: group
	type(simulation_cell):: box
	integer:: md_step

	if (mod(md_step,nl%update_period)==0) then
		call find_gr_nearest_neighbors(nl_nn,nl)
	else
		call find_neighbour_distances(nl_nn,atoms,group,group,box)
	endif
	
end subroutine update_nearest_neighbours_in_graphene

end module graphenenorm