module md_interactions
use md_general
use md_neighbours
use LennardJones
use LennardJones_1g
!use Morse
use graphenenorm
use LennardJonesCosine
use MorseCosine
use RosatoGuillopeLegrand
use TersoffBrenner
use REBOsolidcarbon
implicit none

type interaction_parameters
	type(LennardJones_parameters),allocatable			:: LJ(:)
	type(LennardJones1g_parameters),allocatable			:: LJ1g(:)
	!type(Morse_parameters),allocatable					:: Morse(:)
	type(LennardJonesCosine_parameters),allocatable		:: LJC(:)
	type(MorseCosine_parameters),allocatable			:: MorseC(:)
	type(RosatoGuillopeLegrand_parameters),allocatable	:: RJL(:)
	type(TersoffBrenner_parameters),allocatable			:: TB(:)
	type(REBOsc_parameters),allocatable					:: REBOsc(:)
end type interaction_parameters

type interaction
	integer:: nl_n,neib_order
	integer,allocatable:: group_nums(:)
	type(neighbour_list),allocatable:: nl(:)
	real:: energy
	character(len=32):: interaction_name,parameters_file
	type(interaction_parameters) :: parameters
	logical:: numerical_force
end type interaction	

contains

subroutine create_groups(groups,file_id,out_id,atoms)
	type(particle_group),allocatable:: groups(:)
	type(particles):: atoms
	character(len=32),allocatable::	type_names(:,:)
	character(len=128):: str,frmt
	integer:: i,particle_types_num,groups_num,file_id,out_id
	
	read(file_id,*) str,particle_types_num
	write(out_id,'(A32,i12)') str,particle_types_num
	read(file_id,*) str,groups_num
	write(out_id,'(A32,i12)') str,groups_num
	allocate(groups(groups_num),type_names(particle_types_num,groups_num))
	write(frmt,'("(i6,A,",i0,"A12,i9)")') particle_types_num
	do i=1,groups_num
		read(file_id,*) str,type_names(:,i)
		call create_particle_group(groups(i),type_names(:,i),atoms)
		write(out_id,frmt) i,' ',type_names(:,i),groups(i)%N
	enddo
	
end subroutine create_groups

subroutine create_interactions(interactions,groups,file_id,out_id,input_path)
	type(interaction),allocatable:: interactions(:)
	type(particle_group):: groups(:)
	integer:: file_id,out_id,i,j,inter_n
	character(len=128):: input_path,str
	
	read(file_id,*) str,inter_n
	write(out_id,'(A32,i12)') str,inter_n
	allocate(interactions(inter_n))
	do i=1,inter_n
		read(file_id,*) interactions(i)%interaction_name,interactions(i)%parameters_file
		select case (interactions(i)%interaction_name)
		case('lj')
			allocate(interactions(i)%parameters%LJ(1))
			call read_LJ_parameters(interactions(i)%parameters%LJ(1),trim(input_path)//interactions(i)%parameters_file)
			interactions(i)%nl_n = 2
			interactions(i)%neib_order = 0
			interactions(i)%numerical_force = .false.
		case('lj1g')
			allocate(interactions(i)%parameters%LJ1g(1))
			call read_LJ1g_parameters(interactions(i)%parameters%LJ1g(1),trim(input_path)//interactions(i)%parameters_file)
			interactions(i)%nl_n = 1
			interactions(i)%neib_order = 0
			interactions(i)%numerical_force = .false.
		!case('morse')
		!	allocate(interactions(i)%parameters%Morse(1))
		!	call read_Morse_parameters(interactions(i)%parameters%Morse(1),trim(input_path)//interactions(i)%parameters_file)
		!	interactions(i)%nl_n = 2
		!	interactions(i)%neib_order = 0
		!	interactions(i)%numerical_force = .false.
		case('ljc')
			allocate(interactions(i)%parameters%LJC(1))
			call read_LJC_parameters(interactions(i)%parameters%LJC(1),trim(input_path)//interactions(i)%parameters_file)
			interactions(i)%nl_n = 3
			interactions(i)%neib_order = 0
			interactions(i)%numerical_force = .false.
		case('morsec')
			allocate(interactions(i)%parameters%MorseC(1))
			call read_MorseC_parameters(interactions(i)%parameters%MorseC(1),trim(input_path)//interactions(i)%parameters_file)
			interactions(i)%nl_n = 3
			interactions(i)%neib_order = 0
			interactions(i)%numerical_force = .false.
		case('tb')
			allocate(interactions(i)%parameters%TB(1))
			call read_TB_parameters(interactions(i)%parameters%TB(1),trim(input_path)//interactions(i)%parameters_file)
			interactions(i)%nl_n = 1
			interactions(i)%neib_order = 2
			interactions(i)%numerical_force = .false.
		case('rebosc')
			allocate(interactions(i)%parameters%REBOsc(1))
			call read_REBOsc_parameters(interactions(i)%parameters%REBOsc(1),trim(input_path)//interactions(i)%parameters_file)
			interactions(i)%nl_n = 1
			interactions(i)%neib_order = 3
			interactions(i)%numerical_force = .true.			
		case('rjl')
			allocate(interactions(i)%parameters%RJL(1))
			call read_RJL_parameters(interactions(i)%parameters%RJL(1),trim(input_path)//interactions(i)%parameters_file)
			interactions(i)%nl_n = 1
			interactions(i)%neib_order = 2
			interactions(i)%numerical_force = .false.
		case default
			print*,'error: unknown interaction name'
		end select
		write(out_id,'(2A32,i6)') interactions(i)%interaction_name,interactions(i)%parameters_file,interactions(i)%nl_n
		allocate(interactions(i)%group_nums(2*interactions(i)%nl_n))
		allocate(interactions(i)%nl(interactions(i)%nl_n))
		do j=1,interactions(i)%nl_n
			read(file_id,*) interactions(i)%group_nums(j*2-1:j*2),&
			interactions(i)%nl(j)%neighb_num_max,interactions(i)%nl(j)%r_cut,interactions(i)%nl(j)%update_period
			write(out_id,'(3i6,f16.6,i9)') interactions(i)%group_nums(j*2-1:j*2),&
			interactions(i)%nl(j)%neighb_num_max,interactions(i)%nl(j)%r_cut,interactions(i)%nl(j)%update_period
			interactions(i)%nl(j)%N = groups(interactions(i)%group_nums(j*2-1))%N
			call create_neighbour_list(interactions(i)%nl(j))
		enddo
	enddo
	call allocate_graphene_norm(interactions)
	
end subroutine create_interactions

subroutine update_interactions_neighbour_lists(md_step,interactions,atoms,groups,cell,exe_time_nlsearch,exe_time_nldistance)
	type(interaction):: interactions(:)
	type(particles):: atoms
	type(particle_group):: groups(:)
	type(simulation_cell):: cell
	integer:: i,j,md_step
	real:: exe_time_nlsearch,exe_time_nldistance
	
	do i=1,size(interactions)
	
		call update_neighbour_list(md_step,interactions(i)%nl(1),atoms,&
		groups(interactions(i)%group_nums(1)),groups(interactions(i)%group_nums(2)),cell,exe_time_nlsearch,exe_time_nldistance)
		
		select case (interactions(i)%interaction_name)
		case('lj','morse','ljc','morsec')
			call converce_neighbour_list(interactions(i)%nl(2),&
			groups(interactions(i)%group_nums(2)),interactions(i)%nl(1))
		end select
		
		select case (interactions(i)%interaction_name)
		case('ljc','morsec')
			do j=1,size(interactions)
				select case (interactions(j)%interaction_name)
				case('tb'); exit
				case('rebosc'); exit
				end select
			enddo
			if (j<=size(interactions)) then
				call update_nearest_neighbours_in_graphene(md_step,interactions(i)%nl(3),interactions(j)%nl(1),atoms,&
				groups(interactions(j)%group_nums(1)),cell)
			else
				call update_neighbour_list(md_step,interactions(i)%nl(3),atoms,&
				groups(interactions(i)%group_nums(5)),groups(interactions(i)%group_nums(6)),cell,exe_time_nlsearch,exe_time_nldistance)
			endif
		end select
		
	enddo
	
	call update_norm_in_graphene(interactions)
	
end subroutine update_interactions_neighbour_lists

subroutine allocate_graphene_norm(interactions)
	type(interaction):: interactions(:)
	integer:: i
	
	do i=1,size(interactions)
		select case (interactions(i)%interaction_name)
		case('ljc')
			allocate(interactions(i)%parameters%LJC(1)%gr_norm(3,interactions(i)%nl(3)%N))
		case('morsec')
			allocate(interactions(i)%parameters%MorseC(1)%gr_norm(3,interactions(i)%nl(3)%N))
		end select
	enddo
	
end subroutine allocate_graphene_norm

subroutine update_norm_in_graphene(interactions)
	type(interaction):: interactions(:)
	integer:: i
	
	do i=1,size(interactions)
		select case (interactions(i)%interaction_name)
		case('ljc')
			call find_norm_in_graphene(interactions(i)%parameters%LJC(1)%gr_norm,interactions(i)%nl(3)%dr)
		case('morsec')
			call find_norm_in_graphene(interactions(i)%parameters%MorseC(1)%gr_norm,interactions(i)%nl(3)%dr)
		end select
	enddo

end subroutine update_norm_in_graphene

subroutine calculate_forces(atoms,interactions)
	type(interaction):: interactions(:)
	type(particles):: atoms
	integer:: i	
		
	do i=1,size(interactions)
		if(.not. interactions(i)%numerical_force) then
			select case (interactions(i)%interaction_name)
			case('lj')
				call LJ_forces(atoms,interactions(i)%nl(1),interactions(i)%parameters%LJ(1))
				call LJ_forces(atoms,interactions(i)%nl(2),interactions(i)%parameters%LJ(1))
			case('lj1g')
				call LJ1g_forces(atoms,interactions(i)%nl(1),interactions(i)%parameters%LJ1g(1))
			!case('morse')
			!	call Morse_forces(atoms,interactions(i)%nl(1),interactions(i)%parameters%Morse(1))
			!	call Morse_forces(atoms,interactions(i)%nl(2),interactions(i)%parameters%Morse(1))
			case('ljc')
				call LJC_forces_for_graphene(atoms,interactions(i)%nl(1),interactions(i)%nl(3),interactions(i)%parameters%LJC(1))
				call LJC_forces_for_other_atoms(atoms,interactions(i)%nl(2),interactions(i)%parameters%LJC(1))
			case('morsec')
				call MorseC_forces_for_graphene(atoms,interactions(i)%nl(1),interactions(i)%nl(3),interactions(i)%parameters%MorseC(1))
				call MorseC_forces_for_other_atoms(atoms,interactions(i)%nl(2),interactions(i)%parameters%MorseC(1))
			case('tb')
				call TB_forces(atoms,interactions(i)%nl(1),interactions(i)%parameters%TB(1))
			case('rebosc')
			!	
			case('rjl')
				call RJL_forces(atoms,interactions(i)%nl(1),interactions(i)%parameters%RJL(1))
			end select
		endif
	enddo
	
end subroutine calculate_forces

subroutine energy(inter_name,e,nl,p)
	type(interaction_parameters):: p
	type(neighbour_list):: nl
	character(len=32):: inter_name
	real:: e
	
	select case (inter_name)
	case('lj');		call LJ_energy(e,nl,p%LJ(1))
	case('lj1g');	call LJ1g_energy(e,nl,p%LJ1g(1))
	!case('morse');
	case('ljc');	call LJC_energy(e,nl,p%LJC(1))
	case('morsec');	call MorseC_energy(e,nl,p%MorseC(1))
	case('tb');		call TB_energy(e,nl,p%TB(1))
	case('rebosc');	call REBOsc_energy(e,nl,p%REBOsc(1))
	case('rjl');	call RJL_energy(e,nl,p%RJL(1))
	end select
	
end subroutine energy

subroutine calculate_potential_energies(interactions)
	type(interaction):: interactions(:)
	integer:: i
	
	do i=1,size(interactions)
		call energy(interactions(i)%interaction_name,interactions(i)%energy,interactions(i)%nl(1),interactions(i)%parameters)
	enddo
		
end subroutine calculate_potential_energies

subroutine calculate_forces_numerically(atoms,interactions)
	type(interaction):: interactions(:)
	type(particles):: atoms
	type(neighbour_list):: tnl
	real:: e1,e2
	integer:: i,k,inl,j
	real,parameter:: dx = 10.**(-6)

	do i=1,size(interactions)
		if(interactions(i)%numerical_force) then
			select case (interactions(i)%interaction_name)
			case('ljc','morsec')!nl(1) - calc_tnl(0),add nearest neibls,shift_nl,shift_gr_norm !nl(2) - calc_tnl(0),shift_nl
			case('lj','morse','tb','rebosc','rjl')
				do j=1,interactions(i)%nl_n
					!$OMP PARALLEL firstprivate(k,inl,e1,e2) private(tnl)
					call create_truncated_nl(tnl,interactions(i)%nl(j))
					!$OMP DO
					do inl=1,interactions(i)%nl(j)%N
						call calculate_truncated_nl(tnl,interactions(i)%nl(j),inl,interactions(i)%neib_order)
						do k=1,3
							call shift_drs(tnl,inl,k,interactions(i)%nl_n,-dx)
							call energy(interactions(i)%interaction_name,e1,tnl,interactions(i)%parameters)
							call shift_drs(tnl,inl,k,interactions(i)%nl_n,2*dx)
							call energy(interactions(i)%interaction_name,e2,tnl,interactions(i)%parameters)
							if(k/=3) call shift_drs(tnl,inl,k,interactions(i)%nl_n,-dx)
							atoms%forces(k,interactions(i)%nl(j)%particle_index(inl)) =&
							atoms%forces(k,interactions(i)%nl(j)%particle_index(inl))+(e1-e2)/2/dx
						enddo
					enddo
					!$OMP END DO
					call destroy_truncated_nl(tnl)
					!$OMP END PARALLEL
				enddo
			end select
		endif
	enddo
	
end subroutine

subroutine create_truncated_nl(tnl,nl)
	type(neighbour_list):: tnl,nl
	
	tnl%N = nl%N
	tnl%neighb_num_max = nl%neighb_num_max
	call create_neighbour_list(tnl)
	
end subroutine create_truncated_nl

subroutine destroy_truncated_nl(tnl)
	type(neighbour_list):: tnl
	
	tnl%N = 0
	tnl%neighb_num_max = 0
	if(allocated(tnl%particle_index)) deallocate(tnl%particle_index)
	if(allocated(tnl%nnum)) deallocate(tnl%nnum)
	if(allocated(tnl%nlist)) deallocate(tnl%nlist)
	if(allocated(tnl%dr)) deallocate(tnl%dr)
	if(allocated(tnl%moddr)) deallocate(tnl%moddr)
	
end subroutine destroy_truncated_nl

subroutine calculate_truncated_nl(tnl,nl,i,n)
	type(neighbour_list):: tnl,nl
	integer:: n,i,ni,p,q,j
	
	tnl%particle_index = 0; tnl%nnum = 0; !tnl%nlist = 0!tnl%dr = 0.!tnl%moddr = 0.
	
	tnl%N = nl%N
	tnl%neighb_num_max = nl%neighb_num_max
	tnl%nlist(:nl%nnum(i),i) = nl%nlist(:nl%nnum(i),i)
	tnl%dr(:,:nl%nnum(i),i) = nl%dr(:,:nl%nnum(i),i)
	tnl%moddr(:nl%nnum(i),i) = nl%moddr(:nl%nnum(i),i)
	tnl%nnum(i) = nl%nnum(i)
	tnl%particle_index(i) = nl%particle_index(i)
	do ni=1,n
		if(ni<n) then
			do j=1,tnl%N
				do p=1,tnl%nnum(j)
					q = tnl%nlist(p,j)
					if(tnl%nnum(q)==0) then
						tnl%nlist(:nl%nnum(q),q) = nl%nlist(:nl%nnum(q),q)
						tnl%dr(:,:nl%nnum(q),q) = nl%dr(:,:nl%nnum(q),q)
						tnl%moddr(:nl%nnum(q),q) = nl%moddr(:nl%nnum(q),q)
						tnl%particle_index(q) = nl%particle_index(q)
					endif
				enddo
			enddo
			do j=1,tnl%N
				if(tnl%particle_index(j)/=0) tnl%nnum(j) = nl%nnum(j)
			enddo
		else
			do j=1,tnl%N
				do p=1,tnl%nnum(j)
					q = tnl%nlist(p,j)
					if(tnl%nnum(q)==0) then
						tnl%nnum(q) = tnl%nnum(q)+1
						tnl%nlist(tnl%nnum(q),q) = j
						tnl%dr(:,tnl%nnum(q),q) = -nl%dr(:,p,j)
						tnl%moddr(tnl%nnum(q),q) = nl%moddr(p,j)
						tnl%particle_index(q) = nl%particle_index(q)
					endif
				enddo
			enddo
		endif
	enddo
	
end subroutine calculate_truncated_nl

subroutine shift_drs(tnl,inl,k,nl_n,dx)
	type(neighbour_list):: tnl
	real:: dx
	integer:: p,k,inl,q,j,nl_n
	
	do p=1,tnl%nnum(inl)
		tnl%dr(k,p,inl) = tnl%dr(k,p,inl)-dx
		tnl%moddr(p,inl) = sqrt(sum(tnl%dr(:,p,inl)**2))
	enddo
	if(nl_n==1) then
		do j=1,tnl%N
			do p=1,tnl%nnum(j)
				q = tnl%nlist(p,j)
				if(q==inl) then
					tnl%dr(k,p,j) = tnl%dr(k,p,j)+dx
					tnl%moddr(p,j) = sqrt(sum(tnl%dr(:,p,j)**2))
				endif
			enddo
		enddo
	endif

end subroutine shift_drs

subroutine shift_gr_norm(gr_norm,nl_nn,inl,k,dx)
	real:: gr_norm(:,:)
	type(neighbour_list):: nl_nn
	real:: dx,drj12(3),drj31(3)
	integer:: k,inl,p,q,j,l
	
	do p=1,nl_nn%nnum(inl)
		j = nl_nn%nlist(p,inl)
		do q=1,nl_nn%nnum(j)
			if(nl_nn%nlist(q,j)==inl) exit
		enddo	
		drj12 = nl_nn%dr(:,mod(q+1,3)+1,j)-nl_nn%dr(:,mod(q,3)+1,j)
		drj31 = nl_nn%dr(:,mod(q,3)+1,j)-nl_nn%dr(:,q,j)
		drj31(k) = drj31(k)+dx
		do l=1,3
			gr_norm(l,j) = (drj12(mod(l,3)+1))*(drj31(mod(l+1,3)+1))-(drj12(mod(l+1,3)+1))*(drj31(mod(l,3)+1))
		enddo
		if (gr_norm(3,j)<0.) gr_norm(:,j)=-gr_norm(:,j)
		gr_norm(:,j) = gr_norm(:,j)/sqrt(sum(gr_norm(:,j)**2))
	enddo
	
end subroutine shift_gr_norm

subroutine nlists_load(out_id,interactions)
	type(interaction):: interactions(:)
	integer:: i,j,out_id
	 
	do i=1,size(interactions)
		do j=1,interactions(i)%nl_n
			write(out_id,*) trim(interactions(i)%interaction_name),' ',j,' neib lists load:',&
			maxval(interactions(i)%nl(j)%nnum),'/',interactions(i)%nl(j)%neighb_num_max
		enddo
	enddo
	
end subroutine nlists_load

end module md_interactions