module md_general
use IFPORT
implicit none

type time_steps
	real	ts(4),simulation_time
end type time_steps

type simulation_cell
	real:: box_size(3),half_box_size(3)
end type simulation_cell

type particles
	integer:: N
	real,allocatable::	positions(:,:),velocities(:,:),masses(:),forces(:,:)
	character(len=32),allocatable::	atom_types(:)
end type particles

type particle_group
	integer:: N
	integer,allocatable:: indexes(:)
end type particle_group

type nose_hoover_chain
	integer:: M,L
	real,allocatable:: x(:),v(:),q(:)
	real::	temperature,s,e
end type nose_hoover_chain

type neighbour_list
	integer:: N,neighb_num_max,update_period
	real:: r_cut
	integer,allocatable:: nlist(:,:),nnum(:),lessnnum(:),particle_index(:)
	real,allocatable:: dr(:,:,:),moddr(:,:)
end type neighbour_list

type integrator_params
	integer:: l,period_snapshot,period_log
	real:: dt
	character(len=32):: int_name
end type integrator_params

contains!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine init_time_steps(dt,delta_t)
	type(time_steps):: dt
	integer:: i
	real::	delta_t

	do i=1,size(dt%ts)
		dt%ts(i) = delta_t/2**(i-1)
	enddo

	return
end subroutine init_time_steps

subroutine create_particle_group(group,type_names,atoms)
	type(particles)::	atoms
	type(particle_group)::	group
	character(len=32):: type_names(:)
	integer:: i,j,ind

	group%N = 0
	do j=1,size(type_names)
		group%N = group%N+count(atoms%atom_types==type_names(j))
	enddo
	allocate(group%indexes(group%N))

	ind=0
	do j=1,size(type_names)
		do i=1,atoms%N
			if (atoms%atom_types(i)==type_names(j)) then
				ind=ind+1
				group%indexes(ind) = i
			endif
		enddo
	enddo

	return
end subroutine create_particle_group

subroutine change_particle_group_N(group,md_step,change_ts1,change_ts2,change_frec,init_group)
	type(particle_group)::	group,init_group
	integer:: md_step,change_ts1,change_ts2,change_frec

	if(md_step<=change_ts1) then
		group%N = init_group%N
		if(md_step==change_ts1) group%N = group%N+1
	else
		if(md_step<change_ts2 .and. mod(md_step-change_ts1,change_frec)==0) group%N = group%N+1
	endif
	if(group%N>size(group%indexes)) group%N = size(group%indexes)
	
end subroutine change_particle_group_N

subroutine scale_velocities(atoms,group,s)
	type(particles)::	atoms
	type(particle_group):: group
	real:: s
	integer::		i,ind

	!$OMP PARALLEL firstprivate(i,ind)
	!$OMP DO
	do ind=1,group%N
		i = group%indexes(ind)
		atoms%velocities(:,i) = atoms%velocities(:,i)*s
	enddo
	!$OMP END DO
	!$OMP END PARALLEL

	return
end subroutine scale_velocities

subroutine random_velocities(atoms,group,rand_seed)
	type(particles)::	atoms
	type(particle_group):: group
	real a1,a2,b
	integer::		i,k,ind,rand_seed
	
	!$OMP PARALLEL firstprivate(i,ind,k,a1,a2,b)
	call srand(rand_seed)
	!$OMP DO
	do ind=1,group%N
		i = group%indexes(ind)
		do k=1,3,1
			b = 2.
			do while (b>=1.)
				a1 = 2.*rand()-1.
				a2 = 2.*rand()-1.
				b = a1**2+a2**2
			enddo
			atoms%velocities(k,i) = a1*sqrt(-2.*log(b)/b)
		enddo
	enddo
	!$OMP END DO
	!$OMP END PARALLEL

	return
end subroutine random_velocities

subroutine random_momenta(atoms,group,rand_seed)
	type(particles)::	atoms
	type(particle_group):: group
	real,parameter:: coef = 1.3806488/1.6605389217*10.**(-6)
	integer::		i,ind,rand_seed

	call random_velocities(atoms,group,rand_seed)

	!$OMP PARALLEL firstprivate(i,ind)
	!$OMP DO
	do ind=1,group%N
		i = group%indexes(ind)
		atoms%velocities(:,i) = atoms%velocities(:,i)*sqrt(coef/atoms%masses(i))
	enddo
	!$OMP END DO
	!$OMP END PARALLEL

	return
end subroutine random_momenta

subroutine calculate_kinetic_energy(ke,atoms,group)
	type(particles)::	atoms
	type(particle_group):: group
	real:: ke,ke_priv
	real,parameter:: mass_coef=1.6605389217/1.6021765654*10.**(2)
	integer::		i,ind

	ke = 0.
	ke_priv = 0.
	!$OMP PARALLEL firstprivate(ke_priv,i,ind)
	!$OMP DO
	do ind=1,group%N
		i = group%indexes(ind)
		ke_priv = ke_priv+atoms%masses(i)*sum(atoms%velocities(:,i)**2)/2*mass_coef
	enddo
	!$OMP END DO
	!$OMP ATOMIC
		ke = ke+ke_priv
	!$OMP END PARALLEL

	return
end subroutine calculate_kinetic_energy

subroutine calculate_mass_center(mc,atoms,group)
	type(particles)::	atoms
	type(particle_group):: group
	real:: mc(3),mc_priv(3),totm
	integer::		i,ind,k

	mc = 0.
	mc_priv = 0.
	!$OMP PARALLEL firstprivate(mc_priv,i,ind,k)
	!$OMP DO
	do ind=1,group%N
		i = group%indexes(ind)
		mc_priv = mc_priv+atoms%masses(i)*atoms%positions(:,i)
	enddo
	!$OMP END DO
	do k=1,3
		!$OMP ATOMIC
			mc(k) = mc(k)+mc_priv(k)
	enddo
	!$OMP END PARALLEL
	call calculate_masses_sum(totm,atoms,group)
	mc = mc/totm

	return
end subroutine calculate_mass_center

subroutine calculate_mass_center_velosity(mcv,atoms,group)
	type(particles)::	atoms
	type(particle_group):: group
	real:: mcv(3),mcv_priv(3),totm
	integer::		i,ind,k

	mcv = 0.
	mcv_priv = 0.
	!$OMP PARALLEL firstprivate(mcv_priv,i,ind,k)
	!$OMP DO
	do ind=1,group%N
		i = group%indexes(ind)
		mcv_priv = mcv_priv+atoms%masses(i)*atoms%velocities(:,i)
	enddo
	!$OMP END DO
	do k=1,3
		!$OMP ATOMIC
			mcv(k) = mcv(k)+mcv_priv(k)
	enddo
	!$OMP END PARALLEL
	call calculate_masses_sum(totm,atoms,group)
	mcv = mcv/totm

	return
end subroutine calculate_mass_center_velosity

subroutine zero_momentum(atoms,group)
	type(particles)::	atoms
	type(particle_group):: group
	integer::		i,ind
	real:: mcv(3)
	
	call calculate_mass_center_velosity(mcv,atoms,group)
	!$OMP PARALLEL firstprivate(mcv,i,ind)
	!$OMP DO
	do ind=1,group%N
		i = group%indexes(ind)
		atoms%velocities(:,i) = atoms%velocities(:,i)-mcv
	enddo
	!$OMP END DO
	!$OMP END PARALLEL
	
	return
end subroutine

subroutine calculate_masses_sum(totm,atoms,group)
	type(particles)::	atoms
	type(particle_group):: group
	real:: totm,totm_priv
	integer::		i,ind

	totm = 0.
	totm_priv = 0.
	!$OMP PARALLEL firstprivate(totm_priv,i,ind)
	!$OMP DO
	do ind=1,group%N
		i = group%indexes(ind)
		totm_priv = totm_priv+atoms%masses(i)
	enddo
	!$OMP END DO
	!$OMP ATOMIC
		totm = totm+totm_priv
	!$OMP END PARALLEL

	return
end subroutine calculate_masses_sum

subroutine calculate_force_sum(fs,atoms,group)
	type(particles)::	atoms
	type(particle_group):: group
	real:: fs(3),fs_priv(3)
	integer::		i,ind,k

	fs = 0.
	fs_priv = 0.
	!$OMP PARALLEL firstprivate(fs_priv,i,ind,k)
	!$OMP DO
	do ind=1,group%N
		i = group%indexes(ind)
		fs_priv = fs_priv+atoms%forces(:,i)
	enddo
	!$OMP END DO
	do k=1,3
		!$OMP ATOMIC
			fs(k) = fs(k)+fs_priv(k)
	enddo
	!$OMP END PARALLEL

	return
end subroutine calculate_force_sum

subroutine calculate_temperature(temp,ke,atoms,group)
	type(particles)::	atoms
	type(particle_group):: group
	real,parameter::	kt_a_degree =1.3806488/1.6021765654*10.**(-4)
	real:: temp,ke
	
	call calculate_kinetic_energy(ke,atoms,group)
	temp = 2*ke/kt_a_degree/(3*(group%N))
	
	return
end subroutine calculate_temperature

subroutine set_new_temperature(atoms,group,temp,rand_seed)
	type(particles)::	atoms
	type(particle_group):: group
	real:: temp,ke,temperature
	integer:: rand_seed
	
	call random_momenta(atoms,group,rand_seed)
	call zero_momentum(atoms,group)
	call calculate_temperature(temperature,ke,atoms,group)
	call scale_velocities(atoms,group,sqrt(temp/temperature))

end subroutine set_new_temperature

subroutine check_positions(out_id,atoms,box)
	type(particles)::	atoms
	type(simulation_cell):: box
	integer:: out_id,i,k,p
	real:: tolerance=0.0000001
	
	!$OMP PARALLEL private(i,k,p)
	p=0
	!$OMP DO
	do i=1,atoms%N
		do k=1,3
			if( .not.(atoms%positions(k,i)>(0.-tolerance) .and. atoms%positions(k,i)<(box%box_size(k)+tolerance)) ) then
				write(out_id,*) i,' particle out of cell ',atoms%positions(:,i)
				p = p+1
			endif
		enddo
	enddo
	!$OMP END DO
	!$OMP BARRIER
	if(p>0) stop
	!$OMP END PARALLEL
	
end subroutine check_positions

subroutine check_velocities(out_id,atoms)
	type(particles)::	atoms
	integer:: out_id,i
	real:: maxvel2,v2

	maxvel2 = -1.	
	!$OMP PARALLEL DO private(i,v2) REDUCTION(max:maxvel2)
	do i=1,atoms%N
		v2 = sum(atoms%velocities(:,i)**2)
		if(maxvel2<v2) maxvel2 = v2
	enddo
	!$OMP END PARALLEL DO
	write(out_id,'(A,f16.8,A)') ' max velocity: ',sqrt(maxvel2),' A/fs '
	
end subroutine check_velocities

subroutine invert_z_velocities(atoms,z_low_border,z_high_border)
	type(particles)::	atoms
	integer:: i
	real:: z_low_border,z_high_border

	!$OMP PARALLEL DO private(i)
	do i=1,atoms%N
		if(	(atoms%positions(3,i)>z_low_border .and.&
			atoms%positions(3,i)<(z_low_border+z_high_border)/2 .and.&
			atoms%velocities(3,i)>0.) .or.&
			(atoms%positions(3,i)<z_high_border .and.&
			atoms%positions(3,i)>(z_low_border+z_high_border)/2 .and.&
			atoms%velocities(3,i)<0.)	)	atoms%velocities(3,i) = -atoms%velocities(3,i)
	enddo
	!$OMP END PARALLEL DO
	
end subroutine invert_z_velocities

subroutine position_analysis(av,mi,ma,atoms,group,direction,minimum,maximum)
	type(particles)::	atoms
	type(particle_group):: group
	real:: minimum,maximum,av,mi,ma
	integer::		i,ind,direction,k
	
	k = 0
	av = 0.
	mi = maximum
	ma = minimum
	do ind=1,group%N
		i = group%indexes(ind)
		if (atoms%positions(direction,i)<maximum .and. atoms%positions(direction,i)>minimum) then
			av = av+atoms%positions(direction,i)
			k = k+1
			if (atoms%positions(direction,i)>ma) ma = atoms%positions(direction,i)
			if (atoms%positions(direction,i)<mi) mi = atoms%positions(direction,i)
 		endif
	enddo
	av = av/k
	
end subroutine position_analysis

pure subroutine find_distance(dr,dr2,vec1,vec2,box)
	type(simulation_cell), intent (in) :: box
	real, intent (in) :: vec1(3),vec2(3)
	real, intent (out) :: dr(3),dr2

	dr(1) = vec2(1)-vec1(1)
	dr(2) = vec2(2)-vec1(2)
	dr(3) = vec2(3)-vec1(3)
	!this part consumes a lot of time (approx. 50% of neighbour search and distance)
	!replace with vectorizable function?
	!replace with no pbc function?
	dr(1) = dr(1)-box%half_box_size(1)*(sign(1.,dr(1)-box%half_box_size(1))+sign(1.,dr(1)+box%half_box_size(1)))
	dr(2) = dr(2)-box%half_box_size(2)*(sign(1.,dr(2)-box%half_box_size(2))+sign(1.,dr(2)+box%half_box_size(2)))
	dr(3) = dr(3)-box%half_box_size(3)*(sign(1.,dr(3)-box%half_box_size(3))+sign(1.,dr(3)+box%half_box_size(3)))
	!
	dr2 = dr(1)*dr(1)+dr(2)*dr(2)+dr(3)*dr(3)
	
	return
end subroutine find_distance

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module md_general