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

type neibour_list
	integer:: N,neib_num_max,update_period
	real:: r_cut
	integer,allocatable:: nlist(:,:),nnum(:),lessnnum(:),particle_index(:)
	real,allocatable:: dr(:,:,:),moddr(:,:)
end type neibour_list

type integrator_params
	integer:: l,period_snapshot,period_log
	real:: dt
	character(len=32):: int_name
end type integrator_params

contains!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine read_integrator_params(integr,file_id)
	type(integrator_params):: integr
	integer:: file_id

	read(file_id,*) integr%int_name,integr%dt,integr%l,integr%period_snapshot,integr%period_log

end subroutine read_integrator_params

subroutine init_time_steps(dt,delta_t)
	type(time_steps):: dt
	integer:: i
	real::	delta_t

	do i=1,size(dt%ts)
		dt%ts(i) = delta_t/2**(i-1)
	enddo

	return
end subroutine init_time_steps

subroutine read_box_size(box,filename)
	type(simulation_cell):: box
	real:: boxmatrix(9)
	character(*):: 	filename
	character(len=32)::	c

	open(1,file=filename)
	read(1,*)
	read(1,*) c,boxmatrix
	close(1)
	box%box_size(1)=boxmatrix(1)
	box%box_size(2)=boxmatrix(5)
	box%box_size(3)=boxmatrix(9)
	box%half_box_size = 0.5*box%box_size
	
	return
end subroutine read_box_size

subroutine read_particles(atoms,filename)
	type(particles)::	atoms
	character(*):: 	filename
	integer::		i

	open(1,file=filename)
	read(1,*) atoms%N
	allocate(atoms%positions(3,atoms%N))
	allocate(atoms%velocities(3,atoms%N))
	allocate(atoms%forces(3,atoms%N))
	allocate(atoms%masses(atoms%N))
	allocate(atoms%atom_types(atoms%N))
	read(1,*)
	do i=1,atoms%N
		read(1,*) atoms%positions(:,i),atoms%velocities(:,i),atoms%masses(i),atoms%atom_types(i)
	enddo
	close(1)

	return
end subroutine read_particles

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

subroutine write_particle_group(filename,atoms,group,box)
	type(particles)::	atoms
	type(particle_group)::	group
	type(simulation_cell):: box
	character(*):: 	filename
	integer::		i,ind

	open(2,file=filename)
	write(2,*) group%N
	write(2,'(A,9f16.6,A)') 'Lattice="',box%box_size(1), 0., 0., 0., box%box_size(2), 0., 0., 0., box%box_size(3),&
			' " Properties=pos:R:3:vel:R:3:mass:R:1:species:S:1'
	do ind=1,group%N
		i = group%indexes(ind)
		write(2,'(7f27.16,A,A)') atoms%positions(:,i),atoms%velocities(:,i),atoms%masses(i),'    ',atoms%atom_types(i)
	enddo
	close(2)

	return
end subroutine write_particle_group

subroutine integrate_verlet_positions(atoms,group,steps,box)
	type(particles)::	atoms
	type(particle_group):: group
	type(time_steps):: steps
	type(simulation_cell):: box
	integer::		i,k,ind

	!$OMP PARALLEL firstprivate(i,ind,k)
	!$OMP DO
	do ind=1,group%N
		i = group%indexes(ind)
		do k=1,3
			atoms%positions(k,i) = atoms%positions(k,i)+atoms%velocities(k,i)*steps%ts(1)
			if (atoms%positions(k,i)>box%box_size(k)) then
				atoms%positions(k,i) = atoms%positions(k,i)-box%box_size(k)
			elseif (atoms%positions(k,i)<0.) then
				atoms%positions(k,i) = atoms%positions(k,i)+box%box_size(k)
			endif
		enddo
	enddo
	!$OMP END DO
	!$OMP END PARALLEL

	return
end subroutine integrate_verlet_positions

subroutine integrate_verlet_velocities(atoms,group,steps)
	type(particles)::	atoms
	type(particle_group):: group
	type(time_steps):: steps
	real,parameter:: mass_coef=1.6605389217/1.6021765654*10.**(2)
	integer::		i,ind,k

	!$OMP PARALLEL firstprivate(i,ind,k)
	!$OMP DO
	do ind=1,group%N
		i = group%indexes(ind)
		do k=1,3
			atoms%velocities(k,i) = atoms%velocities(k,i)+atoms%forces(k,i)/atoms%masses(i)/mass_coef*steps%ts(2)
		enddo
	enddo
	!$OMP END DO
	!$OMP END PARALLEL

	return
end subroutine integrate_verlet_velocities

subroutine molecular_static_velocities(atoms,group)
	type(particles)::	atoms
	type(particle_group):: group
	real:: fv,ff
	integer::		i,ind

	!$OMP PARALLEL firstprivate(i,ind,fv,ff)
	!$OMP DO
	do ind=1,group%N
		i = group%indexes(ind)
		fv = sum(atoms%forces(:,i)*atoms%velocities(:,i))
		ff = sum(atoms%forces(:,i)**2)
		if (fv>0. .and. ff>10.**(-12) ) then
			atoms%velocities(:,i) = fv/ff*atoms%forces(:,i)
		else
			atoms%velocities(:,i) = 0.
		endif
	enddo
	!$OMP END DO
	!$OMP END PARALLEL

	return
end subroutine molecular_static_velocities

subroutine zero_forces(atoms,group)
	type(particles)::	atoms
	type(particle_group):: group
	integer:: i,ind,k

	!$OMP PARALLEL firstprivate(i,ind,k)
	!$OMP DO
	do ind=1,group%N
		i = group%indexes(ind)
		do k=1,3
			atoms%forces(k,i) = 0.
		enddo
	enddo
	!$OMP END DO
	!$OMP END PARALLEL
	
end subroutine zero_forces

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

subroutine random_velocities(atoms,group)
	type(particles)::	atoms
	type(particle_group):: group
	real a1,a2,b
	integer::		i,k,ind
	integer omp_get_thread_num
	
	!$OMP PARALLEL firstprivate(i,ind,k,a1,a2,b)
	call srand(omp_get_thread_num())
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

subroutine random_momenta(atoms,group)
	type(particles)::	atoms
	type(particle_group):: group
	real,parameter:: coef = 1.3806488/1.6605389217*10.**(-6)
	integer::		i,ind

	call random_velocities(atoms,group)

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

subroutine set_new_temperature(atoms,group,temp)
	type(particles)::	atoms
	type(particle_group):: group
	real:: temp,ke,temperature
	
	call random_momenta(atoms,group)
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
	integer:: k

	do k=1,3
		dr(k) = vec2(k)-vec1(k)
		if (dr(k)>box%half_box_size(k)) then
			dr(k) = dr(k)-box%box_size(k)
		elseif (dr(k)<-box%half_box_size(k)) then
			dr(k) = dr(k)+box%box_size(k)
		endif
	enddo
	dr2 = dr(1)*dr(1)+dr(2)*dr(2)+dr(3)*dr(3)
	
	return
end subroutine find_distance

subroutine create_nose_hoover_chain(nhc)
	type(nose_hoover_chain):: nhc

	allocate(nhc%x(nhc%M))
	allocate(nhc%v(nhc%M))
	allocate(nhc%q(nhc%M))
	nhc%x = 0.
	nhc%v = 0.
	nhc%q = 0.
	nhc%e = 0.
	nhc%s = 1.

	return
end subroutine create_nose_hoover_chain

subroutine set_nose_hoover_chain(nhc,temp,q1,l)
	type(nose_hoover_chain):: nhc
	integer:: l,i
	real::	q1,temp

	nhc%L = l
	nhc%temperature = temp
	nhc%q(1) = q1
	do i=2,nhc%M,1
		nhc%q(i) = nhc%q(1)/(3.d0*nhc%L)
	enddo

	return
end subroutine set_nose_hoover_chain

subroutine integrate_nose_hoover_chain(nhc,atoms,group,dt)
	type(particles)::	atoms
	type(particle_group):: group
	type(time_steps):: dt
	type(nose_hoover_chain):: nhc
	real:: kedif,ke,kt,b
	integer::		i

	call calculate_kinetic_energy(ke,atoms,group)
	kt=1.3806488/1.6021765654*10.**(-4)*nhc%temperature
	kedif = 2.*ke-3.*nhc%L*kt

	if (nhc%M==1) then
		nhc%v(1)=nhc%v(1)+kedif/nhc%q(1)*dt%ts(3)
	else
		nhc%v(nhc%M)=nhc%v(nhc%M)+(nhc%q(nhc%M-1)*nhc%v(nhc%M-1)**2-kt)/nhc%q(nhc%M)*dt%ts(3)
		do i=nhc%M-1,2,-1
			b=exp(-nhc%v(i+1)*dt%ts(4))
			nhc%v(i)=nhc%v(i)*b**2+(nhc%q(i-1)*nhc%v(i-1)**2-kt)/nhc%q(i)*dt%ts(3)*b
		enddo
		b=exp(-nhc%v(2)*dt%ts(4))
		nhc%v(1)=nhc%v(1)*b**2+kedif/nhc%q(1)*dt%ts(3)*b
	endif

	nhc%s=exp(-nhc%v(1)*dt%ts(2))
	call scale_velocities(atoms,group,nhc%s)
	kedif = 2.*ke*nhc%s**2-3.*nhc%L*kt
	do i=1,nhc%M,1
		nhc%x(i)=nhc%x(i)+nhc%v(i)*dt%ts(2)
	enddo

	if (nhc%M==1) then
		nhc%v(1)=nhc%v(1)+kedif/nhc%q(1)*dt%ts(3)
	else
		nhc%v(1)=nhc%v(1)*b**2+kedif/nhc%q(1)*dt%ts(3)*b
		do i=2,nhc%M-1,1
			b=exp(-nhc%v(i+1)*dt%ts(4))
			nhc%v(i)=nhc%v(i)*b**2+(nhc%q(i-1)*nhc%v(i-1)**2-kt)/nhc%q(i)*dt%ts(3)*b
		enddo
		nhc%v(nhc%M)=nhc%v(nhc%M)+(nhc%q(nhc%M-1)*nhc%v(nhc%M-1)**2-kt)/nhc%q(nhc%M)*dt%ts(3)
	endif

	return
end subroutine integrate_nose_hoover_chain

subroutine calculate_nose_hoover_chain_energy(nhc)
	type(nose_hoover_chain):: nhc
	real:: kt
	integer::		i

	kt=1.3806488/1.6021765654*10.**(-4)*nhc%temperature

	nhc%e = nhc%q(1)/2*nhc%v(1)**2+3.*nhc%L*kt*nhc%x(1)
	do i=2,nhc%M,1
		nhc%e = nhc%e+nhc%q(i)/2*nhc%v(i)**2+kt*nhc%x(i)
	enddo

	return
end subroutine calculate_nose_hoover_chain_energy

subroutine create_neibour_list(nl)
	type(neibour_list):: nl

	if (nl%N>0) then
		allocate(nl%particle_index(nl%N))
		allocate(nl%nnum(nl%N))
		allocate(nl%lessnnum(nl%N))
		allocate(nl%nlist(nl%neib_num_max,nl%N))
		allocate(nl%dr(3,nl%neib_num_max,nl%N))
		allocate(nl%moddr(nl%neib_num_max,nl%N))
		nl%particle_index = 0
		nl%nnum = 0
		nl%lessnnum = 0
		nl%nlist = 0
		nl%dr = 0.
		nl%moddr = 0.
	endif
	
	return
end subroutine create_neibour_list

subroutine update_neibour_list(md_step,nl,atoms,group1,group2,box)
	type(particles)::	atoms
	type(particle_group):: group1,group2
	type(simulation_cell):: box
	type(neibour_list):: nl
	integer:: md_step

	if (group1%N>0 .and. group2%N>0 .and. nl%N>0) then
		if (mod(md_step,nl%update_period)==0) then
			call find_neibours(nl,atoms,group1,group2,box)
		else
			call find_neibour_distances(nl,atoms,group1,group2,box)
		endif
	endif
	
end subroutine update_neibour_list

subroutine find_neibours(nl,atoms,group1,group2,box)
	type(particles)::	atoms
	type(particle_group):: group1,group2
	type(simulation_cell):: box
	type(neibour_list):: nl
	real:: dr(3),dr2
	integer:: i,j,k,ind,jnd,nnumind,lessnnumind

	if (group1%N/=nl%N) then; write(*,*) 'error: group%N/=nl%N',group1%N,nl%N; stop; endif
	!$OMP PARALLEL firstprivate(i,j,k,ind,jnd,dr,dr2,nnumind,lessnnumind)
	!$OMP DO
	do ind=1,group1%N
		i = group1%indexes(ind)
		nl%particle_index(ind) = group1%indexes(ind)
		!nl%nlist(:nl%nnum(ind),ind) = 0
		nnumind = 0
		lessnnumind = -1
		do jnd=1,group2%N
			j = group2%indexes(jnd)
			if (i/=j) then
				call find_distance(dr,dr2,atoms%positions(:,i),atoms%positions(:,j),box)
				if (dr2<nl%r_cut*nl%r_cut) then
					if (lessnnumind==-1 .and. i<j) lessnnumind = nnumind
					nnumind = nnumind+1
					if (nnumind>nl%neib_num_max) then; write(*,*) 'error: too many neibours',ind,nnumind; stop; endif
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
end subroutine find_neibours

subroutine find_neibour_distances(nl,atoms,group1,group2,box)
	type(particles)::	atoms
	type(particle_group):: group1,group2
	type(simulation_cell):: box
	type(neibour_list):: nl
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
end subroutine find_neibour_distances

subroutine converce_neigbour_list(cnl,group,nl)
	type(neibour_list):: nl,cnl
	type(particle_group):: group
	integer::		i,j,p

	if (group%N>0 .and. nl%N>0 .and. cnl%N>0) then
		if (group%N/=cnl%N) then; write(*,*) 'error: group%N/=cnl%N'; stop; endif
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
				if (cnl%nnum(j)>cnl%neib_num_max) then; write(*,*) 'error: too many neibours',j,cnl%nnum(j); stop; endif
				cnl%nlist(cnl%nnum(j),j) = i
				cnl%dr(:,cnl%nnum(j),j) = -nl%dr(:,p,i)
				cnl%moddr(cnl%nnum(j),j) = nl%moddr(p,i)
			enddo
		enddo
	endif
	
	return
end subroutine converce_neigbour_list

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module md_general